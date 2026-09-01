# -*- coding: utf-8 -*-
"""56: Early-pregnancy DMR-gene model external validation.

Discovery/train:  GSE85307 (GPL6244, HuGene 1.0 ST, PAXgene whole blood,
                  10-18 wk gestation, 47 PE vs 110 CT)
External test:    GSE86200 (GPL10558, Illumina HT-12 V4, VDAART,
                  enrollment [early pregnancy] 6 PE vs 24 CT,
                  and 32-38 wk 6 PE vs 24 CT)

Gene set: unique gene symbols of the 166 cfDNA DMRs
          (results/GSE282512_dmr_final_v2.csv, column `symbol`).

Models:
  A) Elastic-net logistic regression (5-fold stratified CV for C selection,
     out-of-fold probabilities for internal AUC), applied unchanged to GSE86200.
  B) Unsupervised signature score: per-gene z (train stats) * sign(train
     Welch t) summed over mapped genes (sensitivity check).

Bootstrap (stratified, 1000x) 95% CI for every AUC.
Outputs: results/GSE85307_GSE86200_gene_validation_{auc,genes}.csv,
         results/GSE85307_GSE86200_gene_validation_summary.txt,
         figures/GSE85307_GSE86200_gene_roc.png
"""
import gzip, csv, os, sys, warnings
import numpy as np
import pandas as pd
from collections import Counter
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VPY = sys.executable

RNG = np.random.default_rng(42)
NBOOT = 1000


def get_meta(meta, target):
    """Find a characteristics key ignoring space/underscore differences."""
    tgt = target.replace(' ', '').replace('_', '').lower()
    for k, v in meta.items():
        if k.replace(' ', '').replace('_', '').lower() == tgt:
            return [x.split(':', 1)[1].strip() if ':' in x else x for x in v]
    raise KeyError(target)


# ---------------------------------------------------------------- IO helpers
def read_series(path):
    """Return (samples, meta{key: vals}, expr DataFrame)."""
    meta, samples = {}, None
    rows, hdr, in_tab = [], None, False
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            s = line.rstrip('\n')
            if s.startswith('!Sample_geo_accession'):
                samples = [x.strip('"') for x in s.split('\t')[1:]]
            elif s.startswith('!Sample_characteristics_ch1'):
                vals = [x.strip('"') for x in s.split('\t')[1:]]
                if vals and ':' in vals[0]:
                    key = vals[0].split(':', 1)[0].strip()
                    meta.setdefault(key, vals)
            elif s.startswith('!Sample_source_name_ch1'):
                meta.setdefault('source', [x.strip('"') for x in s.split('\t')[1:]])
            elif s.startswith('!series_matrix_table_begin'):
                hdr = [c.strip('"') for c in fh.readline().rstrip('\n').split('\t')]
                in_tab = True
            elif s.startswith('!series_matrix_table_end'):
                in_tab = False
            elif in_tab and s:
                rows.append([c.strip('"') for c in s.split('\t')])
    df = pd.DataFrame(rows, columns=hdr).set_index(hdr[0])
    df = df.apply(pd.to_numeric, errors='coerce')
    assert samples is not None and len(samples) == df.shape[1]
    df.columns = samples
    return samples, meta, df


def read_annot(path):
    """GEO .annot.gz -> {probe_id: GENE SYMBOL (upper)}."""
    m = {}
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
        cols = None
        for line in fh:
            s = line.rstrip('\n')
            if s.startswith('!platform_table_begin'):
                cols = fh.readline().rstrip('\n').split('\t')
                idx, sym = cols.index('ID'), cols.index('Gene symbol')
                continue
            if s.startswith('!platform_table_end'):
                break
            if cols is not None and s and not s.startswith('#'):
                p = s.split('\t')
                if len(p) > sym and p[sym].strip():
                    m[p[idx].strip('"')] = p[sym].strip('"').upper()
    return m


def dmr_genes():
    genes = []
    with open(os.path.join(ROOT, 'results', 'GSE282512_dmr_final_v2.csv'),
              encoding='utf-8-sig') as fh:
        for r in csv.DictReader(fh):
            sym = (r.get('symbol') or '').strip().strip('"').upper()
            if sym:
                genes.append(sym)
    return sorted(set(genes))


def gene_matrix(expr, probe2sym, gene_set):
    """Collapse probes -> gene symbols (mean of matched probes)."""
    hit = {p: s for p, s in probe2sym.items() if s in gene_set and p in expr.index}
    if not hit:
        return pd.DataFrame(columns=sorted(gene_set))
    sub = expr.loc[list(hit)].rename(index=hit)
    g = sub.groupby(level=0).mean()
    g = g.reindex(columns=expr.columns)
    return g


# ---------------------------------------------------------------- statistics
def boot_auc(y, p, nboot=NBOOT, seed=0):
    y = np.asarray(y); p = np.asarray(p)
    rng = np.random.default_rng(seed)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    aucs = []
    for _ in range(nboot):
        ii = np.concatenate([rng.choice(pos, len(pos), True),
                             rng.choice(neg, len(neg), True)])
        yy, pp = y[ii], p[ii]
        if yy.min() == yy.max():
            continue
        aucs.append(roc_auc_score(yy, pp))
    if len(aucs) < 200:
        return np.nan, np.nan, np.nan
    return np.percentile(aucs, 2.5), np.percentile(aucs, 97.5), np.mean(aucs)


def main():
    out = []
    genes = dmr_genes()
    print(f'DMR gene set: {len(genes)} unique symbols')

    # ---------- load cohorts
    gs85307, meta85307, ex85 = read_series(os.path.join(
        ROOT, 'data', 'geo_transcriptome', 'GSE85307-GPL6244_series_matrix.txt.gz'))
    cond85 = get_meta(meta85307, 'pregnancy condition')
    y85 = np.array([1 if 'Preeclampsia' in c else 0 for c in cond85])
    print(f'GSE85307: n={len(y85)}, PE={y85.sum()}, CT={(1 - y85).sum()}')

    gs86200, meta86200, ex86 = read_series(os.path.join(
        ROOT, 'data', 'geo_transcriptome', 'GSE86200_series_matrix.txt.gz'))
    cond86 = get_meta(meta86200, 'pregnancy condition')
    src86 = meta86200['source']
    y86 = np.array([1 if 'Preeclampsia' in c else 0 for c in cond86])
    early = np.array(['enrollment' in s for s in src86])
    late = np.array(['32-38 weeks' in s for s in src86])
    print(f'GSE86200: n={len(y86)}, PE={y86.sum()}; early={early.sum()}, late={late.sum()}')

    # ---------- annotations
    an85 = read_annot(os.path.join(ROOT, 'data', 'annot', 'GPL6244.annot.gz'))
    an86 = read_annot(os.path.join(ROOT, 'data', 'annot', 'GPL10558.annot.gz'))
    cov85 = len(set(genes) & set(an85.values()))
    cov86 = len(set(genes) & set(an86.values()))
    print(f'gene coverage: GSE85307 {cov85}/{len(genes)}, GSE86200 {cov86}/{len(genes)}')

    G85 = gene_matrix(ex85, an85, set(genes))
    G86 = gene_matrix(ex86, an86, set(genes))
    shared = sorted(set(G85.index) & set(G86.index))
    print(f'genes with data in BOTH cohorts: {len(shared)}')

    # ---------- train on GSE85307 (shared genes only for transferability)
    use = shared if len(shared) >= 10 else sorted(G85.index)
    X85 = G85.loc[use].T.values            # samples x genes
    keep = ~np.isnan(X85).any(axis=1)
    X85, y85k = X85[keep], y85[keep]
    mu, sd = X85.mean(0), X85.std(0, ddof=1)
    sd[sd == 0] = 1.0
    Z85 = (X85 - mu) / sd
    print(f'train matrix: {Z85.shape[0]} samples x {Z85.shape[1]} genes')

    # Model B direction (Welch t sign)
    from scipy.stats import ttest_ind
    tstat, _ = ttest_ind(Z85[y85k == 1], Z85[y85k == 0], axis=0, equal_var=False)
    tsign = np.sign(tstat)

    def sig_score(Z):
        return (Z * tsign).mean(1)

    # ---------- Model A: elastic-net logistic + 5-fold OOF
    skf = StratifiedKFold(5, shuffle=True, random_state=42)
    Cs = np.logspace(-2, 2, 25)
    best_c, best_auc = None, -1
    for C in Cs:
        aucs = []
        for tr, te in skf.split(Z85, y85k):
            lr = LogisticRegression(penalty='elasticnet', solver='saga', C=C,
                                    l1_ratio=0.5, max_iter=6000)
            lr.fit(Z85[tr], y85k[tr])
            aucs.append(roc_auc_score(y85k[te], lr.predict_proba(Z85[te])[:, 1]))
        m = float(np.mean(aucs))
        if m > best_auc:
            best_auc, best_c = m, C
    print(f'elastic-net best C={best_c:.3g} inner mean AUC={best_auc:.3f}')

    oof = np.zeros(len(y85k))
    for tr, te in skf.split(Z85, y85k):
        lr = LogisticRegression(penalty='elasticnet', solver='saga', C=best_c,
                                l1_ratio=0.5, max_iter=6000)
        lr.fit(Z85[tr], y85k[tr])
        oof[te] = lr.predict_proba(Z85[te])[:, 1]
    auc85_enet = roc_auc_score(y85k, oof)
    lo, hi, _ = boot_auc(y85k, oof)
    print(f'GSE85307 elastic-net OOF AUC = {auc85_enet:.3f} [{lo:.3f},{hi:.3f}]')
    out.append(('GSE85307 elastic-net OOF (train cohort)', auc85_enet, lo, hi,
                len(y85k), int(y85k.sum())))

    auc85_sig = roc_auc_score(y85k, sig_score(Z85))
    lo, hi, _ = boot_auc(y85k, sig_score(Z85))
    print(f'GSE85307 signature-score AUC = {auc85_sig:.3f} [{lo:.3f},{hi:.3f}]')
    out.append(('GSE85307 signature-score (train cohort)', auc85_sig, lo, hi,
                len(y85k), int(y85k.sum())))

    # final full-train models
    lr_full = LogisticRegression(penalty='elasticnet', solver='saga', C=best_c,
                                 l1_ratio=0.5, max_iter=6000).fit(Z85, y85k)
    coef = lr_full.coef_.ravel()

    # ---------- external validation on GSE86200
    X86all = G86.loc[use].T.values
    col86 = np.array(G86.columns)
    res_roc = {}
    for name, mask in [('GSE86200 enrollment (early pregnancy)', early),
                       ('GSE86200 32-38 wk (late pregnancy)', late)]:
        Xs = X86all[mask]
        ys = y86[mask]
        cols = col86[mask]
        ok = ~np.isnan(Xs).any(axis=1)
        Xs, ys, cols = Xs[ok], ys[ok], cols[ok]
        if ys.min() == ys.max():
            print(f'{name}: single class, skip'); continue
        Zs = (Xs - mu) / sd
        p_enet = lr_full.predict_proba(Zs)[:, 1]
        p_sig = sig_score(Zs)
        for tag, p in [('elastic-net', p_enet), ('signature-score', p_sig)]:
            a = roc_auc_score(ys, p)
            lo, hi, _ = boot_auc(ys, p, seed=1)
            out.append((f'{name} {tag}', a, lo, hi, len(ys), int(ys.sum())))
            print(f'{name} {tag}: AUC={a:.3f} [{lo:.3f},{hi:.3f}] (n={len(ys)}, PE={ys.sum()})')
        res_roc[name] = (ys, p_enet, p_sig)

    # ---------- single-gene screen (external, enrollment) for report
    sing = []
    if 'GSE86200 enrollment (early pregnancy)' in res_roc:
        ys, _, _ = res_roc['GSE86200 enrollment (early pregnancy)']
        mask_ok = early.copy()
        ok = ~np.isnan(X86all[early]).any(axis=1)
        idx_ok = np.where(early)[0][ok]
        for j, g in enumerate(use):
            xg = X86all[idx_ok, j]
            if np.nanstd(xg) == 0:
                continue
            a = roc_auc_score(ys, xg)
            sing.append((g, a, 1.0 if a >= 0.5 else -1.0))
    sing_df = pd.DataFrame(sing, columns=['gene', 'auc_single_gse86200_early', 'side'])
    sing_df = sing_df.sort_values('auc_single_gse86200_early', ascending=False)

    # ---------- save
    pd.DataFrame([{'setting': s, 'auc': round(a, 4),
                   'ci_low': round(lo, 4), 'ci_high': round(hi, 4),
                   'n': n, 'n_pe': np_} for s, a, lo, hi, n, np_ in out]
                 ).to_csv(os.path.join(ROOT, 'results',
                 'GSE85307_GSE86200_gene_validation_auc.csv'), index=False)

    pd.DataFrame({'gene': use, 'coef_enet': coef,
                  't_sign_train': tsign,
                  'train_mean_pe': Z85[y85k == 1].mean(0),
                  'train_mean_ct': Z85[y85k == 0].mean(0)}
                 ).merge(sing_df, how='left', left_on='gene', right_on='gene'
                 ).to_csv(os.path.join(ROOT, 'results',
                 'GSE85307_GSE86200_gene_validation_genes.csv'), index=False)

    # ---------- figure: ROC curves
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3))
    panels = [('GSE85307 train (OOF, elastic-net)', y85k, oof),
              ('GSE86200 early pregnancy (external)', None, None),
              ('GSE86200 32-38 wk (external)', None, None)]
    for ax, (title, yy, pp) in zip(axes, panels):
        if yy is None:
            key = ('GSE86200 enrollment (early pregnancy)' if 'early' in title
                   else 'GSE86200 32-38 wk (late pregnancy)')
            if key in res_roc:
                yy, pp, p2 = res_roc[key]
                for lab, pv, col in [('elastic-net', pp, '#C0392B'),
                                     ('signature-score', p2, '#2471A3')]:
                    fpr, tpr, _ = roc_curve(yy, pv)
                    a = roc_auc_score(yy, pv)
                    ax.plot(fpr, tpr, color=col, lw=2,
                            label=f'{lab}: AUC={a:.2f}')
        else:
            fpr, tpr, _ = roc_curve(yy, pp)
            ax.plot(fpr, tpr, color='#C0392B', lw=2,
                    label=f'elastic-net: AUC={roc_auc_score(yy, pp):.2f}')
            f2, t2, _ = roc_curve(yy, sig_score(Z85))
            ax.plot(f2, t2, color='#2471A3', lw=2, ls='--',
                    label=f'signature: AUC={roc_auc_score(yy, sig_score(Z85)):.2f}')
        ax.plot([0, 1], [0, 1], color='grey', ls=':', lw=1)
        ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
        ax.set_title(title, fontsize=11)
        ax.legend(loc='lower right', fontsize=9)
        ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    fig.suptitle('166-DMR gene model: early-pregnancy external validation '
                 '(GSE85307 train / GSE86200 external)', fontsize=12)
    fig.tight_layout()
    figp = os.path.join(ROOT, 'figures', 'GSE85307_GSE86200_gene_roc.png')
    fig.savefig(figp, dpi=200); plt.close(fig)
    print('figure saved:', figp)

    # ---------- summary txt
    with open(os.path.join(ROOT, 'results',
              'GSE85307_GSE86200_gene_validation_summary.txt'),
              'w', encoding='utf-8') as fh:
        fh.write('Early-pregnancy gene-model external validation (scripts/56)\n')
        fh.write('=' * 72 + '\n')
        fh.write(f'DMR gene set: {len(genes)} unique symbols from 166 DMRs\n')
        fh.write(f'Gene coverage: GSE85307 {cov85}, GSE86200 {cov86}, '
                 f'both {len(shared)} (genes used in model: {len(use)})\n\n')
        for s, a, lo, hi, n, np_ in out:
            fh.write(f'{s:55s} AUC={a:.3f} [{lo:.3f},{hi:.3f}] n={n} PE={np_}\n')
        fh.write('\nTop single genes in GSE86200 early pregnancy:\n')
        for _, r in sing_df.head(10).iterrows():
            fh.write(f'  {r.gene:15s} AUC={r.auc_single_gse86200_early:.3f}\n')
    print('DONE')


if __name__ == '__main__':
    main()
