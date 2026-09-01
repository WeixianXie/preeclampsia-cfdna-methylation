# -*- coding: utf-8 -*-
"""56b: Positive control for the DMR-gene validation pipeline.

Question: is the weak DMR-gene performance (AUC ~0.52-0.62 in GSE85307) due to
the pipeline or to the gene set? Train the SAME elastic-net pipeline on the
top-MAD genome-wide genes of GSE85307 (excluding DMR genes) as positive
control, and on matched random gene sets of equal size as negative controls.
"""
import gzip, os, sys, warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore')
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def read_series(path):
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
                    meta.setdefault(vals[0].split(':', 1)[0].strip(), vals)
            elif s.startswith('!series_matrix_table_begin'):
                hdr = [c.strip('"') for c in fh.readline().rstrip('\n').split('\t')]
                in_tab = True
            elif s.startswith('!series_matrix_table_end'):
                in_tab = False
            elif in_tab and s:
                rows.append([c.strip('"') for c in s.split('\t')])
    df = pd.DataFrame(rows, columns=hdr).set_index(hdr[0])
    df = df.apply(pd.to_numeric, errors='coerce')
    df.columns = samples
    return samples, meta, df


def read_annot(path):
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


def main():
    gs, meta, ex = read_series(os.path.join(
        ROOT, 'data', 'geo_transcriptome', 'GSE85307-GPL6244_series_matrix.txt.gz'))
    cond = [v.split(':', 1)[1].strip() for v in meta['pregnancy_condition']]
    y = np.array([1 if 'Preeclampsia' in c else 0 for c in cond])

    an = read_annot(os.path.join(ROOT, 'data', 'annot', 'GPL6244.annot.gz'))
    vals = ex.values
    med = np.nanmedian(vals, axis=1)
    mad = pd.Series(np.nanmedian(np.abs(vals - med[:, None]), axis=1),
                    index=ex.index)  # per-probe MAD over samples
    sym = ex.index.map(lambda p: an.get(p, np.nan))
    ok = sym.notna() & (mad > mad.quantile(0.2))
    probes = ex.index[ok]
    symbols = sym[ok]

    # collapse probes -> gene (max-MAD probe per gene)
    gg = pd.DataFrame({'sym': symbols.values, 'mad': mad[probes].values},
                      index=probes)
    top_probe = gg.sort_values('mad', ascending=False).groupby('sym').head(1)
    X = ex.loc[top_probe.index].T
    X.columns = top_probe['sym'].values

    # DMR genes to exclude
    import csv
    dmr = set()
    with open(os.path.join(ROOT, 'results', 'GSE282512_dmr_final_v2.csv'),
              encoding='utf-8-sig') as fh:
        for r in csv.DictReader(fh):
            s = (r.get('symbol') or '').strip().strip('"').upper()
            if s:
                dmr.add(s)
    X = X.loc[:, [c for c in X.columns if c not in dmr]]

    Xv_med = np.nanmedian(X.values, axis=0)
    Xmad = np.nanmedian(np.abs(X.values - Xv_med[None, :]), axis=0)
    order = X.columns[np.argsort(-Xmad)]
    Xv = X[order[:500]].values
    keep = ~np.isnan(Xv).any(axis=1)
    Xv, yk = Xv[keep], y[keep]
    mu, sd = Xv.mean(0), Xv.std(0, ddof=1)
    sd[sd == 0] = 1
    Z = (Xv - mu) / sd
    print('positive-control matrix:', Z.shape)

    skf = StratifiedKFold(5, shuffle=True, random_state=42)
    oof = np.zeros(len(yk))
    for tr, te in skf.split(Z, yk):
        lr = LogisticRegression(penalty='elasticnet', solver='saga', C=0.1,
                                l1_ratio=0.5, max_iter=8000)
        lr.fit(Z[tr], yk[tr])
        oof[te] = lr.predict_proba(Z[te])[:, 1]
    auc_pos = roc_auc_score(yk, oof)
    print(f'POSITIVE CONTROL genome-wide top-500 OOF AUC = {auc_pos:.3f}')

    # matched random gene sets (size 58, excl. DMR) x 50
    rng = np.random.default_rng(7)
    pool = X.columns.tolist()
    aucs = []
    for b in range(50):
        pick = rng.choice(pool, 58, replace=False)
        Xr = X[pick].values
        keep = ~np.isnan(Xr).any(axis=1)
        Xr, yrr = Xr[keep], y[keep]
        mu, sd = Xr.mean(0), Xr.std(0, ddof=1)
        sd[sd == 0] = 1
        Zr = (Xr - mu) / sd
        oofr = np.zeros(len(yrr))
        for tr, te in skf.split(Zr, yrr):
            lr = LogisticRegression(penalty='elasticnet', solver='saga', C=14.7,
                                    l1_ratio=0.5, max_iter=8000)
            lr.fit(Zr[tr], yrr[tr])
            oofr[te] = lr.predict_proba(Zr[te])[:, 1]
        aucs.append(roc_auc_score(yrr, oofr))
    aucs = np.array(aucs)
    print(f'RANDOM-58 null: mean AUC = {aucs.mean():.3f} '
          f'[{np.percentile(aucs, 2.5):.3f},{np.percentile(aucs, 97.5):.3f}]')

    with open(os.path.join(ROOT, 'results',
              'GSE85307_GSE86200_gene_validation_summary.txt'), 'a',
              encoding='utf-8') as fh:
        fh.write('\nPositive / null controls (scripts/56b, GSE85307, 5-fold OOF):\n')
        fh.write(f'  genome-wide top-500 MAD genes: AUC = {auc_pos:.3f}\n')
        fh.write(f'  50x random 58-gene sets (excl. DMR): mean AUC = {aucs.mean():.3f} '
                 f'[{np.percentile(aucs, 2.5):.3f}, {np.percentile(aucs, 97.5):.3f}]\n')
    print('DONE')


if __name__ == '__main__':
    main()
