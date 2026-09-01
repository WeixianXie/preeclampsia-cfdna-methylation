# -*- coding: utf-8 -*-
"""Compose Figures 2-8 for submission: multi-panel layout, 300 dpi TIFF (LZW) + PNG.

Row-based composition: each figure is a stack of rows; images in a row share the
same pixel height; every row is scaled so the widest row matches the target
pixel width. Panel letters (bold) drawn top-left of each panel.
"""
import os
from PIL import Image, ImageDraw, ImageFont

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")

DPI = 300
FULL_W_MM = 183.0          # double-column width (Clinical Epigenetics / BMC)
FULL_W_PX = int(round(FULL_W_MM / 25.4 * DPI))   # ~2161 px
PAD = 24                   # white padding between panels (px)
LABEL_SZ = 96              # panel letter font size (px on 300dpi canvas)

# matplotlib bundled font (guaranteed present with matplotlib)
import matplotlib
FONT_BOLD = os.path.join(matplotlib.get_data_path(), "fonts", "ttf", "DejaVuSans-Bold.ttf")


def load(name):
    return Image.open(os.path.join(FIGDIR, name + ".png")).convert("RGB")


def compose(figure_name, rows, target_w=None, labels=None, label_color=(0, 0, 0)):
    """rows: list of list of source stems (no extension).
    labels: matching nested list of panel letters (or None per cell)."""
    target_w = target_w or FULL_W_PX
    imgs = [[load(stem) for stem in row] for row in rows]
    # per-row height so that each row, at equal-height panels, fits target_w
    row_heights = []
    for row in imgs:
        # width of row if height = h: sum(w_i/h * h)= sum(w_i)?? per-image: w_i / h_i
        # scale factor s applied to all: row width = sum(w_i * s), keep aspect per image:
        # choose s such that sum(w_i)*s + PAD*(n-1) = target_w  -> but images keep own aspect,
        # equal height h means row width = h * sum(r_i) where r_i = w_i/h_i
        ratios = [im.size[0] / im.size[1] for im in row]
        n = len(row)
        h = (target_w - PAD * (n - 1)) / sum(ratios)
        row_heights.append(h)
    H = sum(int(round(h)) for h in row_heights) + PAD * (len(rows) - 1) + 2 * PAD
    W = target_w + 2 * PAD
    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype(FONT_BOLD, LABEL_SZ)
    y = PAD
    for ri, row in enumerate(imgs):
        h = int(round(row_heights[ri]))
        x = PAD
        for ci, im in enumerate(row):
            w = int(round(im.size[0] * h / im.size[1]))
            im_s = im.resize((w, h), Image.LANCZOS)
            canvas.paste(im_s, (x, y))
            lab = labels[ri][ci] if labels else None
            if lab:
                # letter just above-left of the panel, on white margin
                bbox = draw.textbbox((0, 0), lab, font=font)
                ty = y - (bbox[3] - bbox[1]) - 6
                if ty < 2:
                    ty = y + 8
                draw.text((x + 4, ty), lab, font=font, fill=label_color)
            x += w + PAD
        y += h + PAD
    for ext, kw in (("tif", dict(compression="tiff_lzw", dpi=(DPI, DPI))),
                    ("png", dict(dpi=(DPI, DPI)))):
        out = os.path.join(FIGDIR, "%s.%s" % (figure_name, ext))
        canvas.save(out, **kw)
        print("saved", out, canvas.size)


if __name__ == "__main__":
    # Figure 2: DMR discovery and robustness (A-D)
    compose("Figure2_dmr_discovery",
            rows=[["dmr_volcano", "dmr_manhattan"],
                  ["dmr_sensitivity_scatter", "dmr_site_top6"]],
            labels=[["A", "B"], ["C", "D"]])

    # Figure 3: Classifier evaluation (A-D)
    compose("Figure3_classifier",
            rows=[["roc_loocv", "calibration"],
                  ["dca", "shap_importance"]],
            labels=[["A", "B"], ["C", "D"]])

    # Figure 4: Placental origin testing (single wide panel, no letters needed)
    compose("Figure4_placental_origin",
            rows=[["GSE282512_placenta_direct"]])

    # Figure 5: Longitudinal trajectories (A-C; C on its own row, capped width)
    compose("Figure5_trajectories",
            rows=[["GSE154378_typeI_trajectory", "GSE154378_dmr_trajectory"],
                  ["GSE37722_GSE154378_cell_fraction"]],
            labels=[["A", "B"], ["C"]])

    # Figure 6: Subgroup change-based deconvolution (A-B)
    compose("Figure6_subgroup_deconv",
            rows=[["GSE154378_subgroup_deltaf_forest", "GSE154378_subgroup_deltaf_trajectory"]],
            labels=[["A", "B"]])

    # Figure 7: MR for PE/GH - tall forest, cap width at 120 mm to stay under page height
    compose("Figure7_mr_forest",
            rows=[["GSE282512_mr_forest"]],
            target_w=int(round(120 / 25.4 * DPI)))

    # Figure 8: Organ-wide MR - cap width at 150 mm
    compose("Figure8_organ_mr",
            rows=[["GSE282512_organ_mr_forest"]],
            target_w=int(round(150 / 25.4 * DPI)))

    print("done")
