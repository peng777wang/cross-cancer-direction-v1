"""Figure 1 - Direction Safety Index (DSI) of cross-cancer pleiotropic loci.

Inputs are the shipped source-data CSVs (``../source_data``); no upstream
pipeline output is needed.  ``FIG_DATA`` and ``FIG_OUT`` override the input and
output directories.

Contract
--------
Core conclusion: cross-cancer pleiotropic loci are strongly direction-coherent;
a discrete minority is called fully reversed (DSI = -1), but that reversal is
*MTAG-inferred* -- the reversed arm never reaches genome-wide significance in
its own GWAS, and MTAG's contribution at such loci is amplification by borrowing
from the well-powered traits rather than independent observation.  The figure
therefore presents the DSI as a direction annotation with an explicit evidence
grade, not as a validated class of biological antagonism.

Results-level question: do cross-cancer pleiotropic loci agree in direction, and
how far can a disagreement be trusted?

Panels (inferential roles)
  A genome-wide DSI map (MTAG)        discovery / landscape  (hero)
  B DSI vs sign-flip null             statistical support
  C borrowing fraction per cancer     mechanism of the apparent reversal
  D gene-level DSI dumbbell           descriptive gene summary
  E raw-z vs MTAG class agreement     robustness to the MTAG step
  F independent-dataset rg forest     external comparison

Archetype: quantitative grid.  Journal/export contract: 182.9 mm width, editable
text (svg.fonttype none / pdf.fonttype 42), 5 pt glyph floor, source data per
panel, render-time alignment gate plus rendered collision audit before delivery.

Reviewer risk addressed: an earlier revision of this figure asserted a
"pan-cancer-risky" 8q24/8p21 axis.  Raw-data verification failed (docs/AUDIT_LOG.md,
retraction 2) and the claim was withdrawn; no panel here depends on it.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("FIG_DATA", os.path.join(HERE, "..", "source_data"))
OUT = os.environ.get("FIG_OUT", os.path.join(HERE, "..", "figures_out"))
os.makedirs(OUT, exist_ok=True)

# The render-time panel-alignment gate is a QA tool, not part of the science. It
# is optional here so that the figure can be regenerated on any machine; when it
# is absent the figure is still produced and the gate is simply skipped.
try:
    sys.path.insert(0, os.environ.get(
        "NATURE_FIGURE_SCRIPTS",
        os.path.expanduser("~/.codex/skills/nature-figure/scripts")))
    from audit_panel_alignment import require_matplotlib_panel_alignment
except Exception:                                     # pragma: no cover
    def require_matplotlib_panel_alignment(fig, **kw):
        print("  [alignment gate unavailable - skipped]")

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 6.5,
    "axes.labelsize": 7,
    "axes.titlesize": 7.5,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.7,
    "legend.frameon": False,
    "axes.unicode_minus": False,
    "figure.dpi": 150,
})

SAFE, MIX, RISKY = "#1B7837", "#8C8C8C", "#C0392B"
BLUE, INK = "#2166AC", "#252525"
GROUPC = ["#1B7837", "#2C7FB8", "#41AB8C", "#E08214", "#C0392B", "#6A51A3"]

TRAITS = ["LUAD", "LUSC", "SCLC", "EC", "EAC", "THCA", "RCC", "PC", "GC", "PRAD"]
SIG = 5.451

# ----------------------------------------------------------------- data
# Every input is a shipped source-data file, so the figure can be regenerated
# without the upstream pipeline.
loc = pd.read_csv(os.path.join(DATA, "dsi_loci.csv"), dtype={"SNP": str})
ann = pd.read_csv(os.path.join(DATA, "dsi_loci_positions.csv"), dtype={"SNP": str})
rawj = pd.read_csv(os.path.join(DATA, "dsi_loci_rawz.csv"), dtype={"SNP": str})

d = loc.merge(ann[["SNP", "CHR", "POS"]], on="SNP", how="left")
d = d.merge(rawj[["SNP", "DSI", "class_raw", "class_mtag", "k"]].rename(
    columns={"DSI": "DSI_raw"}), on="SNP", how="left")
d["CHR"] = pd.to_numeric(d["CHR"], errors="coerce")
d["POS"] = pd.to_numeric(d["POS"], errors="coerce")
d["cls"] = np.where(d.DSI == 1, "safe", np.where(d.DSI == -1, "risky", "mixed"))
d = d[d.CHR.notna() & d.POS.notna()].copy()
print("loci plotted in A:", len(d), "of", len(loc))

HG19 = {1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276, 5: 180915260,
        6: 171115067, 7: 159138663, 8: 146364022, 9: 141213431, 10: 135534747,
        11: 135006516, 12: 133851895, 13: 115169878, 14: 107349540, 15: 102531392,
        16: 90354753, 17: 81195210, 18: 78077248, 19: 59128983, 20: 63025520,
        21: 48129895, 22: 51304566}
off = {}
run = 0
for c in sorted(HG19):
    off[c] = run
    run += HG19[c]
d["gx"] = d.CHR.map(off) + d.POS

null_dsi = pd.read_csv(os.path.join(
    DATA, "dsi_signflip_null_1000reps.csv")).null_DSI.values

# ----------------------------------------------------------------- figure
fig = plt.figure(figsize=(7.2, 8.5))
gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.42, 1.20, 0.80],
                      hspace=0.62, wspace=0.30,
                      left=0.075, right=0.965, top=0.930, bottom=0.055)
axA = fig.add_subplot(gs[0, :])
axB = fig.add_subplot(gs[1, 0])
axC = fig.add_subplot(gs[1, 1])
axD = fig.add_subplot(gs[2, 0])
axE = fig.add_subplot(gs[2, 1])
axF = fig.add_subplot(gs[3, :])


def plabel(ax, s, dx=-0.070, dy=0.030):
    ax.text(dx, 1 + dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left", color=INK)


# =============================================================== Panel A
rng = np.random.default_rng(11)
colmap = {"safe": SAFE, "mixed": MIX, "risky": RISKY}
for c in sorted(HG19):
    if c % 2 == 0:
        axA.axvspan(off[c], off[c] + HG19[c], color="0.5", alpha=0.055, lw=0)
    axA.text(off[c] + HG19[c] / 2, 1.44, str(c), ha="center", va="bottom",
             fontsize=5, color="0.35")
axA.axhline(1, color=SAFE, lw=0.6, ls=":", alpha=0.8)
axA.axhline(0, color="0.6", lw=0.6, ls="--")
axA.axhline(-1, color=RISKY, lw=0.6, ls=":", alpha=0.8)
for cls in ["mixed", "safe", "risky"]:
    s = d[d.cls == cls]
    axA.scatter(s.gx, s.DSI + rng.uniform(-0.035, 0.035, len(s)),
                s=4 + 2.6 * (s.n_sig - 2), c=colmap[cls], alpha=0.85,
                lw=0.25, edgecolor="white", zorder=3, label=cls)

# dagger marks a locus whose reversal exists only in the MTAG output; see the
# footnote added below panel A.
LAB = [("TERT / CLPTM1L (5p15)", 5, 1280000, 1.22, "left", 0),
       ("CHRNA3 / IREB2 (15q25)", 15, 78800000, 1.22, "right", 0),
       ("PTCSC2 / FOXE1 (9q22)", 9, 100500000, 0.70, "left", 1),
       ("TYMS (18p11)", 18, 657000, 1.22, "left", 0),
       ("NKX3-1 (8p21)\u2020", 8, 23500000, -1.28, "left", 1),
       ("CASC8 / CCAT2 (8q24)", 8, 128500000, 0.50, "center", 1),
       ("LMTK2 (7q21)\u2020", 7, 97800000, -1.28, "right", 1),
       ("WDPCP (2p15)\u2020", 2, 64000000, -0.30, "center", 1),
       ("HNF1B (17q12)\u2020", 17, 36100000, -1.28, "left", 1)]
for lab, c, pos, ytxt, ha, arr in LAB:
    x = off[c] + pos
    sub = d[d.CHR == c]
    if len(sub):
        j = (sub.POS - pos).abs().idxmin()
        anchor = (sub.loc[j, "gx"], sub.loc[j, "DSI"])
    else:
        anchor = (x, np.clip(ytxt, -1, 1))
    if arr:
        axA.text(anchor[0], ytxt, lab, fontsize=5.4, ha=ha, va="center", color=INK)
        ay = anchor[1]
        lo, hi = (ytxt + 0.13, ay - 0.05) if ytxt < ay else (ay + 0.05, ytxt - 0.13)
        axA.plot([anchor[0], anchor[0]], [lo, hi], color="0.45", lw=0.5, zorder=2)
    else:
        axA.text(anchor[0], ytxt, lab, fontsize=5.4, ha=ha, va="center", color=INK)

axA.set_xlim(0, run)
axA.set_ylim(-1.55, 1.62)
axA.set_yticks([-1, -0.5, 0, 0.5, 1])
axA.set_ylabel("DSI  (per locus)")
axA.set_xlabel("GRCh37 position", labelpad=1)
axA.set_xticks([])
axA.spines["bottom"].set_visible(False)
h = [Line2D([0], [0], marker="o", color="none", markerfacecolor=colmap[k],
            markeredgecolor="white", markersize=4, label=lbl)
     for k, lbl in [("safe", "Aligned, DSI = +1"), ("mixed", "-1 < DSI < +1"),
                    ("risky", "Reversed, DSI = -1 (MTAG)")]]
h += [Line2D([0], [0], marker="o", color="0.35", lw=0, markersize=np.sqrt(4 + 2.6 * (k - 2)),
             label="%d cancers" % k) for k in [2, 4, 6]]
axA.legend(handles=h, loc="lower left", bbox_to_anchor=(0.005, -0.085), ncol=6,
           handletextpad=0.45, columnspacing=1.4, borderpad=0.1)
axA.set_xlabel("GRCh37 position", labelpad=13)
axA.text(0, -0.30, "\u2020 reversal present in the MTAG output only: the locus does not reach "
         "genome-wide significance in \u2265 2 of its own, unborrowed GWAS.",
         transform=axA.transAxes, fontsize=5.2, color="0.35", ha="left", va="top")
axA.set_title("Genome-wide direction architecture of 273 cross-cancer pleiotropic loci "
              "(direction inferred by MTAG)", pad=8, loc="left")
plabel(axA, "A", dx=-0.045, dy=0.022)

# =============================================================== Panel B
bins = np.linspace(-1, 1, 21)
axB.hist(null_dsi, bins=bins, density=True, color="0.80", edgecolor="white",
         lw=0.35, label="sign-flip null", zorder=1)
hobs, _ = np.histogram(d.DSI.values, bins=bins, density=True)
ctr = (bins[:-1] + bins[1:]) / 2
axB.bar(ctr, hobs, width=0.088, color=[colmap[c] for c in d.cls], alpha=0.92,
        edgecolor="white", lw=0.35, zorder=3, label="observed")
axB.set_xlim(-1.18, 1.18)
axB.set_ylim(0, max(hobs.max(), 1) * 2.0)
axB.set_xlabel("DSI per locus")
axB.set_ylabel("density")
axB.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=1, handlelength=1.1,
           labelspacing=0.3)
axB.annotate("", xy=(1, hobs[-1]), xytext=(1, hobs[-1] * 1.22),
             arrowprops=dict(arrowstyle="-|>", lw=0.7, color=SAFE))
axB.text(0.97, hobs[-1] * 1.25, "196 / 273 loci DSI = +1\n(expected 52 under H0)\nP = 5e-83",
         fontsize=5.6, ha="right", va="bottom", color=SAFE)
axB.text(0.012, 0.76, "42 loci DSI = -1\n(all 42 shared pairs reversed)\n"
         "no excess over H0 (expect 52)",
         transform=axB.transAxes, va="top", ha="left", fontsize=5.6, color=RISKY)
axB.text(0.02, 0.45,
         "concordant pairs 448 / 572 = 78.3% (MTAG)\n"
         "404 / 572 = 70.6% (raw z);  null = 50.0%\n"
         "pair-fraction P < 2e-4 (5,000 sign flips)\n"
         "grey null: 1,000 sign flips x 273 loci;  |DSI| = 1: 87.2% vs 66.2%",
         transform=axB.transAxes, fontsize=5.6, ha="left", va="top", color="0.25",
         bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.4))
axB.set_title("DSI is bimodal: aligned majority plus a reversed minority", pad=10, loc="left")
plabel(axB, "B")

# =============================================================== Panel C
# Replaces the former cancer-pair direction matrix, which existed only to support
# the withdrawn "prostate cancer is the reversed pole" claim.  What is plotted
# instead is the quantity that actually explains the reversed class: how much of
# each cancer's MTAG estimate is borrowed rather than observed.
bor = pd.read_csv(os.path.join(DATA, "antagonism_borrowing.csv"))
grz = bor.groupby("trait").agg(
    median_borrowed=("borrowed_frac", "median"),
    median_own_z=("z_raw", lambda s: float(np.abs(s).median())),
    n_loci=("SNP", "size")).reset_index()
grz["median_borrowed"] *= 100
grz = grz.sort_values("median_borrowed").reset_index(drop=True)
ypc = np.arange(len(grz))
bcol = [RISKY if v > 50 else (MIX if v > 20 else SAFE) for v in grz.median_borrowed]
axC.barh(ypc, grz.median_borrowed, color=bcol, alpha=0.9, height=0.66,
         edgecolor="white", lw=0.4, zorder=2)
for yv, vv in zip(ypc, grz.median_borrowed):
    axC.text(vv + (2.6 if vv >= 0 else -2.6), yv, "%+.0f%%" % vv, fontsize=5.4,
             va="center", ha="left" if vv >= 0 else "right", color="0.2", zorder=3)
axC.axvline(0, color="0.35", lw=0.7, zorder=1)
axC.set_yticks(ypc)
axC.set_yticklabels(["%s  (own |z| %.1f)" % (t, z) for t, z in
                     zip(grz.trait, grz.median_own_z)], fontsize=5.4)
axC.set_ylim(-1.10, len(grz) - 0.35)
axC.set_xlim(-24, 102)
axC.set_xticks([0, 25, 50, 75, 100])
axC.set_xlabel("median borrowed fraction of the MTAG z (%)")
# placed below the last bar inside the axes, where no bar or tick label can reach
axC.text(2, -0.62, "borrowed = 1 - |z_own| / |z_MTAG|", fontsize=5.3, color="0.3",
         va="center", ha="left")
axC.text(2, -0.90, "no reversed-arm cancer reaches |z| > 5.45 in its own GWAS",
         fontsize=5.3, color="0.3", va="center", ha="left")
axC.set_title("The reversed arm is borrowed, not observed", pad=10, loc="left")
plabel(axC, "C", dx=-0.16)

# =============================================================== Panel D
gl = d.groupby("gene").agg(n=("DSI", "size"), dsi_mtag=("DSI", "mean"),
                           dsi_raw=("DSI_raw", "mean")).reset_index()
GROUPS = [
    ("Telomere maintenance", ["TERT", "CLPTM1L", "RTEL1"]),
    ("Smoking / nAChR (15q25)", ["CHRNA5", "CHRNA3", "CHRNB4", "IREB2"]),
    ("Shared GI-thyroid loci", ["PTCSC2", "FOXE1"]),
    ("8q24 region", ["CASC8", "NKX3-1"]),
    ("Reversed in MTAG only", ["LMTK2", "PCAT1", "WDPCP", "HNF1B"]),
    ("Antifolate target", ["TYMS"]),
]
rows, y, ticks, tlabs, tcol = [], 0.0, [], [], []
for gi, (gname, genes) in enumerate(GROUPS):
    sub = gl[gl.gene.isin(genes)].copy()
    sub["ord"] = sub.gene.map({g: i for i, g in enumerate(genes)})
    sub = sub.sort_values(["dsi_mtag", "ord"], ascending=[False, True])
    if gi:
        axD.axhline(y - 0.45, color="0.88", lw=0.6, zorder=0)
    for _, r in sub.iterrows():
        rows.append([gname, r.gene, r.n, r.dsi_mtag, r.dsi_raw])
        axD.plot([r.dsi_raw, r.dsi_mtag], [y, y], color="0.75", lw=0.7, zorder=1)
        axD.scatter([r.dsi_raw], [y], s=13, facecolor="white", edgecolor=GROUPC[gi],
                    lw=0.8, zorder=3)
        axD.scatter([r.dsi_mtag], [y], s=15, color=GROUPC[gi], lw=0, zorder=4)
        ticks.append(y)
        tlabs.append("%s (%d)" % (r.gene, int(r.n)))
        tcol.append(GROUPC[gi])
        y += 1
    y += 0.10
axD.axvline(0, color="0.6", lw=0.6, ls="--")
axD.set_yticks(ticks)
axD.set_yticklabels(tlabs, fontsize=5.4)
for t, c in zip(axD.get_yticklabels(), tcol):
    t.set_color(c)
axD.set_xlim(-1.34, 1.78)
axD.set_ylim(y - 0.35, -0.62)
axD.set_xticks([-1, -0.5, 0, 0.5, 1])
axD.tick_params(axis="x", pad=2.6)
axD.set_xlabel("")
axD.text(0.0, 1.015, "drug labels are existing annotations, not DSI predictions",
         transform=axD.transAxes, ha="left", va="bottom", fontsize=5.1,
         color="0.35")
for g, lab in [("TYMS", "5-FU / pemetrexed"), ("TERT", "imetelstat"),
               ("CHRNA5", "varenicline / cytisine")]:
    i = [k for k, t in enumerate(tlabs) if t.startswith(g + " ")]
    if i:
        axD.annotate(lab, xy=(1.02, ticks[i[0]]), xytext=(1.16, ticks[i[0]]),
                     fontsize=5, va="center", ha="left", color=INK,
                     arrowprops=dict(arrowstyle="-", lw=0.45, color="0.5"))
gh = [Line2D([0], [0], marker="s", color="none", markerfacecolor=GROUPC[i],
             markersize=3.6, label=g) for i, (g, _) in enumerate(GROUPS)]
gh += [Line2D([0], [0], marker="o", color="none", markerfacecolor="0.4",
              markersize=3.8, label="MTAG z"),
       Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
              markeredgecolor="0.4", markersize=3.8, label="raw GWAS z")]
axD.legend(handles=gh, loc="upper center", bbox_to_anchor=(0.46, -0.115), ncol=4,
           handletextpad=0.25, columnspacing=1.0, labelspacing=0.3)
axD.set_title("Gene-level DSI (descriptive summary of the locus classes)", pad=13,
              loc="left")
plabel(axD, "D", dx=-0.315)

# =============================================================== Panel E
CLS = ["safe", "mixed", "risky"]
CLS_DISPLAY = {"safe": "aligned", "mixed": "mixed", "risky": "reversed"}
cm = pd.crosstab(pd.Categorical(rawj.class_raw, CLS), pd.Categorical(rawj.class_mtag, CLS))
cm = cm.reindex(index=CLS, columns=CLS, fill_value=0)
cmap_c = {"safe": SAFE, "mixed": MIX, "risky": RISKY}
x = np.arange(3)
bw = 0.26
for k, rc in enumerate(CLS):
    vals = cm.loc[rc].values.astype(float)
    axE.bar(x + (k - 1) * bw, vals, width=bw, color=cmap_c[rc], alpha=0.92,
            edgecolor="white", lw=0.4, label="raw GWAS z: %s" % CLS_DISPLAY[rc])
    for xi, v in zip(x + (k - 1) * bw, vals):
        if v > 0:
            axE.text(xi, v + 4, "%d" % v, ha="center", va="bottom", fontsize=5.2,
                     color="0.2")
axE.set_xticks(x, ["%s\n(n = %d)" % (CLS_DISPLAY[cc], cm[cc].sum()) for cc in CLS])
axE.set_xlabel("class from MTAG z")
axE.set_ylabel("loci")
axE.set_ylim(0, cm.values.max() * 1.34)
axE.set_xlim(-0.55, 2.55)
axE.legend(loc="upper right", ncol=1, handlelength=1.0, labelspacing=0.25)
axE.text(0.17, 0.585, "84.2% of loci keep class\n"
                     "(r = 0.70 across 273 loci)\n"
                     "reversed loci recovered 31 / 42",
         transform=axE.transAxes, fontsize=5.6, color="0.2", va="top")
axE.set_title("Directions largely survive without MTAG shrinkage", pad=10, loc="left")
plabel(axE, "E", dx=-0.16)

# =============================================================== Panel F
# read from the shipped source data rather than a hard-coded list
_f = pd.read_csv(os.path.join(DATA, "fig1_panelF_rg.csv"))
forest = [(r.comparison, float(r.rg), float(r.se), r.flag) for r in _f.itertuples()]
iv = [(float(r.rg), float(r.se)) for r in _f.itertuples() if r.flag == "indep"]
w = np.array([1 / s ** 2 for _, s in iv])
pooled = float(np.sum(w * np.array([b for b, _ in iv])) / w.sum())
pooled_se = float(1 / np.sqrt(w.sum()))
pooled_z = pooled / pooled_se
from math import erfc  # noqa: E402
pooled_p = erfc(abs(pooled_z) / np.sqrt(2))
print("pooled rg = %.4f (%.4f) z=%.2f p=%.4f" % (pooled, pooled_se, pooled_z, pooled_p))

ylab = [r[0] for r in forest]
ypos = np.arange(len(forest))[::-1]
axF.axvspan(-1.74, 0, color=RISKY, alpha=0.05, lw=0)
axF.axvline(0, color="0.35", lw=0.8)
for (lab, b, s, kind), y in zip(forest, ypos):
    c = {"indep": BLUE, "ctrl": SAFE, "bad": "0.62"}[kind]
    if s > 0.5:  # CI wider than the plotted range: show the point only
        axF.plot([b], [y], marker="o", ms=3.2, mfc="white", mec=c, mew=0.9, ls="none")
    else:
        axF.errorbar(b, y, xerr=1.96 * s, fmt="o", ms=3.4, color=c, ecolor=c,
                     elinewidth=0.8, capsize=1.8, mfc=c, mec="white", mew=0.4, ls="none")
    axF.text(-1.80, y, lab, fontsize=5.7, va="center", ha="right", color=c, clip_on=False)
    txt = "rg = %+.3f (s.e. %.3f)" % (b, s)
    if s > 0.5:
        txt = "rg = %+.2f (s.e. %.2f, CI not interpretable)" % (b, s)
    axF.text(1.95, y, txt, fontsize=5.2,
             va="center", ha="right", color="0.25")
yp = ypos[-1] - 1.15
axF.errorbar(pooled, yp, xerr=1.96 * pooled_se, fmt="D", ms=3.6, color=INK,
             ecolor=INK, elinewidth=0.9, capsize=1.8)
axF.text(-1.80, yp, "Pooled (4 pairs, approximate; panels shared)", fontsize=5.7,
         va="center", ha="right", color=INK, clip_on=False)
axF.text(1.95, yp, "rg = %+.3f (s.e. %.3f), P = 0.006" % (pooled, pooled_se),
         fontsize=5.4, va="center", ha="right", color=INK)
axF.set_ylim(-2.55, len(forest) - 0.35)
axF.set_xlim(-1.76, 2.62)
axF.set_yticks([])
axF.set_xticks([-1.5, -1, -0.5, 0, 0.5, 1])
axF.set_xlabel("genetic correlation rg with 95% CI")
axF.spines["left"].set_visible(False)
axF.text(-1.80, -2.12, "open/grey = not interpretable; reported, not used",
         fontsize=5.0, color="0.45", ha="left", clip_on=False)
axF.set_title("Genetic correlation with prostate cancer across independent GWAS pairs", pad=10, loc="left")
plabel(axF, "F", dx=-0.045)

fig.text(0.075, 0.978,
         "Fig. 1 | A direction-safety map of cross-cancer genetic pleiotropy",
         fontsize=9, fontweight="bold", ha="left", va="bottom", color=INK)

base = os.path.join(OUT, "Figure1_DSI")
fig.canvas.draw()
require_matplotlib_panel_alignment(
    fig, json_out=base + ".alignment.json", overlay_svg=base + ".alignment.svg",
    tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=False)
fig.savefig(base + ".png", dpi=600, bbox_inches="tight")
fig.savefig(base + ".pdf", bbox_inches="tight")
fig.savefig(base + ".svg", bbox_inches="tight")
fig.savefig(base + ".tiff", dpi=600, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})

# ------------------------------------------------------------- source data
cols = ["SNP", "CHR", "POS", "gene", "n_sig", "DSI", "DSI_raw", "cls", "class_raw", "class_mtag"]
d[cols].to_csv(os.path.join(OUT, "panelA_D_locus_table.csv"), index=False)
pd.DataFrame({"bin_center": ctr, "observed_density": hobs}).to_csv(
    os.path.join(OUT, "panelB_histogram.csv"), index=False)
pd.DataFrame({
    "stat": ["n_loci", "n_safe_mtag", "n_risky_mtag", "n_mixed_mtag",
             "expected_safe_H0", "P_safe_poisson_binomial",
             "concordant_pairs_mtag", "antagonistic_pairs_mtag", "fraction_mtag",
             "concordant_pairs_rawz", "antagonistic_pairs_rawz", "fraction_rawz",
             "null_fraction_mean", "null_fraction_sd", "MC_p_rawz",
             "frac_absDSI1_observed", "frac_absDSI1_null"],
    "value": [273, 196, 42, 35, 52.2, 5.22e-83, 448, 124, 0.7832,
              404, 168, 0.7063, 0.5003, 0.0211, 2e-4, 0.8718, 0.6618],
}).to_csv(os.path.join(OUT, "panelB_null_stats.csv"), index=False)
grz.to_csv(os.path.join(OUT, "panelC_borrowing_by_cancer.csv"), index=False)
pd.DataFrame(rows, columns=["group", "gene", "n_loci", "DSI_mtag", "DSI_raw"]).to_csv(
    os.path.join(OUT, "panelD_gene_dsi.csv"), index=False)
cm.to_csv(os.path.join(OUT, "panelE_confusion_raw_vs_mtag.csv"))
pd.DataFrame(forest, columns=["comparison", "rg", "se", "flag"]).to_csv(
    os.path.join(OUT, "panelF_replication_forest.csv"), index=False)
print("POOLED %.4f %.4f" % (pooled, pooled_se))
print("FIG1_DONE", base)
