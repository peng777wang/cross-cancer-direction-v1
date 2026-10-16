"""Assemble the cross-cancer direction resource (v1.1, 2026-10-13).

Change from v1.0: the z-scores are stored to six decimals and the borrowing
fraction is computed from the stored (rounded) values, so that
1 - |z_raw| / |z_mtag| reproduces the stored fraction for every locus-cancer
pair.  In v1.0 the z-scores were stored to three decimals while the fraction
used full precision, so 156 of 2,730 records did not satisfy the stated
formula when it was recomputed from the released file.

Outputs (into release/resource/ and release/source_data/):
  pleiotropy_resource.json   every locus x cancer: MTAG z, raw z, borrowing
  pleiotropy_resource.csv    the same data flat, one row per locus x cancer
  locus_summary.csv          one row per locus with class and summary numbers
  pleiotropy_explorer.html   single-file offline explorer (built separately)

Definitions are those of the manuscript:
  DSI        sign agreement across the cancers in which the locus is significant
  borrowing  1 - |z_raw| / |z_MTAG|, the share of that cancer's MTAG estimate
             that is not present in that cancer's own GWAS
Nothing is recomputed differently here; this script only reshapes the audited
outputs so they can be queried.
"""
import json
import os

import numpy as np
import pandas as pd

OUT = "outputs"
RES = "release/resource"
SRC = "release/source_data"
TRAITS = ["LUAD", "LUSC", "SCLC", "EC", "EAC", "THCA", "RCC", "PC", "GC", "PRAD"]
GWS_Z = 5.451
Z_DP = 6          # decimals kept for z_mtag / z_raw
BORROW_DP = 6     # decimals kept for the borrowing fraction


def main():
    os.makedirs(RES, exist_ok=True)
    os.makedirs(SRC, exist_ok=True)

    loc = pd.read_csv(os.path.join(OUT, "dsi_loci.csv"), dtype={"SNP": str})
    pos = pd.read_csv(os.path.join(SRC, "dsi_loci_positions.csv"),
                      dtype={"SNP": str})
    raw = pd.read_csv(os.path.join(OUT, "dsi_loci_rawz.csv"), dtype={"SNP": str})
    zt = pd.read_csv(os.path.join(OUT, "mtag_z_pleio.tsv"), sep="\t",
                     dtype={"SNP": str})
    zr = pd.read_csv(os.path.join(OUT, "raw_z_pleio.tsv"), sep="\t",
                     dtype={"SNP": str})
    for d in (zt, zr):
        d.set_index("SNP", inplace=True)

    cyto = json.load(open("work/cytoBand_hg19.json", encoding="utf-8"))

    loc = loc.merge(pos[["SNP", "CHR", "POS"]], on="SNP", how="left")
    loc = loc.merge(raw[["SNP", "class_raw", "class_mtag", "k"]], on="SNP",
                    how="left")
    loc["cls"] = np.where(loc.DSI == 1, "aligned",
                          np.where(loc.DSI == -1, "reversed (MTAG-inferred)",
                                   "mixed"))

    rows = []
    worst = 0.0
    n_bad = 0
    for r in loc.itertuples():
        snp = r.SNP
        if snp not in zt.index or snp not in zr.index:
            continue
        rec = {"snp": snp, "chr": None if pd.isna(r.CHR) else int(r.CHR),
               "pos": None if pd.isna(r.POS) else int(r.POS),
               "gene": r.gene, "k": int(r.k), "n_conc": int(r.n_conc),
               "n_ant": int(r.n_ant),
               "dsi_mtag": round(float(r.DSI), 4),
               "dsi_raw": round(float(raw.loc[raw.SNP == snp, "DSI"].iloc[0]), 4),
               "cls": r.cls, "cls_raw": r.class_raw,
               "trait": {}}
        for t in TRAITS:
            zm_full = float(zt.loc[snp, t]) if t in zt.columns else np.nan
            zr_full = float(zr.loc[snp, t]) if t in zr.columns else np.nan
            if not np.isfinite(zm_full):
                continue
            zm = round(zm_full, Z_DP)
            zr_ = None if not np.isfinite(zr_full) else round(zr_full, Z_DP)
            borrow = None
            if zr_ is not None and abs(zm) > 1e-9:
                borrow = round(1 - abs(zr_) / abs(zm), BORROW_DP)
                back = 1 - abs(zr_) / abs(zm)
                worst = max(worst, abs(back - borrow))
                if abs(back - borrow) > 5e-7:
                    n_bad += 1
            rec["trait"][t] = {
                "z_mtag": zm, "z_raw": zr_, "borrowed": borrow,
                "mtag_sig": bool(abs(zm) > GWS_Z),
                "own_sig": bool(zr_ is not None and abs(zr_) > GWS_Z)}
        med = [v["borrowed"] for v in rec["trait"].values()
               if v["borrowed"] is not None and v["mtag_sig"]
               and v["z_mtag"] * v["z_raw"] < 0]
        rec["median_borrowed_if_antagonistic"] = (
            round(float(np.median(med)), 3) if med else None)
        rec["own_gwas_supports"] = bool(
            sum(v["own_sig"] for v in rec["trait"].values()) >= 2)
        rows.append(rec)

    loci = sorted(rows, key=lambda x: (x["chr"] or 99, x["pos"] or 0))
    payload = {
        "meta": {
            "title": "Cross-cancer direction and MTAG-borrowing explorer",
            "version": "1.2",
            "date": "2026-10-13",
            "z_decimals": Z_DP,
            "n_loci": len(loci),
            "traits": TRAITS,
            "gws_z": GWS_Z,
            "cytobands": cyto,
            "cytoband_source": "UCSC hg19 cytoBand (GRCh37), matching the "
                               "coordinates used for the loci",
            "definitions": {
                "DSI": "sign agreement of a locus across the cancers in which "
                       "it is genome-wide significant; +1 = all pairs "
                       "concordant, -1 = all pairs antagonistic",
                "borrowed": "1 - |z_own| / |z_MTAG|; the share of a cancer's "
                            "MTAG estimate that is not present in that cancer's "
                            "own GWAS. z-scores are stored to six decimals and "
                            "the fraction is computed from the stored values, "
                            "so it reproduces the formula exactly.",
                "own_gwas_supports": "the locus reaches genome-wide "
                                     "significance in at least two of its own, "
                                     "unborrowed GWAS",
            },
            "caveats": [
                "Direction classes are inferred from MTAG output. In the "
                "discovery-panel reversed class every locus is a two-cancer "
                "locus, and the reversed-arm cancer never reaches genome-wide "
                "significance in its own GWAS (median borrowing 73.9-79.3%) "
                "while the concordant-arm cancer always does, so the "
                "reversed-arm direction is a model prediction rather than an "
                "observation. The replication panel contains 18 reversed loci "
                "of a different kind, in which the opposing arms are each "
                "supported by their own data.",
                "Borrowing is not a defect of MTAG; it is what MTAG is for. It "
                "is reported so that a reader can see which estimates are "
                "observations.",
                "Loci were selected for cross-cancer pleiotropy, so the excess "
                "of aligned loci is an upper bound on concordance, not an "
                "unbiased estimate; the reversed class shows no excess over "
                "the sign-flip null.",
                "Effect sizes are on the standardised (z) scale and are not "
                "comparable across cancers with different sample sizes.",
                "Significance flags are computed from the stored six-decimal "
                "z-scores. A value within 5e-7 of the genome-wide threshold "
                "can therefore differ from the flag in the manuscript tables, "
                "which use full precision.",
            ],
            "source": "public GWAS summary statistics; see the accompanying "
                      "manuscript and code release",
        },
        "loci": loci,
    }

    for d in (RES, SRC):
        with open(os.path.join(d, "pleiotropy_resource.json"), "w",
                  encoding="utf-8") as fh:
            json.dump(payload, fh, separators=(",", ":"))

    flat = []
    for r in loci:
        for t, v in r["trait"].items():
            flat.append({"snp": r["snp"], "gene": r["gene"], "chr": r["chr"],
                         "pos": r["pos"], "class_mtag": r["cls"],
                         "DSI_mtag": r["dsi_mtag"], "DSI_raw": r["dsi_raw"],
                         "locus_k": r["k"], "trait": t,
                         "z_mtag": v["z_mtag"], "z_raw": v["z_raw"],
                         "borrowed_fraction": v["borrowed"],
                         "mtag_significant": v["mtag_sig"],
                         "own_gwas_significant": v["own_sig"]})
    flat = pd.DataFrame(flat)
    flat.to_csv(os.path.join(SRC, "pleiotropy_resource.csv"), index=False)
    flat.to_csv(os.path.join(RES, "pleiotropy_resource.csv"), index=False)

    summ = pd.DataFrame([{k: v for k, v in r.items() if k != "trait"}
                         for r in loci])
    summ.to_csv(os.path.join(SRC, "locus_summary.csv"), index=False)
    summ.to_csv(os.path.join(RES, "locus_summary.csv"), index=False)

    print("loci:", len(loci), "| locus x trait rows:", len(flat))
    print("aligned %d | reversed %d | mixed %d"
          % (sum(1 for r in loci if r["cls"].startswith("aligned")),
             sum(1 for r in loci if r["cls"].startswith("reversed")),
             sum(1 for r in loci if r["cls"] == "mixed")))
    print("loci whose own GWAS supports them (>=2):",
          sum(1 for r in loci if r["own_gwas_supports"]))
    print("borrowing reproducibility: max deviation %.3e, rows >5e-7: %d"
          % (worst, n_bad))
    print("JSON size: %.0f KB" % (os.path.getsize(
        os.path.join(RES, "pleiotropy_resource.json")) / 1024))


if __name__ == "__main__":
    main()
