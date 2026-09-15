"""Rebuild the resource JSON/CSV from the shipped source data.

Run from anywhere: paths are resolved relative to this file.

    python release/resource/build_resource.py
    python release/resource/build_explorer.py     # then re-inline into the HTML

Inputs : ../source_data/dsi_loci.csv, dsi_loci_rawz.csv, dsi_loci_positions.csv,
         mtag_z/raw_z per locus per cancer (shipped in source_data), and the
         bundled UCSC hg19 cytoband table
Outputs: pleiotropy_resource.json/.csv, locus_summary.csv in this directory
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "source_data")
TRAITS = ["LUAD", "LUSC", "SCLC", "EC", "EAC", "THCA", "RCC", "PC", "GC", "PRAD"]
GWS_Z = 5.451


def main():
    loc = pd.read_csv(os.path.join(DATA, "dsi_loci.csv"), dtype={"SNP": str})
    pos = pd.read_csv(os.path.join(DATA, "dsi_loci_positions.csv"), dtype={"SNP": str})
    raw = pd.read_csv(os.path.join(DATA, "dsi_loci_rawz.csv"), dtype={"SNP": str})
    zt = pd.read_csv(os.path.join(DATA, "mtag_z_pleio.tsv"), sep="\t",
                     dtype={"SNP": str}).set_index("SNP")
    zr = pd.read_csv(os.path.join(DATA, "raw_z_pleio.tsv"), sep="\t",
                     dtype={"SNP": str}).set_index("SNP")
    cyto = json.load(open(os.path.join(HERE, "cytoband_hg19.json"), encoding="utf-8"))

    loc = loc.merge(pos[["SNP", "CHR", "POS"]], on="SNP", how="left")
    loc = loc.merge(raw[["SNP", "class_raw", "class_mtag", "k"]], on="SNP", how="left")
    loc["cls"] = np.where(loc.DSI == 1, "aligned",
                          np.where(loc.DSI == -1, "reversed (MTAG-inferred)", "mixed"))

    rows = []
    for r in loc.itertuples():
        snp = r.SNP
        if snp not in zt.index or snp not in zr.index:
            continue
        rec = {"snp": snp, "chr": None if pd.isna(r.CHR) else int(r.CHR),
               "pos": None if pd.isna(r.POS) else int(r.POS), "gene": r.gene,
               "k": int(r.k), "n_conc": int(r.n_conc), "n_ant": int(r.n_ant),
               "dsi_mtag": round(float(r.DSI), 4),
               "dsi_raw": round(float(raw.loc[raw.SNP == snp, "DSI"].iloc[0]), 4),
               "cls": r.cls, "cls_raw": r.class_raw, "trait": {}}
        for t in TRAITS:
            zm = float(zt.loc[snp, t]) if t in zt.columns else np.nan
            zq = float(zr.loc[snp, t]) if t in zr.columns else np.nan
            if not np.isfinite(zm):
                continue
            borrow = None
            if np.isfinite(zq) and abs(zm) > 1e-9:
                borrow = round(float(1 - abs(zq) / abs(zm)), 4)
            rec["trait"][t] = {"z_mtag": round(zm, 3),
                               "z_raw": None if not np.isfinite(zq) else round(zq, 3),
                               "borrowed": borrow, "mtag_sig": bool(abs(zm) > GWS_Z),
                               "own_sig": bool(np.isfinite(zq) and abs(zq) > GWS_Z)}
        med = [v["borrowed"] for v in rec["trait"].values()
               if v["borrowed"] is not None and v["mtag_sig"]
               and v["z_mtag"] * v["z_raw"] < 0]
        rec["median_borrowed_if_antagonistic"] = round(float(np.median(med)), 3) if med else None
        rec["own_gwas_supports"] = bool(sum(v["own_sig"] for v in rec["trait"].values()) >= 2)
        rows.append(rec)

    loci = sorted(rows, key=lambda x: (x["chr"] or 99, x["pos"] or 0))
    meta = json.load(open(os.path.join(HERE, "resource_meta.json"), encoding="utf-8"))
    meta.update({"n_loci": len(loci), "traits": TRAITS, "gws_z": GWS_Z,
                 "cytobands": cyto})
    payload = {"meta": meta, "loci": loci}

    for path in [os.path.join(HERE, "pleiotropy_resource.json"),
                 os.path.join(DATA, "pleiotropy_resource.json")]:
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, separators=(",", ":"))

    flat = [{"snp": r["snp"], "gene": r["gene"], "chr": r["chr"], "pos": r["pos"],
             "class_mtag": r["cls"], "DSI_mtag": r["dsi_mtag"], "DSI_raw": r["dsi_raw"],
             "locus_k": r["k"], "trait": t, "z_mtag": v["z_mtag"], "z_raw": v["z_raw"],
             "borrowed_fraction": v["borrowed"], "mtag_significant": v["mtag_sig"],
             "own_gwas_significant": v["own_sig"]}
            for r in loci for t, v in r["trait"].items()]
    pd.DataFrame(flat).to_csv(os.path.join(HERE, "pleiotropy_resource.csv"), index=False)
    pd.DataFrame(flat).to_csv(os.path.join(DATA, "pleiotropy_resource.csv"), index=False)
    summ = pd.DataFrame([{k: v for k, v in r.items() if k != "trait"} for r in loci])
    summ.to_csv(os.path.join(HERE, "locus_summary.csv"), index=False)
    summ.to_csv(os.path.join(DATA, "locus_summary.csv"), index=False)

    print("loci %d | locus x trait %d | aligned %d | reversed %d | mixed %d"
          % (len(loci), len(flat),
             sum(1 for r in loci if r["cls"] == "aligned"),
             sum(1 for r in loci if r["cls"].startswith("reversed")),
             sum(1 for r in loci if r["cls"] == "mixed")))


if __name__ == "__main__":
    main()
