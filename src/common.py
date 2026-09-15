"""Shared, individually-verified primitives used by every analysis script.

Keeping these in one place matters: allele harmonisation and PLINK parsing were
each re-implemented several times during the project, and each re-implementation
was a chance for the sign convention to drift. Every function here has a
corresponding verification in docs/AUDIT_LOG.md.
"""
import gzip
import os

import numpy as np
import pandas as pd

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


# --------------------------------------------------------------------------- alleles
def same_allele(a, b):
    """True if two allele strings denote the same allele (strand-agnostic)."""
    a, b = str(a).upper(), str(b).upper()
    return a == b or a == COMPLEMENT.get(b, "?")


def harmonise(q1, q2, r1, r2):
    """Map a query effect allele pair onto a reference allele pair.

    Returns (flip, ambiguous, ok):
      flip      -- True if the query effect allele equals the reference *second* allele
      ambiguous -- True for A/T and C/G pairs, where strand cannot be resolved
                   by allele matching alone (resolve with allele frequency, or drop)
      ok        -- False if the alleles do not match the reference at all
    """
    q1, q2, r1, r2 = (str(x).upper() for x in (q1, q2, r1, r2))
    noflip = same_allele(q1, r1) and same_allele(q2, r2)
    flip = same_allele(q1, r2) and same_allele(q2, r1)
    if noflip and not flip:
        return False, False, True
    if flip and not noflip:
        return True, False, True
    if noflip and flip:
        return None, True, True
    return None, False, False


def harmonise_z(z, effect_allele, other_allele, ref_a1, ref_a2, freq=None,
                ref_freq=None, freq_tol=0.15):
    """Harmonise a signed z-score onto the reference A1 allele.

    Ambiguous variants are resolved with allele frequency when both the query
    frequency and the reference frequency are supplied; otherwise they are
    dropped (returns None).
    """
    flip, ambiguous, ok = harmonise(effect_allele, other_allele, ref_a1, ref_a2)
    if not ok:
        return None
    if ambiguous:
        if freq is None or ref_freq is None:
            return None
        d_no = abs(ref_freq - freq)
        d_fl = abs(ref_freq - (1 - freq))
        flip = d_fl < d_no
    if np.isfinite(ref_freq) and np.isfinite(freq) and not ambiguous:
        # frequency must agree in the chosen orientation
        expect = (1 - freq) if flip else freq
        if abs(ref_freq - expect) > freq_tol:
            return None
    return float(z) * (-1.0 if flip else 1.0)


# --------------------------------------------------------------------------- PLINK
_LUT = np.zeros((256, 4), dtype=np.uint8)
for _v in range(256):
    _LUT[_v] = [(_v >> 6) & 3, (_v >> 4) & 3, (_v >> 2) & 3, _v & 3]


def read_plink_bed_snps(panel, snp_indices, n_samples=None):
    """Return dosages (n_snps x n_samples) counting the .bim A1 allele.

    Verified convention (docs/AUDIT_LOG.md, item 2): byte code 00 is homozygous
    for the *first* allele in the .bim file (A1), 10 is heterozygous, 11 is
    homozygous for A2, 01 is missing. Dosages are returned as the A1 count.
    """
    fam = panel + ".fam"
    if n_samples is None:
        n_samples = sum(1 for _ in open(fam))
    bps = (n_samples + 3) // 4
    idx = np.asarray(snp_indices)
    lo, hi = idx.min(), idx.max()
    with open(panel + ".bed", "rb") as fh:
        fh.seek(3 + lo * bps)
        raw = np.frombuffer(fh.read((hi - lo + 1) * bps), dtype=np.uint8)
    raw = raw.reshape(hi - lo + 1, bps)
    codes = _LUT[raw].reshape(raw.shape[0], -1)[:, :n_samples][idx - lo]
    dos = np.full(codes.shape, np.nan)
    dos[codes == 0] = 2.0
    dos[codes == 2] = 1.0
    dos[codes == 3] = 0.0
    return dos


def ld_matrix(dos, maf_min=0.01, missing_max=0.05):
    """Pearson correlation matrix of dosages, with standard QC filters."""
    maf = np.nanmean(dos, axis=1) / 2.0
    maf = np.minimum(maf, 1 - maf)
    miss = np.isnan(dos).mean(axis=1)
    keep = (maf > maf_min) & (miss < missing_max)
    d = dos[keep]
    mu = np.nanmean(d, axis=1, keepdims=True)
    x = np.where(np.isnan(d), mu, d)
    x = x - x.mean(axis=1, keepdims=True)
    x = x / np.where(x.std(axis=1, keepdims=True) == 0, 1, x.std(axis=1, keepdims=True))
    return ((x @ x.T) / x.shape[1]).astype(np.float32), keep


# --------------------------------------------------------------------------- numerics
def logsumexp(x):
    m = np.max(x)
    if not np.isfinite(m):
        return m
    return m + np.log(np.exp(x - m).sum())


def logdiff(a, b):
    """log(exp(a) - exp(b)) for a >= b."""
    d = b - a
    return -np.inf if d >= 0 else a + np.log1p(-np.exp(d))


# --------------------------------------------------------------------------- io
def read_eqtl_region(path, chrom, start, end, cols=None):
    """Read a GTEx/eQTL-Catalogue region through tabix (requires pysam)."""
    import pysam
    default = ["variant", "r2", "pvalue", "mt_obj", "mt_id", "maf", "gene_id",
               "median_tpm", "beta", "se", "an", "ac", "chrom", "position",
               "ref", "alt", "type", "rsid"]
    tf = pysam.TabixFile(path, index=path + ".tbi")
    rows = [r.split("\t") for r in tf.fetch(str(chrom), int(start), int(end))]
    tf.close()
    d = pd.DataFrame(rows, columns=cols or default)
    for c in ["pvalue", "maf", "beta", "se", "position"]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d["ref"] = d.ref.str.upper()
    d["alt"] = d.alt.str.upper()
    return d


def read_gwas_region(clean_gz, chrom, start, end):
    """Read the project's clean GWAS format for one region."""
    d = pd.read_csv(clean_gz, sep="\t",
                    usecols=["SNP", "CHR", "POS", "A1", "A2", "BETA", "SE", "P", "FRQ", "N"],
                    dtype={"SNP": str, "A1": str, "A2": str})
    d = d[(d.CHR == int(chrom)) & (d.POS >= start) & (d.POS <= end)]
    return d.drop_duplicates("SNP")
