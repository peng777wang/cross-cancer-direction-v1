# Cross-cancer direction and MTAG-borrowing explorer

A single self-contained HTML file that lets a reader query every locus in the
manuscript's cross-cancer direction map, see the per-cancer evidence behind each
locus, and see how much of every MTAG estimate is borrowed rather than observed.

## What is in it

| file | what it is |
|---|---|
| `pleiotropy_explorer.html` | the explorer. Open it in any browser; it works offline, has no dependencies and needs no server (340 KB) |
| `pleiotropy_resource.json` | the same data machine-readable: 273 loci with, for every cancer, the MTAG z, the cancer's own GWAS z and the borrowed fraction (plus the GRCh37 cytoband table) |
| `pleiotropy_resource.csv` | flat form, one row per locus x cancer (2,730 rows) |
| `locus_summary.csv` | one row per locus: class, DSI under MTAG and under the own GWAS, pair counts, median borrowing at antagonistic loci, whether the locus is supported by its own GWAS |
| `build_explorer.py`, `explorer_template.html` | how the HTML is produced from the JSON |
| `test_explorer.js` | 21 functional checks over the real page script (see below) |

## Using it

Open `pleiotropy_explorer.html`. You can

* search by rsID, gene, region (`15:78000000-79500000`) or cytogenetic band
  (`8q24`, `8q24.21`, `15q25`, `8q`);
* filter by direction class, by how many cancers a locus is significant in, and by
  chromosome;
* sort any column, and click a locus to see its per-cancer MTAG z against that
  cancer's own GWAS z, with the borrowed share of each estimate;
* download the current view as CSV.

Cytogenetic bands are resolved through the bundled UCSC hg19 table, not by
assuming the band number is a megabase: 8q24 spans 117.7-146.4 Mb in GRCh37.

## What it does and does not claim

The page states its own limits at the top, and they are not decoration:

* the direction classes are **inferred from MTAG**. In the reversed class no cancer
  reaches genome-wide significance in its own GWAS and the median borrowing is
  76-82%, so the reversal is a model prediction, not established biological
  antagonism;
* the **borrowed** column is `1 - |z_own| / |z_MTAG|`. It is not a defect of MTAG
  but the quantity that says which estimates are observations and which are
  predictions. It can be **negative**, which means that cancer's own GWAS
  estimate is larger than its MTAG estimate — MTAG shrinks as well as amplifies;
* the loci were selected for cross-cancer pleiotropy, so the excess of aligned loci
  is an **upper bound** on concordance. The reversed class shows no excess over the
  sign-flip null;
* z-scores are on the standardised scale and are not comparable across cancers with
  different sample sizes.

## Hosting it

The file is static, so any of these works and all are free:

* commit it to a GitHub repository and turn on Pages; the URL becomes
  `https://<user>.github.io/<repo>/pleiotropy_explorer.html`;
* drop it on any static host (Zenodo, OSF, an institutional page);
* or just ship it as supplementary material, since it opens from the filesystem.

Journals usually ask for a stable URL for a resource; the GitHub Pages route is
the least effort for that.

## How it was verified

No browser was available in the build environment, so the page's own script was
exercised against a minimal DOM stub covering the behaviour that can actually
break:

```bash
node test_explorer.js        # 21 checks, all passing
```

The checks cover the initial render, rsID / gene / region / cytoband search, the
class and k filters, numeric sorting, the row-click detail panel, reset, and the
CSV export. The cytoband expectations are hard-coded from the UCSC table rather
than read back from the page, so the test is independent of the code it tests.

Also checked: the HTML has balanced tags, contains no external references of any
kind (`http`, `img`, `script src`, `link`), and the embedded JSON parses.

What this does **not** verify: fonts, colours and layout in a real browser. Those
need one visual pass by eye; the markup is plain and dependency-free, so the risk
is low, but it should be looked at once before submission.
