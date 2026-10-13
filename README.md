# PGL: gene lists and sampled human brain regions

This repository contains the **R scientific analysis used by PGL Online** and compact catalogs of the genes and brain regions represented in its six-donor atlas. It is intended to let researchers inspect the method and query the atlas coverage, including with an AI assistant. The analysis code is in [`R/PGL.R`](R/PGL.R); the catalogs are in [`catalogs/`](catalogs/README.md).

## What the analysis does

1. Read a CSV with one exact gene symbol per row and no header. Symbols absent from the atlas do not contribute probes.
2. Map symbols to microarray probes using the Donor1 probe annotation. Label those probes as the target group and all other probes as the comparison group.
3. For each of six donor brains, compare target and comparison probe expression in each sampled region using a Wilcoxon test. Apply Bonferroni correction and select results with adjusted p-value ≤ 0.05.
4. Attach anatomical region names, keep up to ten selected regions per donor, and combine donor rankings. The exported `finalranking` uses the donor `BrainCompleteness` weights. `brainFreq` counts the donors in which each structure was selected.

The output CSV has `structurename`, `finalranking`, and `brainFreq`. These are **relative region rankings**, not gene-level effect sizes or evidence of causality. Read [`KNOWN_LIMITATIONS.md`](KNOWN_LIMITATIONS.md) before interpreting results: this reference copy preserves known behavior of the current implementation.

The R expressions in [`R/PGL.R`](R/PGL.R) were checked against the deployed `pgl_core.R` and are identical. Formatting and explanatory comments were added for readability; the scientific calculations were not changed. [`R/run_job.R`](R/run_job.R) is a small command-line entry point for this repository.

## Atlas coverage

The same **29,130 gene symbols and 58,692 reference probes** are available in all six donor expression matrices. This was verified from the exact `HumanBrainAtlas.rda` used by PGL Online; the probe annotations also match across donors. The sampled regions differ:

| Donor | Samples (`well_id`) | Distinct annotated regions |
| --- | ---: | ---: |
| Donor1 | 946 | 327 |
| Donor2 | 893 | 348 |
| Donor3 | 363 | 158 |
| Donor4 | 529 | 185 |
| Donor5 | 470 | 187 |
| Donor6 | 501 | 200 |

Browse the [combined region-presence table](catalogs/regions_presence.csv) to see all **414 distinct sampled structures** and which donor brains contain each one. [`catalogs/README.md`](catalogs/README.md) explains the files, columns, and example AI queries. **Sampled** means expression data is present for the region; it does not imply that PGL found the region significant for any particular gene list.

The same CSV catalogs are available as [direct static files](https://cognify-science.github.io/PGL/) for browsers and AI tools that cannot read GitHub's interactive CSV viewer. The static links serve the CSV data itself; use the `sample_count` column in each donor's region file for sampling questions.

## Run the R analysis

Install R packages `sqldf` and `dplyr`, obtain the matching `HumanBrainAtlas.rda`, and run:

```sh
Rscript R/run_job.R input.csv /path/to/HumanBrainAtlas.rda result.csv result.Rdata
```

The atlas is not included in this repository. Recreate the catalogs with:

```sh
Rscript scripts/export_catalogs.R /path/to/HumanBrainAtlas.rda catalogs
```

The catalog script uses base R. Its output and source-file fingerprint are recorded in [`catalogs/SOURCE.md`](catalogs/SOURCE.md).

## Source and attribution

The underlying microarray data and sample metadata come from the [Allen Human Brain Atlas](https://human.brain-map.org/static/download), which provides data for six brains. This repository distributes compact coverage lists and region annotations, not the full expression matrices or atlas file. Use and cite the source according to the [Allen Institute citation policy](https://alleninstitute.org/citation-policy/) and [terms of use](https://alleninstitute.org/terms-of-use/). The source file fingerprint in `catalogs/SOURCE.md` identifies the exact local atlas from which these lists were exported.
