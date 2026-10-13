# Atlas catalogs

These CSV files are derived from the six-donor atlas used by the R code in this repository. They describe **coverage**, not results of a PGL analysis. All CSVs are UTF-8 and include a header row.

For direct CSV downloads without GitHub's interactive viewer, use the [static catalog links](https://cognify-science.github.io/PGL/).

## Genes

- [`genes/Donor1.csv`](genes/Donor1.csv) through [`genes/Donor6.csv`](genes/Donor6.csv): one row per symbol available for analysis in that donor. `gene_symbol` is the exact case-sensitive symbol used by the R code; `probe_count` is the number of distinct Donor1 reference probe IDs that occur in that donor's normalized expression matrix. All six current files contain the same 29,130 symbols and 58,692 probes.
- [`genes_common_all_donors.csv`](genes_common_all_donors.csv): one `gene_symbol` per row in the intersection of the six lists. In this atlas, the intersection contains all 29,130 symbols.

"Available" here means an input symbol maps to at least one reference probe that occurs in the donor's expression matrix. It does not guarantee that the symbol is unique, passes a statistical test, or produces a region in a particular analysis. Gene aliases are not expanded.

## Brain regions

- [`regions_presence.csv`](regions_presence.csv): one row per distinct sampled structure across all six brains. It gives `structure_id`, `structure_acronym`, `structure_name`, six `Donor1`–`Donor6` columns (`TRUE` if sampled in that donor), and `donor_count`.
- [`regions/Donor1_samples.csv`](regions/Donor1_samples.csv) through [`regions/Donor6_samples.csv`](regions/Donor6_samples.csv): one row per sampled `well_id`, with its `structure_id`, `structure_acronym`, and `structure_name`.
- [`regions/Donor1.csv`](regions/Donor1.csv) through [`regions/Donor6.csv`](regions/Donor6.csv): one row per distinct annotated structure in a donor, with `sample_count` giving the number of sampled wells assigned to it.
- [`summary.csv`](summary.csv): donor-level gene, probe, sample, and region counts. `same_probe_annotation_as_donor1` checks both probe IDs and gene symbols against Donor1's annotation.

The same named structure may have several sampled wells. Region names and identifiers come from `SampleAnnot`. A region's presence in these files means it was **sampled**, not that it was statistically selected by PGL. The catalog covers the sampled structures in this atlas, not every anatomical region that exists in a human brain.

## Questions to ask with an AI assistant

- "Is `APOE` listed for all six donors? Show the matching rows and probe counts."
- "Which sampled structures contain 'hippocampus' in Donor3? Distinguish unique structures from individual wells."
- "Which donors contain this `structure_name`? Use `regions_presence.csv` and report its `structure_id`."
- "Compare the number of sampled wells for a named structure across the six donor files."

Ask the assistant to cite the specific CSV files and rows it used. For a biological interpretation of PGL results, also provide the input gene list and the result CSV; coverage files alone cannot establish significance.
