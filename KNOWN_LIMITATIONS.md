# Known limitations of the current scientific code

This repository presents the R analysis currently used by PGL Online. Formatting and explanatory comments do not change its parsed expressions. The points below describe behavior that should be considered when interpreting results; they are not silently corrected in this reference copy.

- In `WilcoxonFS`, the comparison label is created as `target` but the code attempts to remove `Target`. The different capitalization means the label remains in the table of tested columns. This needs a separate scientific correction and validation against reference results.
- A donor with no significant regions can trigger a downstream error when the analysis attaches donor information or constructs the cross-donor summary. An empty result should therefore not yet be interpreted as evidence that no regions are associated with the input list.
- Symbols not matched to the atlas probe table contribute no probes. Matching is exact; gene aliases and case variants are not resolved automatically.
- The exported `finalranking` is a relative ranking assembled from donor-level selections and donor weights. It is not an effect size or a probability of biological association.

These limitations concern interpretation of the scientific analysis. The catalogs are extracted directly from the atlas and do not run the statistical method.
