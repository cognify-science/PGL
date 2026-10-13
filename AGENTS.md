# PGL data access and interpretation

Use the original CSV catalogs to answer atlas-coverage questions. This file contains instructions, not replacement data. Follow the user's requested scope.

## Choose the access method

1. If the repository is checked out or connected through a repository tool, read the CSV files directly. Prefer local file reads or repository APIs to web previews.
2. If only web browsing is available, use the explicit GitHub plain-view links below. Preserve `?plain=1`; the default GitHub CSV preview may expose only the header to a browser that does not execute JavaScript.
3. For tools that can fetch CSV responses, use the static CSV links below. Raw GitHub URLs can also work, but an individual tool may block that host.
4. Verify that you received data records, not just column names or GitHub's displayed line count. Do not treat a failed fetch as an empty dataset or zero samples. Try another listed access method before reporting failure.
5. If access still fails, report the exact URL opened, whether its query was preserved, which access methods failed, and whether any actual CSV record was received. Do not claim the repository is private, empty, or globally blocked based on one tool's failure. Request the relevant CSV attachment or a connected checkout when necessary.

GitHub plain views of the original region CSVs:

- [Donor1](https://github.com/cognify-science/PGL/blob/main/catalogs/regions/Donor1.csv?plain=1)
- [Donor2](https://github.com/cognify-science/PGL/blob/main/catalogs/regions/Donor2.csv?plain=1)
- [Donor3](https://github.com/cognify-science/PGL/blob/main/catalogs/regions/Donor3.csv?plain=1)
- [Donor4](https://github.com/cognify-science/PGL/blob/main/catalogs/regions/Donor4.csv?plain=1)
- [Donor5](https://github.com/cognify-science/PGL/blob/main/catalogs/regions/Donor5.csv?plain=1)
- [Donor6](https://github.com/cognify-science/PGL/blob/main/catalogs/regions/Donor6.csv?plain=1)
- [Region presence](https://github.com/cognify-science/PGL/blob/main/catalogs/regions_presence.csv?plain=1)

Static CSV responses, published from the same repository branch:

- [Donor1 CSV](https://cognify-science.github.io/PGL/catalogs/regions/Donor1.csv)
- [Donor2 CSV](https://cognify-science.github.io/PGL/catalogs/regions/Donor2.csv)
- [Donor3 CSV](https://cognify-science.github.io/PGL/catalogs/regions/Donor3.csv)
- [Donor4 CSV](https://cognify-science.github.io/PGL/catalogs/regions/Donor4.csv)
- [Donor5 CSV](https://cognify-science.github.io/PGL/catalogs/regions/Donor5.csv)
- [Donor6 CSV](https://cognify-science.github.io/PGL/catalogs/regions/Donor6.csv)
- [Region presence CSV](https://cognify-science.github.io/PGL/catalogs/regions_presence.csv)

A web tool may retain an older page. If the user supplies a commit-specific URL, keep that revision when following relative links. Do not remove the revision or `?plain=1`. A new instruction file cannot invalidate a browser's cache or override its access policy.

## Interpret and join the CSVs

- `catalogs/regions/DonorN.csv` has one row per annotated structure in donor N. `sample_count` counts sampled wells (`well_id`), not donor brains. A missing structure means zero samples only after reading the complete donor catalog.
- `catalogs/regions/DonorN_samples.csv` has one row per sampled well. Use it to audit sample counts.
- `catalogs/regions_presence.csv` records which donors sampled each structure. `donor_count` counts donors with at least one sample.
- PGL output `brainFreq` counts donors in which PGL selected a structure for the particular gene list. It is distinct from both `sample_count` and atlas `donor_count`.
- Match a result's `structurename` to the catalogs' `structure_name` exactly. Preserve left/right and anatomical subdivisions. Check for duplicate or ambiguous keys. Do not silently use fuzzy matching, replace a missing count with zero after a failed join, or substitute another atlas or Excel table.
- When adding sampling information to a result, preserve the existing scientific result columns. State whether the new count is per donor or the sum across all six donors; keep unmatched rows explicitly identified.
- For queries such as "subiculum", search `structure_name` case-insensitively and report each matching hemisphere separately. Use counts from CSV records and cite the source files, matched names or IDs, and revision when known.
- Gene coverage files are `catalogs/genes/DonorN.csv` and `catalogs/genes_common_all_donors.csv`. Read these as CSV data rather than inferring availability from descriptions. Gene symbols are case-sensitive in the analysis; aliases are not expanded.
- Read `catalogs/README.md` for column definitions and `KNOWN_LIMITATIONS.md` before explaining scientific results. The full atlas is not required to answer coverage questions from the exported catalogs.

## Verify that instructions are loaded

Repository agents may automatically load supported instruction filenames. Ordinary web chats do not necessarily load this file merely because the user supplies the repository URL. The root README links to it for discovery; web-only users may need to explicitly request that it be read. Instructions guide retrieval but cannot grant network access.
