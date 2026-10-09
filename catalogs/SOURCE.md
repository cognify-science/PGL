# Catalog provenance

- Source: `HumanBrainAtlas.rda` used by the current PGL Online analysis.
- Local file size: **3,160,925,869 bytes**.
- SHA-256: `6f84875bd5275f6ea2a7e2f47780ab379bc9317c363933a6542f64389fe27f7d`.
- Verification: on **2026-10-09**, both the local source and the private production atlas object had this size and SHA-256. The private atlas itself is not distributed in this repository.
- Exporter: [`scripts/export_catalogs.R`](../scripts/export_catalogs.R), using the Donor1 probe annotation, each donor's normalized expression `ProbeID` values, and each donor's `SampleAnnot`.
- Upstream resource: [Allen Human Brain Atlas microarray downloads](https://human.brain-map.org/static/download). See the [Allen Institute citation policy](https://alleninstitute.org/citation-policy/) and [terms of use](https://alleninstitute.org/terms-of-use/).

These lists are a derived view of the exact R atlas used by PGL Online. The SHA-256 identifies the container file, but the repository does not assert that the container is an unmodified copy of an Allen download; it may include local preparation steps. The exported CSVs contain symbols and sample annotations, not expression measurements.
