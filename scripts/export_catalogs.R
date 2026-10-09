# Extract small, reproducible catalogs from the atlas used by PGL Online.
# Run: Rscript scripts/export_catalogs.R /path/to/HumanBrainAtlas.rda catalogs
# This script exports metadata only; it does not modify the atlas or run PGL.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Expected atlas path and output directory")
}

atlas_path <- args[[1L]]
output_dir <- args[[2L]]
if (!file.exists(atlas_path)) stop("Atlas file does not exist")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "genes"), showWarnings = FALSE)
dir.create(file.path(output_dir, "regions"), showWarnings = FALSE)

atlas_environment <- new.env(parent = emptyenv())
loaded_names <- load(atlas_path, envir = atlas_environment)
if (!"HumanBrainAtlas" %in% loaded_names) stop("Missing HumanBrainAtlas object")
atlas <- atlas_environment$HumanBrainAtlas
donor_names <- paste0("Donor", 1:6)
if (!all(donor_names %in% names(atlas))) stop("Expected six donors")

# GetROIs maps input symbols to Donor1 probes and joins those probe IDs to
# each donor's normalized expression matrix. Use this same definition of
# 'available' rather than listing all annotated genes without expression.
reference <- atlas$Donor1$Probe[, c("probe_id", "gene_symbol")]
reference$probe_id <- as.character(reference$probe_id)
reference$gene_symbol <- as.character(reference$gene_symbol)
reference <- unique(reference[!is.na(reference$gene_symbol) &
  nzchar(reference$gene_symbol) & !is.na(reference$probe_id), ])

gene_sets <- list()
region_catalogs <- list()
summary_rows <- list()
for (donor_name in donor_names) {
  donor <- atlas[[donor_name]]
  expression <- donor$Microarray$Standarized
  if (!"ProbeID" %in% names(expression)) stop("Missing ProbeID in ", donor_name)
  expressed_probes <- as.character(expression$ProbeID)
  matched <- reference[reference$probe_id %in% expressed_probes, ]
  if (nrow(matched) == 0L) stop("No matched probes in ", donor_name)
  probe_counts <- aggregate(probe_id ~ gene_symbol, matched,
    FUN = function(ids) length(unique(ids)))
  names(probe_counts) <- c("gene_symbol", "probe_count")
  probe_counts <- probe_counts[order(probe_counts$gene_symbol), ]
  rownames(probe_counts) <- NULL
  write.csv(probe_counts, file.path(output_dir, "genes", paste0(donor_name, ".csv")),
    row.names = FALSE, fileEncoding = "UTF-8")
  gene_sets[[donor_name]] <- probe_counts$gene_symbol

  # A sampled region is represented by a well_id column in expression and a
  # matching row in SampleAnnot. Keep individual samples as well as a summary.
  samples <- donor$SampleAnnot[, c("well_id", "structure_id",
    "structure_acronym", "structure_name")]
  if (anyNA(samples) || anyDuplicated(samples$well_id)) {
    stop("Incomplete or duplicate sample annotations in ", donor_name)
  }
  expression_wells <- names(expression)[names(expression) != "ProbeID"]
  if (!setequal(as.character(samples$well_id), expression_wells)) {
    stop("Sample annotations do not match expression columns in ", donor_name)
  }
  samples <- samples[order(samples$well_id), ]
  rownames(samples) <- NULL
  write.csv(samples, file.path(output_dir, "regions", paste0(donor_name, "_samples.csv")),
    row.names = FALSE, fileEncoding = "UTF-8")

  regions <- aggregate(well_id ~ structure_id + structure_acronym + structure_name,
    samples, FUN = function(ids) length(unique(ids)))
  names(regions)[names(regions) == "well_id"] <- "sample_count"
  regions <- regions[order(regions$structure_name, regions$structure_id), ]
  rownames(regions) <- NULL
  if (sum(regions$sample_count) != nrow(samples)) {
    stop("Region counts do not add up in ", donor_name)
  }
  write.csv(regions, file.path(output_dir, "regions", paste0(donor_name, ".csv")),
    row.names = FALSE, fileEncoding = "UTF-8")
  region_catalogs[[donor_name]] <- regions

  same_probe_annotation <- identical(
    as.character(donor$Probe$probe_id),
    as.character(atlas$Donor1$Probe$probe_id)
  ) && identical(
    as.character(donor$Probe$gene_symbol),
    as.character(atlas$Donor1$Probe$gene_symbol)
  )
  summary_rows[[donor_name]] <- data.frame(
    donor = donor_name, gene_count = nrow(probe_counts),
    matched_probe_count = length(unique(matched$probe_id)),
    sample_count = nrow(samples), region_count = nrow(regions),
    same_probe_annotation_as_donor1 = same_probe_annotation
  )
}

# One row per distinct sampled structure across the atlas. TRUE means that
# donor has at least one annotated sample for that structure.
region_fields <- c("structure_id", "structure_acronym", "structure_name")
all_regions <- unique(do.call(rbind, lapply(region_catalogs,
  function(regions) regions[, region_fields])))
all_regions <- all_regions[order(all_regions$structure_name,
  all_regions$structure_id), ]
rownames(all_regions) <- NULL
region_key <- function(regions) {
  paste(regions$structure_id, regions$structure_acronym,
    regions$structure_name, sep = "\u001f")
}
all_keys <- region_key(all_regions)
for (donor_name in donor_names) {
  all_regions[[donor_name]] <- all_keys %in% region_key(region_catalogs[[donor_name]])
}
all_regions$donor_count <- rowSums(all_regions[, donor_names])
write.csv(all_regions, file.path(output_dir, "regions_presence.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")

common_genes <- sort(Reduce(intersect, gene_sets))
write.csv(data.frame(gene_symbol = common_genes),
  file.path(output_dir, "genes_common_all_donors.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(do.call(rbind, summary_rows), file.path(output_dir, "summary.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")
cat("Exported", length(common_genes), "genes common to all donors\n")
cat("Exported", nrow(all_regions), "distinct sampled structures\n")
print(do.call(rbind, summary_rows), row.names = FALSE)
