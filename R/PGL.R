# PGL Online scientific analysis. This is the R code used by the current service.
# Formatting and comments have been added without changing parsed expressions.

#' Find brain regions associated with an input list of gene symbols.
#' @param x Path to a one-column CSV, without a header.
#' @param atlas_path Path to HumanBrainAtlas.rda (not distributed here).
#' @param output_csv Destination for the ranked regions CSV.
#' @param output_rdata Destination for the complete Result object.
#' @return Result list; the CSV has structurename, finalranking, brainFreq.
GetROIs <- function(x, atlas_path, output_csv, output_rdata) {
  library(sqldf)
  library(dplyr)

  tryCatch(
    {
      FSFunction <- "W"
      print("FAST PGL")

      # Input symbols are matched exactly against the atlas probe annotation.
      GenesList <- read.csv(x, header = FALSE)

      print("Loading Human Brain Atlas")
      load(file = atlas_path, envir = .GlobalEnv)

      print("GetROIs Started")
      set.seed(1)

      Result <- list()

      # Donor1 supplies the probe-to-gene mapping used for all six donors.
      # Multiple probes may represent the same gene symbol.
      SelectedProbes <- HumanBrainAtlas$Donor1$Probe[HumanBrainAtlas$Donor1$Probe$gene_symbol %in% GenesList[, 1], ]
      `%ni%` <- Negate(`%in%`)

      # Mark probes matched to the submitted genes as target (1); other probes
      # become the comparison group (0).
      Result$ProbeIDTarget <- as.data.frame(
        rbind(
          cbind(HumanBrainAtlas$Donor1$Microarray$Raw[HumanBrainAtlas$Donor1$Microarray$Raw[, 1] %in% SelectedProbes[, 1], ]$ProbeID, 1),
          (cbind(HumanBrainAtlas$Donor1$Microarray$Raw[HumanBrainAtlas$Donor1$Microarray$Raw[, 1] %ni% SelectedProbes[, 1], ]$ProbeID, 0))
        )
      )

      colnames(Result$ProbeIDTarget) <- c("ProbeID", "Target")

      Result$FoundGenes <- as.list(unique(SelectedProbes$gene_symbol))
      print(paste0("Number of Found Genes:", length(Result$FoundGenes)))

      Result$OriginalList <- as.list(GenesList)

      # Test each of the six sampled brains separately.
      for (i in 1:6) {
        Donor <- paste("Donor", i, sep = "")

        print(Donor)

        # Join probe labels to this donor's normalized expression matrix.
        # Each remaining expression column corresponds to one sampled well_id.
        data <- merge(Result$ProbeIDTarget, HumanBrainAtlas[[Donor]]$Microarray$Standarized, by = "ProbeID")
        target <- data$Target
        data$ProbeID <- NULL
        data$Target <- NULL
        colnames(data) <- c(HumanBrainAtlas[[Donor]]$SampleAnnot$well_id)

        # Compare target and non-target probes in each sampled region.
        Result$WilcoxonFS[[Donor]] <- WilcoxonFS(data, target)
        print(paste0("# Regiones Significativas: ", nrow(Result$WilcoxonFS[[Donor]])))
        if (ncol(Result$WilcoxonFS[[Donor]]) > 0) {
          colnames(Result$WilcoxonFS[[Donor]])[1] <- "well_id"
          colnames(Result$WilcoxonFS[[Donor]])[2] <- "p-value"
          colnames(Result$WilcoxonFS[[Donor]])[3] <- "corrected.p-value"
          colnames(Result$WilcoxonFS[[Donor]])[4] <- "conf.int.inf"
          colnames(Result$WilcoxonFS[[Donor]])[5] <- "conf.int.sup"
          colnames(Result$WilcoxonFS[[Donor]])[6] <- "difference"
          colnames(Result$WilcoxonFS[[Donor]])[7] <- "W"
          Result$WilcoxonFS[[Donor]] <- merge(Result$WilcoxonFS[[Donor]], HumanBrainAtlas[[Donor]]$SampleAnnot, by = "well_id")
          Result$WilcoxonFS[[Donor]]$Donor <- i
        }
      }

      # Combine donor-specific significant regions into cross-donor rankings.
      Result <- summarizePGL(Result)
      Result_allregions <- c(
        unique(Result$WilcoxonFS$Donor1$structure_name),
        unique(Result$WilcoxonFS$Donor2$structure_name),
        unique(Result$WilcoxonFS$Donor3$structure_name),
        unique(Result$WilcoxonFS$Donor4$structure_name),
        unique(Result$WilcoxonFS$Donor5$structure_name),
        unique(Result$WilcoxonFS$Donor6$structure_name)
      )

      # brainFreq counts donors in which a structure was selected; each donor
      # contributes that structure at most once.
      Result_allregions_table <- as.data.frame(table(Result_allregions))
      colnames(Result_allregions_table) <- c("structurename", "brainFreq")
      Result$summarized$Merged$final1 <- merge(Result$summarized$Merged$final2, Result_allregions_table, by = "structurename")
      Result$summarized$Merged$final2 <- merge(Result$summarized$Merged$final2, Result_allregions_table, by = "structurename")
      Result$summarized$Merged$final3 <- merge(Result$summarized$Merged$final2, Result_allregions_table, by = "structurename")

      Result$summarized$Merged$final1 <- Result$summarized$Merged$final1[order(Result$summarized$Merged$final1$finalranking), ]
      Result$summarized$Merged$final2 <- Result$summarized$Merged$final2[order(Result$summarized$Merged$final2$finalranking), ]
      Result$summarized$Merged$final3 <- Result$summarized$Merged$final3[order(Result$summarized$Merged$final3$finalranking), ]

      # final2 weights donor ranks by BrainCompleteness and is the published CSV.
      Result_to_Save <- Result$summarized$Merged$final2
      Result_to_Save <- Result_to_Save[order(-Result_to_Save$finalranking), ]
      write.csv(Result_to_Save, file = output_csv, row.names = FALSE)
      save(Result, file = output_rdata)

      return(Result)
    },
    error = function(cond) {
      stop(cond)
    }
  )
}

#' Compare expression for target versus other probes, one donor at a time.
#' @param x Probe-by-region normalized expression matrix.
#' @param y Binary probe label: 1 for a submitted gene, otherwise 0.
#' @return Columns with Bonferroni-adjusted p-value <= 0.05.
WilcoxonFS <- function(x, y) {
  ncolX = ncol(x)
  data <- cbind(x, y)
  colnames(data)[ncolX + 1] <- "target"
  data.class0 <- subset(data, target == 0)
  data.class1 <- subset(data, target == 1)

  # Known limitation: the column above is named 'target' but this line
  # removes 'Target'. The difference is retained to match production.
  data.class0$Target <- NULL
  data.class1$Target <- NULL

  ncoldata <- ncol(data.class0)
  wilcoxtest <- list()
  wilcoxtest.pvalue <- 1:ncoldata
  wilcoxtest.confint1 <- 1:ncoldata
  wilcoxtest.confint2 <- 1:ncoldata
  wilcoxtest.statistic <- 1:ncoldata
  wilcoxtest.difference <- 1:ncoldata

  confinter <- FALSE

  for (i in 1:ncoldata) {
    cat("\r", paste("Area:", i))
    wilcoxtest[[i]] <- try(wilcox.test(data.class1[, i], data.class0[, i], conf.int = confinter, conf.level = 0.95))

    if (class(wilcoxtest[[i]]) == "try-error") {
      print("ERROR")
      wilcoxtest.pvalue[i] <- 0
      wilcoxtest.confint1[i] <- 0
      wilcoxtest.confint2[i] <- 0
      wilcoxtest.difference[i] <- 0
      wilcoxtest.statistic[i] <- 0
    }
    if (confinter == FALSE) {
      wilcoxtest.pvalue[i] <- wilcoxtest[[i]]$p.value
      wilcoxtest.confint1[i] <- 0
      wilcoxtest.confint2[i] <- 0
      wilcoxtest.difference[i] <- 0
      wilcoxtest.statistic[i] <- 0
    } else {
      wilcoxtest.pvalue[i] <- wilcoxtest[[i]]$p.value
      wilcoxtest.confint1[i] <- wilcoxtest[[i]]$conf.int[1]
      wilcoxtest.confint2[i] <- wilcoxtest[[i]]$conf.int[2]
      wilcoxtest.difference[i] <- wilcoxtest[[i]]$estimate
      wilcoxtest.statistic[i] <- wilcoxtest[[i]]$statistic
    }
  }

  # Adjust the collection of tests for this donor using Bonferroni.
  wilcoxtest.bonferroni <- p.adjust(wilcoxtest.pvalue, "bonferroni")

  features.selected <- data.frame()

  Header <- colnames(x)

  for (i in 1:ncoldata) {
    if (wilcoxtest.bonferroni[i] <= 0.05) {
      features.selected <- rbind(features.selected, cbind(as.character(Header[i]), wilcoxtest.pvalue[i], wilcoxtest.bonferroni[i], wilcoxtest.confint1[i], wilcoxtest.confint2[i], wilcoxtest.difference[i], wilcoxtest.statistic[i]))

      features.selected[, 1] <- as.character(features.selected[, 1])
    }
  }

  return(features.selected)
}

#' Summarize significant regions across six donors.
#' The current workflow uses the Wilcoxon branch and exports final2,
#' whose donor scores are weighted by BrainCompleteness.
#' @param Result Donor-specific statistical results.
#' @return Result with donor summaries and merged rankings.
summarizePGL <- function(Result) {
  print("Start Summarizing")
  Result$summarized <- list()
  ndonors <- c(1, 2, 3, 4, 5, 6)
  methods <- vector()
  if (!is.null(Result$WilcoxonFS)) {
    methods <- c(methods, "WilcoxonFS")
  }
  if (!is.null(Result$RandomForestFS)) {
    methods <- c(methods, "RandomForestFS")
  }
  if (!is.null(Result$SVMFS)) {
    methods <- c(methods, "SVMFS")
  }
  if (!is.null(Result$rKNNFS)) {
    methods <- c(methods, "rKNNFS")
  }

  donorlist <- list()
  donorcount <- 0
  for (i in 1:length(ndonors)) {
    Donor <- paste("Donor", ndonors[i], sep = "")
    Result$summarized[[Donor]] <- list()

    for (j in 1:length(methods)) {
      donorcount <- donorcount + 1
      method <- methods[j]
      temp <- Result[[method]][[Donor]]
      if (!is.null(temp)) {
        if (method == "WilcoxonFS") {
          temp <- Result[[method]][[Donor]]
          colnames(temp)[3] <- "correctedpvalue"
          temp[, 3] <- as.numeric(as.character(temp[, 3]))
          colname <- "correctedpvalue"
          Result$summarized[[Donor]][[method]] <- sqldf(paste("SELECT structure_name, COUNT(structure_name), MIN(", colname, ") FROM temp GROUP BY structure_name ORDER BY MIN(", colname, ")", sep = ""))
        }
        if (method == "RandomForestFS") {
          temp <- Result[[method]][[Donor]]$rfvarimp
          colname <- "rfimportance"
          temp[[colname]] <- as.numeric(as.character(temp[[colname]]))
          Result$summarized[[Donor]][[method]] <- sqldf(paste("SELECT structure_name, COUNT(structure_name), MAX(", colname, ") FROM temp GROUP BY structure_name ORDER BY MAX(", colname, ") DESC", sep = ""))
        }
        if (method == "SVMFS") {
          temp <- Result[[method]][[Donor]]$svmvarimp
          colname <- "svmimportance"
          temp[[colname]] <- as.numeric(as.character(temp[[colname]]))
          Result$summarized[[Donor]][[method]] <- sqldf(paste("SELECT structure_name, COUNT(structure_name), MAX(", colname, ") FROM temp GROUP BY structure_name ORDER BY MAX(", colname, ") DESC", sep = ""))
        }
        if (method == "rKNNFS") {
          temp <- Result[[method]][[Donor]]$knnvarimp
          colname <- "knnimportance"
          temp[[colname]] <- as.numeric(as.character(temp[[colname]]))
          Result$summarized[[Donor]][[method]] <- sqldf(paste("SELECT structure_name, COUNT(structure_name), MAX(", colname, ") FROM temp GROUP BY structure_name ORDER BY MAX(", colname, ") DESC", sep = ""))
        }
        colnames(Result$summarized[[Donor]][[method]]) <- c("structurename", "structurenamecount", paste("min", colname, sep = ""))

        Result$summarized[[Donor]][[method]] <- Result$summarized[[Donor]][[method]][1:10, ]
        Result$summarized[[Donor]][[method]] <- Result$summarized[[Donor]][[method]][complete.cases(Result$summarized[[Donor]][[method]]), ]
        Result$summarized[[Donor]][[method]]$ranking <- nrow(Result$summarized[[Donor]][[method]]):1 / 10
        donorlist[[donorcount]] <- as.data.frame(cbind(as.character(Result$summarized[[Donor]][[method]]$structurename), Result$summarized[[Donor]][[method]]$ranking))
        colnames(donorlist[[donorcount]]) <- c("structurename", paste("Donor", i, method, "ranking", sep = ""))
      }
    }
  }

  Result$summarized$Merged <- Reduce(function(x, y) merge(x = x, y = y, by = "structurename", all = TRUE, no.dups = FALSE), donorlist)

  if (ncol(Result$summarized$Merged) == 2) {
    Result$summarized$Merged[, 2] <- as.numeric(Result$summarized$Merged[, 2])
  } else {
    Result$summarized$Merged[, 2:ncol(Result$summarized$Merged)] <- Result$summarized$Merged[, 2:ncol(Result$summarized$Merged)] %>%
      mutate_all(as.numeric)
  }

  Result$summarized$Merged[is.na(Result$summarized$Merged)] <- 0

  Donor1Merged <- Result$summarized$Merged[grepl("Donor1", names(Result$summarized$Merged))]
  Donor2Merged <- Result$summarized$Merged[grepl("Donor2", names(Result$summarized$Merged))]
  Donor3Merged <- Result$summarized$Merged[grepl("Donor3", names(Result$summarized$Merged))]
  Donor4Merged <- Result$summarized$Merged[grepl("Donor4", names(Result$summarized$Merged))]
  Donor5Merged <- Result$summarized$Merged[grepl("Donor5", names(Result$summarized$Merged))]
  Donor6Merged <- Result$summarized$Merged[grepl("Donor6", names(Result$summarized$Merged))]

  Result$summarized$Merged$Mean1 <- rowMeans(Donor1Merged)
  Result$summarized$Merged$Mean2 <- rowMeans(Donor2Merged)
  Result$summarized$Merged$Mean3 <- rowMeans(Donor3Merged)
  Result$summarized$Merged$Mean4 <- rowMeans(Donor4Merged)
  Result$summarized$Merged$Mean5 <- rowMeans(Donor5Merged)
  Result$summarized$Merged$Mean6 <- rowMeans(Donor6Merged)

  meanbind <- as.data.frame(cbind(Result$summarized$Merged$Mean1, Result$summarized$Merged$Mean2, Result$summarized$Merged$Mean3, Result$summarized$Merged$Mean4, Result$summarized$Merged$Mean5, Result$summarized$Merged$Mean6))
  meanbind <- meanbind[colSums(!is.na(meanbind)) > 0]

  Result$summarized$Merged$final1 <- cbind(as.character(Result$summarized$Merged$structurename), rowMeans(meanbind))
  Result$summarized$Merged$final1 <- as.data.frame(Result$summarized$Merged$final1)
  colnames(Result$summarized$Merged$final1) <- c("structurename", "finalranking")
  Result$summarized$Merged$final1$finalranking <- as.numeric(as.character(Result$summarized$Merged$final1$finalranking))
  Result$summarized$Merged$final1 <- Result$summarized$Merged$final1[order(-Result$summarized$Merged$final1$finalranking), c(1, 2)]

  # final2 applies each donor's BrainCompleteness weight.
  Result$summarized$Merged$Mean1 <- rowMeans(Donor1Merged) * HumanBrainAtlas$Donor1$summary$BrainCompleteness
  Result$summarized$Merged$Mean2 <- rowMeans(Donor2Merged) * HumanBrainAtlas$Donor2$summary$BrainCompleteness
  Result$summarized$Merged$Mean3 <- rowMeans(Donor3Merged) * HumanBrainAtlas$Donor3$summary$BrainCompleteness
  Result$summarized$Merged$Mean4 <- rowMeans(Donor4Merged) * HumanBrainAtlas$Donor4$summary$BrainCompleteness
  Result$summarized$Merged$Mean5 <- rowMeans(Donor5Merged) * HumanBrainAtlas$Donor5$summary$BrainCompleteness
  Result$summarized$Merged$Mean6 <- rowMeans(Donor6Merged) * HumanBrainAtlas$Donor6$summary$BrainCompleteness

  Result$summarized$Merged$final2 <- cbind(as.character(Result$summarized$Merged$structurename), rowMeans(cbind(Result$summarized$Merged$Mean1, Result$summarized$Merged$Mean2, Result$summarized$Merged$Mean3, Result$summarized$Merged$Mean4, Result$summarized$Merged$Mean5, Result$summarized$Merged$Mean6)))
  Result$summarized$Merged$final2 <- as.data.frame(Result$summarized$Merged$final2)
  colnames(Result$summarized$Merged$final2) <- c("structurename", "finalranking")
  Result$summarized$Merged$final2$finalranking <- as.numeric(as.character(Result$summarized$Merged$final2$finalranking))
  Result$summarized$Merged$final2 <- Result$summarized$Merged$final2[order(-Result$summarized$Merged$final2$finalranking), c(1, 2)]

  # final3 additionally weights donor ranks by BrainStructureDensity.
  Result$summarized$Merged$Mean1 <- rowMeans(Donor1Merged) * HumanBrainAtlas$Donor1$summary$BrainCompleteness * HumanBrainAtlas$Donor1$summary$BrainStructureDensity
  Result$summarized$Merged$Mean2 <- rowMeans(Donor2Merged) * HumanBrainAtlas$Donor2$summary$BrainCompleteness * HumanBrainAtlas$Donor2$summary$BrainStructureDensity
  Result$summarized$Merged$Mean3 <- rowMeans(Donor3Merged) * HumanBrainAtlas$Donor3$summary$BrainCompleteness * HumanBrainAtlas$Donor3$summary$BrainStructureDensity
  Result$summarized$Merged$Mean4 <- rowMeans(Donor4Merged) * HumanBrainAtlas$Donor4$summary$BrainCompleteness * HumanBrainAtlas$Donor4$summary$BrainStructureDensity
  Result$summarized$Merged$Mean5 <- rowMeans(Donor5Merged) * HumanBrainAtlas$Donor5$summary$BrainCompleteness * HumanBrainAtlas$Donor5$summary$BrainStructureDensity
  Result$summarized$Merged$Mean6 <- rowMeans(Donor6Merged) * HumanBrainAtlas$Donor6$summary$BrainCompleteness * HumanBrainAtlas$Donor6$summary$BrainStructureDensity

  Result$summarized$Merged$final3 <- cbind(as.character(Result$summarized$Merged$structurename), rowMeans(cbind(Result$summarized$Merged$Mean1, Result$summarized$Merged$Mean2, Result$summarized$Merged$Mean3, Result$summarized$Merged$Mean4, Result$summarized$Merged$Mean5, Result$summarized$Merged$Mean6)))
  Result$summarized$Merged$final3 <- as.data.frame(Result$summarized$Merged$final3)
  colnames(Result$summarized$Merged$final3) <- c("structurename", "finalranking")
  Result$summarized$Merged$final3$finalranking <- as.numeric(as.character(Result$summarized$Merged$final3$finalranking))
  Result$summarized$Merged$final3 <- Result$summarized$Merged$final3[order(-Result$summarized$Merged$final3$finalranking), c(1, 2)]

  return(Result)
}
