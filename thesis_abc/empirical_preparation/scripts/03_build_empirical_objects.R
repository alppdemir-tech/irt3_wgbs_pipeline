#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))


read_unit_table <- function(path) {

  if (grepl("\\.gz$", path)) {
    x <- fread(
      cmd = paste("gzip -dc", shQuote(path)),
      na.strings = "NA"
    )
  } else {
    x <- fread(
      path,
      na.strings = "NA"
    )
  }

  as.data.frame(x)
}


build_empirical_objects <- function(meanMeth_units) {

  pops <- c("benthic", "littoral", "riverine")

  if (!identical(names(meanMeth_units), pops)) {
    stop(
      "meanMeth_units must contain populations in order: ",
      paste(pops, collapse = ", ")
    )
  }

  # Structural information must be identical across populations.
  structural_cols <- setdiff(
    names(meanMeth_units$benthic),
    c("coverage", "meth_frac")
  )

  if (!identical(
    meanMeth_units$benthic[, structural_cols],
    meanMeth_units$littoral[, structural_cols]
  )) {
    stop("Benthic and littoral structural skeletons differ.")
  }

  if (!identical(
    meanMeth_units$benthic[, structural_cols],
    meanMeth_units$riverine[, structural_cols]
  )) {
    stop("Benthic and riverine structural skeletons differ.")
  }


  # ------------------------------------------------------------
  # Merge the three populations exactly as in methAnalysisReport.Rmd
  # ------------------------------------------------------------

  benthic <- meanMeth_units$benthic[
    ,
    c("unit_id", "chrom", "start", "end", "type", "meth_frac")
  ]

  littoral <- meanMeth_units$littoral[
    ,
    c("unit_id", "meth_frac")
  ]

  riverine <- meanMeth_units$riverine[
    ,
    c("unit_id", "meth_frac")
  ]

  names(benthic)[names(benthic) == "meth_frac"] <- "meth_benthic"
  names(littoral)[names(littoral) == "meth_frac"] <- "meth_littoral"
  names(riverine)[names(riverine) == "meth_frac"] <- "meth_riverine"

  merged <- merge(
    benthic,
    littoral,
    by = "unit_id"
  )

  merged <- merge(
    merged,
    riverine,
    by = "unit_id"
  )

  merged <- merged[
    order(merged$chrom, merged$start),
  ]


  # ------------------------------------------------------------
  # Keep only units observed in all three populations
  # ------------------------------------------------------------

  merged <- merged[
    !is.na(merged$meth_benthic) &
    !is.na(merged$meth_littoral) &
    !is.na(merged$meth_riverine),
  ]


  # ------------------------------------------------------------
  # Categorize methylation
  # Sarah:
  #   < 0.2       -> 0
  #   [0.2, 0.8]  -> 0.5
  #   > 0.8       -> 1
  # ------------------------------------------------------------

  merged$methCat_benthic <- ifelse(
    merged$meth_benthic < 0.2,
    0,
    ifelse(merged$meth_benthic > 0.8, 1, 0.5)
  )

  merged$methCat_littoral <- ifelse(
    merged$meth_littoral < 0.2,
    0,
    ifelse(merged$meth_littoral > 0.8, 1, 0.5)
  )

  merged$methCat_riverine <- ifelse(
    merged$meth_riverine < 0.2,
    0,
    ifelse(merged$meth_riverine > 0.8, 1, 0.5)
  )


  # ------------------------------------------------------------
  # Region ID
  # ------------------------------------------------------------

  merged$region_id <- sub(
    "_[0-9]+$",
    "",
    merged$unit_id
  )

  merged <- merged[
    order(merged$chrom, merged$start),
  ]


  # ------------------------------------------------------------
  # Build data_list and genome_spatial_str
  # ------------------------------------------------------------

  n_regions <- length(unique(merged$region_id))

  chrom_col <- character(n_regions)
  region_id_col <- character(n_regions)
  n_col <- integer(n_regions)
  globalState_col <- character(n_regions)

  row_i <- 1L

  chroms <- unique(merged$chrom)

  data_list <- setNames(
    vector("list", length(chroms)),
    chroms
  )

  for (chr_n in seq_along(chroms)) {

    chrom <- chroms[chr_n]

    data_list[[chrom]] <- setNames(
      vector("list", length(pops)),
      pops
    )

    merged_chrom <- merged[
      merged$chrom == chrom,
    ]

    regions_chrom <- unique(
      merged_chrom$region_id
    )

    for (pop_n in seq_along(pops)) {

      pop <- pops[pop_n]

      data_list[[chrom]][[pop]] <- setNames(
        vector("list", length(regions_chrom)),
        regions_chrom
      )

      for (reg_n in seq_along(regions_chrom)) {

        region <- regions_chrom[reg_n]

        merged_region <- merged_chrom[
          merged_chrom$region_id == region,
        ]

        col_name <- paste0(
          "methCat_",
          pop
        )

        data_list[[chrom]][[pop]][[region]] <-
          merged_region[[col_name]]

        if (pop_n == 1L) {

          region_types <- unique(
            merged_region$type
          )

          if (length(region_types) != 1L) {
            stop(
              "Region ",
              region,
              " contains multiple structure types."
            )
          }

          chrom_col[row_i] <- chrom
          region_id_col[row_i] <- region
          n_col[row_i] <- nrow(merged_region)

          globalState_col[row_i] <- ifelse(
            region_types == "archipelago",
            "U",
            "M"
          )

          row_i <- row_i + 1L
        }
      }
    }
  }

  genome_spatial_str <- data.frame(
    chrom = chrom_col,
    region_id = region_id_col,
    n = n_col,
    globalState = globalState_col,
    stringsAsFactors = FALSE
  )


  # ------------------------------------------------------------
  # Build sumstats_data
  # ------------------------------------------------------------

  sumstats_data <- list()

  for (chrom in names(data_list)) {

    for (pop in names(data_list[[chrom]])) {

      if (!pop %in% names(sumstats_data)) {
        sumstats_data[[pop]] <- list()
      }

      for (str in names(data_list[[chrom]][[pop]])) {

        str_chrom <- paste0(
          str,
          "*",
          chrom
        )

        sumstats_data[[pop]][[str_chrom]] <-
          data_list[[chrom]][[pop]][[str]]
      }
    }
  }


  # ------------------------------------------------------------
  # Internal consistency checks
  # ------------------------------------------------------------

  if (sum(genome_spatial_str$n) != nrow(merged)) {
    stop(
      "sum(genome_spatial_str$n) does not equal nrow(merged)."
    )
  }

  if (!identical(
    names(sumstats_data),
    pops
  )) {
    stop("Unexpected population order in sumstats_data.")
  }


  list(
    merged = merged,
    data_list = data_list,
    sumstats_data = sumstats_data,
    genome_spatial_str = genome_spatial_str
  )
}


# ------------------------------------------------------------
# Command-line interface
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# When sourced for validation, do not execute the CLI.
if (length(args) > 0) {

  get_arg <- function(flag) {

    i <- match(flag, args)

    if (is.na(i) || i == length(args)) {
      stop("Missing required argument: ", flag)
    }

    args[i + 1L]
  }


  benthic_file <- get_arg("--benthic")
  littoral_file <- get_arg("--littoral")
  riverine_file <- get_arg("--riverine")
  output_file <- get_arg("--output")


  cat("Loading unit-level methylation tables...\n")

  meanMeth_units <- list(
    benthic = read_unit_table(benthic_file),
    littoral = read_unit_table(littoral_file),
    riverine = read_unit_table(riverine_file)
  )


  cat(
    "Input rows per population:",
    nrow(meanMeth_units$benthic),
    "\n"
  )


  result <- build_empirical_objects(
    meanMeth_units
  )

  merged <- result$merged
  data_list <- result$data_list
  sumstats_data <- result$sumstats_data
  genome_spatial_str <- result$genome_spatial_str


  cat("\n=== Final empirical geometry ===\n")

  cat(
    "Common units:",
    format(nrow(merged), big.mark = ","),
    "\n"
  )

  cat(
    "Regions:",
    format(nrow(genome_spatial_str), big.mark = ","),
    "\n"
  )

  cat("\nGlobal states:\n")
  print(
    table(genome_spatial_str$globalState)
  )

  cat(
    "\nSum of n:",
    format(sum(genome_spatial_str$n), big.mark = ","),
    "\n"
  )

  cat(
    "Min n:",
    min(genome_spatial_str$n),
    "\n"
  )

  cat(
    "Max n:",
    max(genome_spatial_str$n),
    "\n"
  )


  cat("\nMethylation categories:\n")

  for (pop in c(
    "benthic",
    "littoral",
    "riverine"
  )) {

    cat("\n", pop, "\n", sep = "")

    print(
      table(
        merged[[paste0("methCat_", pop)]]
      )
    )
  }


  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  save(
    merged,
    data_list,
    sumstats_data,
    genome_spatial_str,
    file = output_file
  )

  cat(
    "\nSaved:",
    output_file,
    "\n"
  )
}
