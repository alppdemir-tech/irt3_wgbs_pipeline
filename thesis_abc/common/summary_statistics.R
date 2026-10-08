# Common summary-statistic functions for cichlid MethEvolSIM simulations.
#
# Cichlid-specific choices follow the existing cichlid workflow:
#   globalState "U" = archipelago/island
#   globalState "M" = interarchipelago/non-island
#   benthic + littoral form the single cherry (cherry_index = 1)
#   minN_sites = 5
#   shore_length = 2
#
# The returned 14 statistics follow the broader ABC statistic set used by Ferit,
# while MeanSiteFChange uses the cichlid-specific single-cherry extraction.

prepare_cichlid_sumstat_data <- function(sim_output, genome_spatial_str) {

  required_cols <- c("chrom", "region_id", "n", "globalState")
  missing_cols <- setdiff(required_cols, names(genome_spatial_str))

  if (length(missing_cols) > 0L) {
    stop(
      "genome_spatial_str is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!is.list(sim_output) || length(sim_output) == 0L) {
    stop("sim_output must be a non-empty chromosome-level list.")
  }

  chr_ids <- names(sim_output)

  if (is.null(chr_ids) || any(chr_ids == "")) {
    stop("All chromosomes in sim_output must be named.")
  }

  unknown_chr <- setdiff(chr_ids, unique(genome_spatial_str$chrom))
  if (length(unknown_chr) > 0L) {
    stop(
      "Chromosomes absent from genome_spatial_str: ",
      paste(unknown_chr, collapse = ", ")
    )
  }

  population_names <- names(sim_output[[1]])

  if (is.null(population_names) || any(population_names == "")) {
    stop("Population lists in sim_output must be named.")
  }

  data_restructured <- setNames(
    vector("list", length(population_names)),
    population_names
  )

  spatial_parts <- vector("list", length(chr_ids))

  for (chr_i in seq_along(chr_ids)) {

    chr <- chr_ids[[chr_i]]

    if (!identical(names(sim_output[[chr]]), population_names)) {
      stop("Population names/order differ between chromosomes: ", chr)
    }

    chr_spatial <- genome_spatial_str[
      genome_spatial_str$chrom == chr,
      ,
      drop = FALSE
    ]

    n_regions <- nrow(chr_spatial)

    if (n_regions == 0L) {
      stop("No spatial-structure rows found for chromosome: ", chr)
    }

    for (pop in population_names) {

      pop_regions <- sim_output[[chr]][[pop]]

      if (length(pop_regions) != n_regions) {
        stop(
          "Region count mismatch for ", chr, " / ", pop,
          ": expected ", n_regions,
          ", observed ", length(pop_regions)
        )
      }

      observed_lengths <- vapply(pop_regions, length, integer(1))

      if (!identical(observed_lengths, as.integer(chr_spatial$n))) {
        stop("Region site-count mismatch for ", chr, " / ", pop)
      }

      region_names <- paste0(chr_spatial$region_id, "*", chr)
      names(pop_regions) <- region_names

      data_restructured[[pop]] <- c(
        data_restructured[[pop]],
        pop_regions
      )
    }

    spatial_parts[[chr_i]] <- chr_spatial
  }

  spatial_used <- do.call(rbind, spatial_parts)
  rownames(spatial_used) <- NULL

  expected_names <- paste0(
    spatial_used$region_id,
    "*",
    spatial_used$chrom
  )

  for (pop in population_names) {
    if (!identical(names(data_restructured[[pop]]), expected_names)) {
      stop(
        "Restructured region order does not match genome_spatial_str for ",
        pop
      )
    }
  }

  island_indices <- which(spatial_used$globalState == "U")
  nonisland_indices <- which(spatial_used$globalState == "M")

  if (length(island_indices) == 0L) {
    stop("No U/archipelago regions are present in the simulated subset.")
  }

  if (length(nonisland_indices) == 0L) {
    stop("No M/interarchipelago regions are present in the simulated subset.")
  }

  list(
    data = data_restructured,
    spatial = spatial_used,
    island_indices = island_indices,
    nonisland_indices = nonisland_indices
  )
}


compute_cichlid_sumstats <- function(
    sim_output,
    tree_1scaled,
    genome_spatial_str,
    L,
    minN_sites = 5L,
    shore_length = 2L,
    cherry_index = 1L,
    p_value_threshold = 0.05) {

  if (!inherits(tree_1scaled, "phylo")) {
    stop("tree_1scaled must be an ape phylo object.")
  }

  if (!is.numeric(L) || length(L) != 1L || !is.finite(L) || L <= 0) {
    stop("L must be one positive finite number.")
  }

  prepared <- prepare_cichlid_sumstat_data(
    sim_output = sim_output,
    genome_spatial_str = genome_spatial_str
  )

  data <- prepared$data
  island_indices <- prepared$island_indices
  nonisland_indices <- prepared$nonisland_indices

  tip_names <- tree_1scaled$tip.label

  if (!setequal(names(data), tip_names)) {
    stop(
      "Population names in simulation output do not match tree tip labels."
    )
  }

  # Put populations in exactly the same order as the tree tips.
  data <- data[tip_names]

  sample_n <- length(tip_names)

  # Scale the biological tree by L, but do NOT add the artificial
  # root branch used during simulation.
  tree_no_root <- tree_1scaled
  tree_no_root$edge.length <- tree_no_root$edge.length * L
  tree_string <- ape::write.tree(tree_no_root)

  cherry_stats <- MethEvolSIM::MeanSiteFChange_cherry(
    data = data,
    categorized_data = TRUE,
    tree = tree_string,
    index_islands = island_indices,
    index_nonislands = nonisland_indices
  )

  if (cherry_index < 1L ||
      cherry_index > length(cherry_stats$island_meanFChange)) {
    stop("cherry_index is outside the available cherry range.")
  }

  raw_fitch <- MethEvolSIM::computeFitch_islandGlbSt(
    index_islands = island_indices,
    data = data,
    tree = tree_string,
    u_threshold = 0.2,
    m_threshold = 0.8,
    testing = FALSE
  )

  # Preserve Ferit's definition: total Fitch score divided by
  # the number of island regions.
  mean_fitch_per_island <- sum(raw_fitch) / length(island_indices)

  stats <- data.frame(
    islandMeanFreqP =
      MethEvolSIM::get_islandMeanFreqP(
        index_islands = island_indices,
        data = data,
        sample_n = sample_n,
        categorized_data = TRUE
      ),

    islandMeanFreqM =
      MethEvolSIM::get_islandMeanFreqM(
        index_islands = island_indices,
        data = data,
        sample_n = sample_n,
        categorized_data = TRUE
      ),

    islandSDFreqM =
      MethEvolSIM::get_islandSDFreqM(
        index_islands = island_indices,
        data = data,
        sample_n = sample_n,
        categorized_data = TRUE
      ),

    islandSDFreqP =
      MethEvolSIM::get_islandSDFreqP(
        index_islands = island_indices,
        data = data,
        sample_n = sample_n,
        categorized_data = TRUE
      ),

    meanCor_i =
      MethEvolSIM::compute_meanCor_i(
        index_islands = island_indices,
        minN_CpG = minN_sites,
        shore_length = shore_length,
        data = data,
        sample_n = sample_n,
        categorized_data = TRUE
      ),

    nonislandMeanFreqM =
      MethEvolSIM::get_nonislandMeanFreqM(
        index_nonislands = nonisland_indices,
        data = data,
        sample_n = sample_n,
        categorized_data = TRUE
      ),

    nonislandSDFreqM =
      MethEvolSIM::get_nonislandSDFreqM(
        index_nonislands = nonisland_indices,
        data = data,
        sample_n = sample_n,
        categorized_data = TRUE
      ),

    nonislandSDFreqP =
      MethEvolSIM::get_nonislandSDFreqP(
        index_nonislands = nonisland_indices,
        data = data,
        sample_n = sample_n,
        categorized_data = TRUE
      ),

    nonislandMeanFreqP =
      MethEvolSIM::get_nonislandMeanFreqP(
        index_nonislands = nonisland_indices,
        data = data,
        sample_n = sample_n,
        categorized_data = TRUE
      ),

    meanCor_ni =
      MethEvolSIM::compute_meanCor_ni(
        index_nonislands = nonisland_indices,
        minN_CpG = minN_sites,
        shore_length = shore_length,
        data = data,
        sample_n = sample_n,
        categorized_data = TRUE
      ),

    # The cichlid tree contains one cherry: benthic + littoral.
    MeanSiteFChange_i =
      cherry_stats$island_meanFChange[cherry_index],

    MeanSiteFChange_ni =
      cherry_stats$nonisland_meanFChange[cherry_index],

    mean_fitch_islands =
      mean_fitch_per_island,

    mean_TreeFreqsChange_i =
      MethEvolSIM::mean_TreeFreqsChange_i(
        tree = tree_string,
        data = data,
        categorized_data = TRUE,
        index_islands = island_indices,
        pValue_threshold = p_value_threshold
      ),

    stringsAsFactors = FALSE
  )

  finite_stats <- vapply(
    stats,
    function(x) length(x) == 1L && is.finite(x),
    logical(1)
  )

  if (any(!finite_stats)) {
    warning(
      "Non-finite summary statistic(s): ",
      paste(names(stats)[!finite_stats], collapse = ", "),
      ". Values are retained and are NOT replaced by zero."
    )
  }

  stats
}
