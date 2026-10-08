theta_to_methevol_params <- function(theta, alpha_Ri = 10) {

  required <- c(
    "E_pI", "phi_pI",
    "E_mI", "phi_mI",
    "E_pNI", "phi_pNI",
    "E_mNI", "phi_mNI",
    "iota", "L", "muL"
  )

  if (is.null(names(theta)) || !all(required %in% names(theta))) {
    stop("theta must be a named vector containing all 11 required parameters.")
  }

  theta <- theta[required]

  E_names   <- c("E_pI", "E_mI", "E_pNI", "E_mNI")
  phi_names <- c("phi_pI", "phi_mI", "phi_pNI", "phi_mNI")

  if (any(theta[E_names] <= 0 | theta[E_names] >= 1)) {
    stop("All E parameters must be between 0 and 1.")
  }

  if (any(theta[phi_names] <= 0)) {
    stop("All phi parameters must be > 0.")
  }

  if (theta["iota"] <= 0 || theta["iota"] > 1) {
    stop("iota must be > 0 and <= 1.")
  }

  if (theta["L"] <= 0 || theta["muL"] <= 0) {
    stop("L and muL must be > 0.")
  }

  params <- data.frame(
    alpha_pI  = theta["E_pI"]  * theta["phi_pI"],
    beta_pI   = (1 - theta["E_pI"])  * theta["phi_pI"],

    alpha_mI  = theta["E_mI"]  * theta["phi_mI"],
    beta_mI   = (1 - theta["E_mI"])  * theta["phi_mI"],

    alpha_pNI = theta["E_pNI"] * theta["phi_pNI"],
    beta_pNI  = (1 - theta["E_pNI"]) * theta["phi_pNI"],

    alpha_mNI = theta["E_mNI"] * theta["phi_mNI"],
    beta_mNI  = (1 - theta["E_mNI"]) * theta["phi_mNI"],

    alpha_Ri = alpha_Ri,
    mu = theta["muL"] / theta["L"],
    iota = theta["iota"]
  )

  rownames(params) <- NULL

  return(params)
}


scale_cichlid_tree <- function(tree_1scaled,
                               L,
                               root_branch_length = 50) {

  if (!inherits(tree_1scaled, "phylo")) {
    stop("tree_1scaled must be an ape phylo object.")
  }

  if (!is.finite(L) || L <= 0) {
    stop("L must be finite and > 0.")
  }

  if (!is.finite(root_branch_length) || root_branch_length <= 0) {
    stop("root_branch_length must be finite and > 0.")
  }

  scaled_tree <- tree_1scaled
  scaled_tree$edge.length <- scaled_tree$edge.length * L

  newick <- ape::write.tree(scaled_tree)

  newick_with_root <- paste0(
    "(",
    substr(newick, 1, nchar(newick) - 1),
    ":",
    root_branch_length,
    ");"
  )

  return(newick_with_root)
}


simulate_cichlid_particle <- function(theta,
                                      tree_1scaled,
                                      genome_spatial_str,
                                      chromosomes = NULL,
                                      seed = 1L,
                                      n_cores = 1L,
                                      alpha_Ri = 10,
                                      root_branch_length = 50,
                                      CFTP = FALSE) {

  if (!all(c("chrom", "n", "globalState") %in% names(genome_spatial_str))) {
    stop("genome_spatial_str must contain chrom, n, and globalState columns.")
  }

  if (!all(genome_spatial_str$globalState %in% c("M", "U"))) {
    stop("globalState must contain only M or U.")
  }

  params <- theta_to_methevol_params(
    theta = theta,
    alpha_Ri = alpha_Ri
  )

  scaled_tree <- scale_cichlid_tree(
    tree_1scaled = tree_1scaled,
    L = unname(theta["L"]),
    root_branch_length = root_branch_length
  )

  all_chr <- unique(genome_spatial_str$chrom)

  if (is.null(chromosomes)) {
    chr_ids <- all_chr
  } else {
    missing_chr <- setdiff(chromosomes, all_chr)

    if (length(missing_chr) > 0) {
      stop(
        "Requested chromosome(s) not present in genome_spatial_str: ",
        paste(missing_chr, collapse = ", ")
      )
    }

    chr_ids <- chromosomes
  }

  seed <- as.integer(seed)
  n_cores <- as.integer(n_cores)

  if (is.na(seed)) {
    stop("seed must be coercible to an integer.")
  }

  if (is.na(n_cores) || n_cores < 1) {
    stop("n_cores must be >= 1.")
  }

  chr_seeds <- (seed + seq_along(chr_ids) - 1L) %% .Machine$integer.max
  chr_seeds[chr_seeds == 0L] <- 1L

  run_one_chromosome <- function(i) {

    chr_id <- chr_ids[i]
    set.seed(chr_seeds[i])

    spatial_str <- genome_spatial_str[
      genome_spatial_str$chrom == chr_id,
      c("n", "globalState"),
      drop = FALSE
    ]

    output <- MethEvolSIM::simulate_evolData(
      infoStr = spatial_str,
      tree = scaled_tree,
      params = params,
      CFTP = CFTP
    )

    tip_data <- lapply(
      seq_along(tree_1scaled$tip.label),
      function(tip) output$data[[1]][[tip]]$seq
    )

    names(tip_data) <- tree_1scaled$tip.label

    expected_regions <- nrow(spatial_str)
    observed_regions <- vapply(tip_data, length, integer(1))

    if (any(observed_regions != expected_regions)) {
      stop(
        "Unexpected number of simulated regions for chromosome ",
        chr_id,
        ". Expected ", expected_regions,
        " regions; observed ",
        paste(observed_regions, collapse = ", "),
        "."
      )
    }

    region_lengths_ok <- vapply(
      tip_data,
      function(x) {
        identical(
          as.integer(vapply(x, length, integer(1))),
          as.integer(spatial_str$n)
        )
      },
      logical(1)
    )

    if (any(!region_lengths_ok)) {
      stop(
        "Simulated region lengths do not match spatial_str$n for chromosome ",
        chr_id,
        "."
      )
    }

    return(tip_data)
  }

  if (n_cores == 1L || length(chr_ids) == 1L) {
    data_allChromosomes <- lapply(
      seq_along(chr_ids),
      run_one_chromosome
    )
  } else {
    data_allChromosomes <- parallel::mclapply(
      seq_along(chr_ids),
      run_one_chromosome,
      mc.cores = n_cores,
      mc.set.seed = FALSE
    )
  }

  names(data_allChromosomes) <- chr_ids

  return(data_allChromosomes)
}
