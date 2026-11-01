# ============================================================================
# Leave-one-out Analysis - Standalone Script
# Processes both partitions ("co_cy_pass" calendar-year and "co_sy_pass"
# shifted-year) in a single sequential run.
# ============================================================================
# Designed to run detached on a cluster (e.g. via `Rscript run_loo.R` under
# SLURM/sbatch/SGE, nohup, or similar).
#
# Saves ONE .RData checkpoint file per (partition, country). At the end, all
# checkpoint files (both partitions) are consolidated into ONE final .RData
# file for easy loading into the downstream R Markdown report.
#
# INPUT: ONE .RData "bundle" file containing BOTH partitions' RAW data
# (NOT pre-standardized -- standardize_columns() is run INSIDE this script,
# per partition, using each partition's own dataset_type), organized as
# THREE named list objects, each keyed by partition name:
#   raw_df_list           -- raw_df_list$co_cy_pass, raw_df_list$co_sy_pass
#                            (each the RAW data frame with original columns:
#                            COUNTRY, COUNTRY_CODE, and either YEAR/WEEK
#                            (calendar) or shifted_year/shifted_week
#                            (shifted), plus inf_all_burden)
#   original_labels_list  -- original_labels_list$co_cy_pass, etc.
#                            (each a named vector, ISO3 -> cluster number --
#                            your chosen final clustering from the FULL,
#                            all-countries analysis, for THAT partition)
#   lat_lookup_list       -- lat_lookup_list$co_cy_pass, etc.
#                            (each a data frame with columns country_code,
#                            latitude -- can be identical content across
#                            partitions, since latitude doesn't depend on
#                            calendar vs. shifted year, but still stored
#                            per-partition for a uniform structure)
#
# NOTE: country_lookup is NOT needed as a separate input -- it is derived
# automatically inside this script by running standardize_columns() on each
# partition's raw_df, exactly as the main pipeline does.
#
# USAGE:
#   Rscript run_loo.R

# Packages
library(dplyr)
library(tidyr)
library(tibble)
library(cluster)
library(factoextra)
library(igraph)
library(kernlab)


# Configuration
CHECKPOINT_DIR   <- "loo_checkpoints"
FINAL_OUTPUT     <- "loo_analysis_final.RData"
LOG_FILE         <- "loo_analysis_log.txt"

INPUT_BUNDLE_PATH <- "loo_input_bundle.RData"

# Partition definitions:
# Each partition has its own dataset_type (controls which raw columns
# standardize_columns() pulls year/week from - "calendar" uses YEAR/WEEK,
# "shifted" uses shifted_year/shifted_week) AND its own algorithm/k/final-
# method settings, since the two datasets may have different chosen settings.
# Data for each partition comes from the single input bundle's
# raw_df_list/original_labels_list/lat_lookup_list, keyed by these same
# partition names.
PARTITIONS <- list(
  co_cy_pass = list(
    dataset_type = "calendar",
    chosen_algo = "kmeans",
    k_final = 2,
    final_algorithm = "spectral",
    n_boot = 1000,
    resample_pct = 0.9
  ),
  co_sy_pass = list(
    dataset_type = "shifted",
    chosen_algo = "kmeans",
    k_final = 3,
    final_algorithm = "spectral",
    n_boot = 1000,
    resample_pct = 0.9
  )
)

WEEK_DENOM <- 52 # matches the WEEK_DENOM used in the original analysis
SEED       <- 1234



### Functions

# Standardize raw columns based on dataset_type:
# Identical logic to the main pipeline's standardize_columns(): maps the raw
# COUNTRY/COUNTRY_CODE/YEAR-or-shifted_year/WEEK-or-shifted_week/inf_all_burden
# columns into the generic country/country_name/year_col/week_col/value shape
# that build_year_matrices() expects. This is THE step that makes calendar
# vs. shifted-year data get read from the correct source columns -- it must
# run BEFORE build_year_matrices() for each partition, using that partition's
# own dataset_type.
standardize_columns <- function(raw_df, dataset_type) {
  year_col <- switch(dataset_type, "calendar" = "YEAR", "shifted" = "shifted_year",
                      stop("dataset_type must be 'calendar' or 'shifted'"))
  week_col <- switch(dataset_type, "calendar" = "WEEK", "shifted" = "shifted_week",
                      stop("dataset_type must be 'calendar' or 'shifted'"))

  required_cols <- c("COUNTRY", "COUNTRY_CODE", year_col, week_col, "inf_all_burden")
  missing_cols <- setdiff(required_cols, names(raw_df))
  if (length(missing_cols) > 0) stop("Raw df is missing expected column(s): ", paste(missing_cols, collapse = ", "))

  df <- raw_df %>%
    transmute(
      country      = COUNTRY_CODE,
      country_name = COUNTRY,
      year_col     = .data[[year_col]],
      week_col     = .data[[week_col]],
      value        = inf_all_burden
    )

  country_lookup <- df %>% distinct(country, country_name) %>% deframe()

  list(df = df, country_lookup = country_lookup)
}

# Build one wide matrix per year from long-format df
build_year_matrices <- function(df, week_denom = WEEK_DENOM) {
  years <- sort(unique(df$year_col))
  year_mats <- list()

  for (y in years) {
    wide <- df %>%
      filter(year_col == y) %>%
      select(country, week_col, value) %>%
      pivot_wider(names_from = week_col, values_from = value) %>%
      arrange(country)

    mat <- as.matrix(wide %>% select(-country))
    rownames(mat) <- wide$country

    complete_rows <- complete.cases(mat)
    if (!all(complete_rows)) {
      warning(sprintf("Year %s: dropping %d countries with incomplete weeks", y, sum(!complete_rows)))
      mat <- mat[complete_rows, , drop = FALSE]
    }

    n_weeks <- ncol(mat)
    year_mats[[as.character(y)]] <- list(
      matrix = mat, n_weeks = n_weeks, weight = n_weeks / week_denom, countries = rownames(mat)
    )
  }
  year_mats
}

# Base clustering dispatch
run_single_clustering <- function(mat, k, algo) {
  d <- dist(mat, method = "euclidean")
  labels <- switch(algo,
    "kmeans"      = kmeans(mat, centers = k, nstart = 10, iter.max = 50)$cluster,
    "hclust_ward" = cutree(hclust(d, method = "ward.D2"), k = k),
    "hkmeans"     = factoextra::hkmeans(mat, k = k)$cluster,
    "pam"         = cluster::pam(d, k = k)$clustering,
    stop("Unknown algorithm: ", algo)
  )
  names(labels) <- rownames(mat)
  labels
}

# Fixed-percentage subsample without replacement
draw_resample_indices <- function(n_year, resample_pct) {
  n_draw <- max(2, round(n_year * resample_pct))
  sample(seq_len(n_year), size = n_draw, replace = FALSE)
}

# Internal validity metrics (needed inside build_consensus_matrices)
compute_wss <- function(mat, labels) {
  sum(sapply(unique(labels), function(cl) {
    members <- mat[labels == cl, , drop = FALSE]
    if (nrow(members) < 2) return(0)
    center <- colMeans(members)
    sum(rowSums(sweep(members, 2, center)^2))
  }))
}

compute_calinski_harabasz <- function(mat, labels) {
  n <- nrow(mat); k <- length(unique(labels))
  if (k < 2 || k >= n) return(NA)
  overall_center <- colMeans(mat)
  wss <- compute_wss(mat, labels)
  bss <- sum(sapply(unique(labels), function(cl) {
    members <- mat[labels == cl, , drop = FALSE]
    nrow(members) * sum((colMeans(members) - overall_center)^2)
  }))
  (bss / (k - 1)) / (wss / (n - k))
}

compute_dunn_index <- function(mat_or_dist, labels, is_dist = FALSE) {
  d <- if (is_dist) as.matrix(mat_or_dist) else as.matrix(dist(mat_or_dist))
  clusters <- unique(labels)
  if (length(clusters) < 2) return(NA)
  inter_dists <- c()
  for (i in seq_along(clusters)) for (j in seq_along(clusters)) if (i < j) {
    idx_i <- which(labels == clusters[i]); idx_j <- which(labels == clusters[j])
    inter_dists <- c(inter_dists, min(d[idx_i, idx_j, drop = FALSE]))
  }
  intra_diams <- sapply(clusters, function(cl) {
    idx <- which(labels == cl)
    if (length(idx) < 2) return(0)
    max(d[idx, idx, drop = FALSE])
  })
  if (max(intra_diams) == 0) return(NA)
  min(inter_dists) / max(intra_diams)
}

# Name-alignment safety check (used in extract_final_clusters_graph)
.verify_name_alignment <- function(source_names_before_overwrite, expected_names, context = "") {
  if (length(source_names_before_overwrite) != length(expected_names)) {
    stop(sprintf("Name alignment check failed (%s): length mismatch (%d vs %d)",
                  context, length(source_names_before_overwrite), length(expected_names)))
  }
  if (!is.null(source_names_before_overwrite) &&
      !identical(as.character(source_names_before_overwrite), as.character(expected_names))) {
    stop(sprintf("Name alignment check failed (%s): output names do not match expected order.", context))
  }
  invisible(TRUE)
}

# Build weighted co-association (consensus) matrices
build_consensus_matrices <- function(year_mats, k_range, algos, n_boot, resample_pct) {
  all_countries <- sort(unique(unlist(lapply(year_mats, `[[`, "countries"))))
  n_all <- length(all_countries)
  knames <- paste0("k", k_range)

  blank_mat <- function() matrix(0, n_all, n_all, dimnames = list(all_countries, all_countries))
  numer <- setNames(lapply(algos, function(a) setNames(lapply(k_range, function(k) blank_mat()), knames)), algos)
  denom <- setNames(lapply(algos, function(a) setNames(lapply(k_range, function(k) blank_mat()), knames)), algos)

  for (yname in names(year_mats)) {
    ymat <- year_mats[[yname]]$matrix
    weight <- year_mats[[yname]]$weight
    year_countries <- rownames(ymat)
    n_year <- length(year_countries)
    if (n_year < 2) next

    for (algo in algos) {
      for (b in seq_len(n_boot)) {
        boot_idx <- draw_resample_indices(n_year, resample_pct)
        boot_mat <- ymat[boot_idx, , drop = FALSE]
        unique_present <- rownames(boot_mat)
        if (length(unique_present) < 2) next

        for (k in k_range) {
          k_eff <- min(k, length(unique_present))
          if (k_eff < 2) next
          kname <- paste0("k", k)

          labels <- run_single_clustering(boot_mat, k_eff, algo)

          present <- unique_present
          for (i in seq_len(length(present) - 1)) {
            ci <- present[i]; li <- labels[ci]
            for (j in (i + 1):length(present)) {
              cj <- present[j]; lj <- labels[cj]
              denom[[algo]][[kname]][ci, cj] <- denom[[algo]][[kname]][ci, cj] + weight
              denom[[algo]][[kname]][cj, ci] <- denom[[algo]][[kname]][cj, ci] + weight
              if (li == lj) {
                numer[[algo]][[kname]][ci, cj] <- numer[[algo]][[kname]][ci, cj] + weight
                numer[[algo]][[kname]][cj, ci] <- numer[[algo]][[kname]][cj, ci] + weight
              }
            }
          }
        }
      }
    }
  }

  consensus <- setNames(lapply(algos, function(algo) {
    setNames(lapply(knames, function(kn) {
      cm <- numer[[algo]][[kn]] / denom[[algo]][[kn]]
      diag(cm) <- 1
      cm[is.nan(cm)] <- NA
      cm
    }), knames)
  }), algos)

  list(consensus = consensus, countries = all_countries)
}

# Graph-based final clustering (louvain/leiden/walktrap/spectral)
extract_final_clusters_graph <- function(cm, k_final, country_lookup = NULL,
                                          algorithm = c("louvain", "leiden", "walktrap", "spectral"),
                                          edge_threshold = 0, resolution_range = seq(0.1, 3, by = 0.05), seed = 1234) {
  algorithm <- match.arg(algorithm)
  set.seed(seed)

  cm_fill <- cm; cm_fill[is.na(cm_fill)] <- 0; diag(cm_fill) <- 0
  cm_fill[cm_fill < edge_threshold] <- 0
  g <- igraph::graph_from_adjacency_matrix(cm_fill, mode = "undirected", weighted = TRUE, diag = FALSE)

  isolated <- igraph::degree(g) == 0
  if (any(isolated)) {
    warning(sprintf("%d countries have no edges above threshold (%g) and are excluded: %s",
                     sum(isolated), edge_threshold, paste(igraph::V(g)$name[isolated], collapse = ", ")))
    g <- igraph::delete_vertices(g, which(isolated))
  }
  active <- igraph::V(g)$name
  n_available <- length(active)

  if (algorithm %in% c("louvain", "leiden")) {
    run_comm <- function(r) {
      if (algorithm == "louvain") igraph::cluster_louvain(g, weights = igraph::E(g)$weight, resolution = r)
      else igraph::cluster_leiden(g, weights = igraph::E(g)$weight, objective_function = "modularity", resolution = r)
    }
    search <- Filter(Negate(is.null), lapply(resolution_range, function(r) {
      comm <- tryCatch(run_comm(r), error = function(e) NULL)
      if (is.null(comm)) return(NULL)
      mem <- igraph::membership(comm)
      mod <- tryCatch(igraph::modularity(g, mem, weights = igraph::E(g)$weight), error = function(e) NA)
      list(resolution = r, n_comm = length(unique(mem)), modularity = mod, comm = comm)
    }))
    matching <- Filter(function(x) x$n_comm == k_final, search)
    if (length(matching) == 0) {
      stop(sprintf("No resolution produced k = %d communities. Achieved: %s",
                    k_final, paste(sort(unique(sapply(search, `[[`, "n_comm"))), collapse = ", ")))
    }
    best <- matching[[which.max(sapply(matching, `[[`, "modularity"))]]
    mem <- igraph::membership(best$comm)
    .verify_name_alignment(names(mem), active, context = sprintf("extract_final_clusters_graph/%s", algorithm))
    labels <- setNames(as.integer(mem), active)

    extra <- list(resolution = best$resolution, modularity = best$modularity, community_obj = best$comm)
    method_str <- sprintf("%s on consensus graph, resolution=%.3f, k=%d", algorithm, best$resolution, k_final)

  } else if (algorithm == "walktrap") {
    if (k_final >= n_available) stop(sprintf("k_final (%d) must be less than available countries (%d)", k_final, n_available))
    comm <- igraph::cluster_walktrap(g, weights = igraph::E(g)$weight)
    mem <- igraph::cut_at(comm, no = k_final)
    if (length(mem) != length(active)) {
      stop(sprintf("Name alignment check failed (walktrap): cut_at() returned %d labels for %d active vertices",
                    length(mem), length(active)))
    }
    labels <- setNames(as.integer(mem), active)
    mod <- tryCatch(igraph::modularity(g, mem, weights = igraph::E(g)$weight), error = function(e) NA)

    extra <- list(modularity = mod, community_obj = comm)
    method_str <- sprintf("walktrap on consensus graph, k=%d", k_final)

  } else { # spectral
    if (k_final >= n_available) stop(sprintf("k_final (%d) must be less than available countries (%d)", k_final, n_available))
    sim_mat <- cm_fill[active, active]
    spec_result <- kernlab::specc(kernlab::as.kernelMatrix(sim_mat), centers = k_final)
    spec_labels <- as.integer(spec_result@.Data)
    if (length(spec_labels) != length(active)) {
      stop(sprintf("Name alignment check failed (spectral): specc() returned %d labels for %d active countries",
                    length(spec_labels), length(active)))
    }
    labels <- setNames(spec_labels, active)
    mod <- tryCatch(igraph::modularity(g, labels, weights = igraph::E(g)$weight), error = function(e) NA)

    extra <- list(modularity = mod, specc_obj = spec_result)
    method_str <- sprintf("spectral (kernlab::specc) on consensus matrix, k=%d", k_final)
  }

  if (any(isolated)) cat(sprintf("(%d isolated countries excluded)\n", sum(isolated)))

  result <- c(list(labels = labels, graph_obj = g, method = method_str), extra)
  if (!is.null(country_lookup)) {
    result$labels_df <- data.frame(country_code = names(labels), country_name = country_lookup[names(labels)],
                                    cluster = unname(labels)) %>% arrange(cluster, country_name)
  }
  result
}

# Single LOO iteration: rebuild ensemble + final-cluster without one country
run_loo_single_country <- function(df, excluded_country, algo, k, n_boot, resample_pct,
                                    final_algorithm = c("louvain", "leiden", "walktrap", "spectral"),
                                    country_lookup = NULL, seed = 1234) {
  final_algorithm <- match.arg(final_algorithm)

  df_excl <- df %>% filter(country != excluded_country)
  year_mats_excl <- build_year_matrices(df_excl)

  result <- build_consensus_matrices(year_mats_excl, k_range = k, algos = algo, n_boot = n_boot, resample_pct = resample_pct)
  cm_excl <- result$consensus[[algo]][[paste0("k", k)]]

  final <- tryCatch(
    extract_final_clusters_graph(cm_excl, k, country_lookup = country_lookup, algorithm = final_algorithm, seed = seed),
    error = function(e) {
      warning(sprintf("LOO for %s: final clustering failed (%s)", excluded_country, conditionMessage(e)))
      NULL
    }
  )

  list(excluded_country = excluded_country, consensus_matrix = cm_excl, final = final, year_mats = year_mats_excl)
}

# Latitude-proxy cluster assignment for the excluded country
assign_by_latitude_proxy <- function(excluded_country, final_labels, lat_lookup) {
  lat_lookup_vec <- setNames(lat_lookup$latitude, lat_lookup$country_code)

  excluded_lat <- lat_lookup_vec[excluded_country]
  if (is.na(excluded_lat)) {
    warning(sprintf("No latitude found for excluded country %s", excluded_country))
    return(list(proxy_cluster = NA, excluded_lat = NA, cluster_mean_lats = NA))
  }

  member_lats <- lat_lookup_vec[names(final_labels)]
  cluster_mean_lats <- tapply(member_lats, final_labels, mean, na.rm = TRUE)

  if (all(is.na(cluster_mean_lats))) {
    warning(sprintf("No latitude data available for any cluster member when excluding %s", excluded_country))
    return(list(proxy_cluster = NA, excluded_lat = excluded_lat, cluster_mean_lats = cluster_mean_lats))
  }

  proxy_cluster <- as.integer(names(cluster_mean_lats)[which.min(abs(cluster_mean_lats - excluded_lat))])
  list(proxy_cluster = proxy_cluster, excluded_lat = excluded_lat, cluster_mean_lats = cluster_mean_lats)
}



# Logging Helper:
log_msg <- function(...) {
  msg <- sprintf(...)
  timestamped <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(timestamped, "\n")
  cat(timestamped, "\n", file = LOG_FILE, append = TRUE)
}

# Loading and fully validating both partitions before starting any processing
# means a malformed bundle for EITHER partition is caught immediately, rather
# than discovering a problem with the second partition's data only after the
# first partition's (potentially many-hour) run has already completed.

if (!dir.exists(CHECKPOINT_DIR)) dir.create(CHECKPOINT_DIR, recursive = TRUE)

log_msg("=== LOO analysis starting (partitions: %s) ===", paste(names(PARTITIONS), collapse = ", "))
log_msg("Loading input bundle from %s", INPUT_BUNDLE_PATH)

if (!file.exists(INPUT_BUNDLE_PATH)) {
  stop(sprintf("Input bundle not found: %s", INPUT_BUNDLE_PATH))
}

bundle_env <- new.env()
loaded_names <- load(INPUT_BUNDLE_PATH, envir = bundle_env)

required_list_objects <- c("raw_df_list", "original_labels_list", "lat_lookup_list")
missing_list_objects <- setdiff(required_list_objects, loaded_names)
if (length(missing_list_objects) > 0) {
  stop(sprintf(
    "Input bundle is missing required object(s): %s\nObjects found in bundle: %s",
    paste(missing_list_objects, collapse = ", "), paste(loaded_names, collapse = ", ")
  ))
}

raw_df_list           <- bundle_env$raw_df_list
original_labels_list  <- bundle_env$original_labels_list
lat_lookup_list       <- bundle_env$lat_lookup_list

# Validate + standardize each partition's data:
# This is where calendar vs. shifted-year actually diverges: each partition's
# RAW df is run through standardize_columns() using THAT partition's own
# dataset_type, so co_cy_pass pulls from YEAR/WEEK and co_sy_pass pulls from
# shifted_year/shifted_week -- exactly mirroring the main pipeline's logic.
validate_partition_data <- function(partition_name) {
  part_cfg <- PARTITIONS[[partition_name]]

  missing_from <- c()
  if (!(partition_name %in% names(raw_df_list))) missing_from <- c(missing_from, "raw_df_list")
  if (!(partition_name %in% names(original_labels_list))) missing_from <- c(missing_from, "original_labels_list")
  if (!(partition_name %in% names(lat_lookup_list))) missing_from <- c(missing_from, "lat_lookup_list")

  if (length(missing_from) > 0) {
    stop(sprintf("Partition '%s' is missing from: %s", partition_name, paste(missing_from, collapse = ", ")))
  }

  raw_df <- raw_df_list[[partition_name]]
  original_labels <- original_labels_list[[partition_name]]
  lat_lookup <- lat_lookup_list[[partition_name]]

  log_msg("[%s] Standardizing raw columns using dataset_type = '%s'...", partition_name, part_cfg$dataset_type)
  std <- standardize_columns(raw_df, part_cfg$dataset_type)
  df <- std$df
  country_lookup <- std$country_lookup

  if (!all(c("country_code", "latitude") %in% names(lat_lookup))) {
    stop(sprintf("[%s] lat_lookup_list entry must have columns: country_code, latitude", partition_name))
  }
  if (is.null(names(original_labels)) || length(original_labels) == 0) {
    stop(sprintf("[%s] original_labels_list entry must be a named vector (names = country codes)", partition_name))
  }

  log_msg("[%s] Validated OK: df=%d rows/%d countries/%d years, original_labels=%d countries/%d clusters, lat_lookup=%d rows",
          partition_name, nrow(df), length(unique(df$country)), length(unique(df$year_col)),
          length(original_labels), length(unique(original_labels)), nrow(lat_lookup))

  list(df = df, country_lookup = country_lookup, original_labels = original_labels, lat_lookup = lat_lookup)
}

bundles <- setNames(lapply(names(PARTITIONS), validate_partition_data), names(PARTITIONS))

log_msg("Input bundle loaded, standardized, and validated successfully for partitions: %s", paste(names(PARTITIONS), collapse = ", "))

# Checkpoint files are tagged with the partition name so both partitions can
# safely share CHECKPOINT_DIR: loo_<partition>_<country>.RData

start_time <- Sys.time()

for (partition_name in names(PARTITIONS)) {
  part_cfg <- PARTITIONS[[partition_name]]
  bundle <- bundles[[partition_name]]

  all_countries <- names(bundle$original_labels)
  n_total <- length(all_countries)

  log_msg("=== Starting partition '%s': %d countries, algo='%s', k=%d, final='%s', n_boot=%d, resample_pct=%.2f ===",
          partition_name, n_total, part_cfg$chosen_algo, part_cfg$k_final,
          part_cfg$final_algorithm, part_cfg$n_boot, part_cfg$resample_pct)

  for (i in seq_along(all_countries)) {
    excluded <- all_countries[i]
    checkpoint_file <- file.path(CHECKPOINT_DIR, sprintf("loo_%s_%s.RData", partition_name, excluded))

    if (file.exists(checkpoint_file)) {
      log_msg("[%s][%d/%d] %s already has a checkpoint -- skipping", partition_name, i, n_total, excluded)
      next
    }

    log_msg("[%s][%d/%d] Starting LOO run for %s...", partition_name, i, n_total, excluded)
    iter_start <- Sys.time()

    loo_run <- tryCatch(
      run_loo_single_country(
        df = bundle$df,
        excluded_country = excluded,
        algo = part_cfg$chosen_algo,
        k = part_cfg$k_final,
        n_boot = part_cfg$n_boot,
        resample_pct = part_cfg$resample_pct,
        final_algorithm = part_cfg$final_algorithm,
        country_lookup = bundle$country_lookup,
        seed = SEED
      ),
      error = function(e) {
        log_msg("[%s][%d/%d] ERROR for %s: %s", partition_name, i, n_total, excluded, conditionMessage(e))
        NULL
      }
    )

    if (is.null(loo_run) || is.null(loo_run$final)) {
      loo_result_entry <- list(
        partition = partition_name,
        excluded_country = excluded,
        succeeded = FALSE,
        error_note = "LOO run or final clustering failed -- see log for details",
        consensus_matrix = NULL,
        final = NULL,
        proxy = NULL
      )
      save(loo_result_entry, file = checkpoint_file)
      log_msg("[%s][%d/%d] FAILED for %s -- failure record saved to %s", partition_name, i, n_total, excluded, checkpoint_file)
      next
    }

    proxy_result <- assign_by_latitude_proxy(excluded, loo_run$final$labels, bundle$lat_lookup)

    loo_result_entry <- list(
      partition = partition_name,
      excluded_country = excluded,
      succeeded = TRUE,
      error_note = NA,
      consensus_matrix = loo_run$consensus_matrix,
      final = loo_run$final,
      proxy = proxy_result,
      original_cluster = unname(bundle$original_labels[excluded])
    )

    save(loo_result_entry, file = checkpoint_file)

    iter_elapsed <- as.numeric(difftime(Sys.time(), iter_start, units = "mins"))
    log_msg("[%s][%d/%d] Completed %s in %.1f minutes -- checkpoint saved to %s",
            partition_name, i, n_total, excluded, iter_elapsed, checkpoint_file)
  }

  log_msg("=== Partition '%s' complete ===", partition_name)
}

total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))
log_msg("=== All partitions complete. Total elapsed: %.2f hours ===", total_elapsed)


# Consolidate to final file
log_msg("Consolidating all checkpoint files into %s...", FINAL_OUTPUT)

checkpoint_files <- list.files(CHECKPOINT_DIR, pattern = "^loo_.*\\.RData$", full.names = TRUE)
log_msg("Found %d checkpoint files across all partitions", length(checkpoint_files))

loo_results <- setNames(vector("list", length(PARTITIONS)), names(PARTITIONS))
for (p in names(PARTITIONS)) loo_results[[p]] <- list()

summary_rows <- list()

for (cf in checkpoint_files) {
  e <- new.env()
  load(cf, envir = e)
  entry <- e$loo_result_entry
  partition_name <- entry$partition
  excluded <- entry$excluded_country

  loo_results[[partition_name]][[excluded]] <- entry

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    partition = partition_name,
    excluded_country = excluded,
    excluded_country_name = if (!is.null(bundles[[partition_name]]$country_lookup)) bundles[[partition_name]]$country_lookup[excluded] else NA,
    original_cluster = if (entry$succeeded) entry$original_cluster else NA,
    loo_rerun_succeeded = entry$succeeded,
    proxy_cluster = if (entry$succeeded) entry$proxy$proxy_cluster else NA,
    excluded_latitude = if (entry$succeeded) entry$proxy$excluded_lat else NA,
    note = entry$error_note
  )
}

loo_summary <- do.call(rbind, summary_rows)
rownames(loo_summary) <- NULL

log_msg("Consolidated summary (both partitions):")
print(loo_summary)

# Per-partition run metadata, bundled together for the downstream Rmd
loo_run_metadata <- setNames(lapply(names(PARTITIONS), function(p) {
  list(
    partition           = p,
    dataset_type        = PARTITIONS[[p]]$dataset_type,
    run_date            = Sys.time(),
    chosen_algo         = PARTITIONS[[p]]$chosen_algo,
    k_final             = PARTITIONS[[p]]$k_final,
    final_algorithm     = PARTITIONS[[p]]$final_algorithm,
    n_boot              = PARTITIONS[[p]]$n_boot,
    resample_pct        = PARTITIONS[[p]]$resample_pct,
    seed                = SEED,
    n_countries         = length(bundles[[p]]$original_labels)
  )
}), names(PARTITIONS))
loo_run_metadata$total_elapsed_hours <- total_elapsed

# country_lookup and lat_lookup carried through per-partition (in case they
# differ slightly between the two, e.g. a country present in one but not
# the other) -- stored as named lists keyed by partition.
country_lookups <- setNames(lapply(names(PARTITIONS), function(p) bundles[[p]]$country_lookup), names(PARTITIONS))
original_labels_list <- setNames(lapply(names(PARTITIONS), function(p) bundles[[p]]$original_labels), names(PARTITIONS))
lat_lookups <- setNames(lapply(names(PARTITIONS), function(p) bundles[[p]]$lat_lookup), names(PARTITIONS))

save(loo_results, loo_summary, loo_run_metadata, country_lookups, original_labels_list, lat_lookups,
     file = FINAL_OUTPUT)

log_msg("=== Final consolidated output (both partitions) saved to %s ===", FINAL_OUTPUT)
log_msg("=== LOO analysis complete ===")