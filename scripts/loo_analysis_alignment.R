# Post-HOC Fix: Realign LOO proxy_cluster
# This script:
#   1. Loads your consolidated results
#   2. Re-aligns each LOO run's raw labels to that partition's original
#      clustering (align_loo_labels_to_original())
#   3. Recomputes proxy_cluster using the NOW-ALIGNED labels
#   4. Overwrites loo_summary and adds aligned_labels to each loo_results
#      entry, then saves a corrected .RData file

INPUT_FILE  <- "loo_analysis_final.RData"
OUTPUT_FILE <- "loo_analysis_final_corrected.RData"

library(dplyr)

# Align LOO run labels to the ORIGINAL clustering:
# (identical logic to the main pipeline / cluster script -- see those for
# full documentation)
align_loo_labels_to_original <- function(loo_labels, original_labels) {
  common <- intersect(names(loo_labels), names(original_labels))
  if (length(common) == 0) return(NULL)

  overlap <- table(loo_labels[common], original_labels[common])

  mapping <- c()
  used_original <- c()
  for (r in rownames(overlap)[order(-apply(overlap, 1, max))]) {
    candidates <- overlap[r, ]
    candidates[names(candidates) %in% used_original] <- -1
    best_original <- names(candidates)[which.max(candidates)]
    mapping[r] <- best_original
    used_original <- c(used_original, best_original)
  }

  aligned <- as.integer(mapping[as.character(loo_labels)])
  names(aligned) <- names(loo_labels)
  aligned
}

# Latitude-proxy assignment (must be called on aligned labels)
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


# Load existing results
cat(sprintf("Loading existing results from %s...\n", INPUT_FILE))

e <- new.env()
loaded_names <- load(INPUT_FILE, envir = e)
cat(sprintf("Loaded objects: %s\n", paste(loaded_names, collapse = ", ")))

required <- c("loo_results", "loo_summary", "loo_run_metadata", "original_labels_list", "lat_lookups")
missing <- setdiff(required, loaded_names)
if (length(missing) > 0) {
  stop(sprintf(
    "Input file is missing required object(s): %s\nFound: %s\n\nNOTE: if your file used different object names (e.g. a single-partition\nstructure from an earlier version of the script), adjust the object\nnames referenced below accordingly.",
    paste(missing, collapse = ", "), paste(loaded_names, collapse = ", ")
  ))
}

loo_results           <- e$loo_results            # nested: loo_results[[partition]][[country]]
loo_summary           <- e$loo_summary
loo_run_metadata       <- e$loo_run_metadata
original_labels_list   <- e$original_labels_list
lat_lookups             <- e$lat_lookups
country_lookups         <- if ("country_lookups" %in% loaded_names) e$country_lookups else NULL


# Re-align
cat("\nRe-aligning labels and recomputing latitude proxy for every entry...\n")

n_fixed <- 0
n_skipped_already_failed <- 0
n_newly_failed_alignment <- 0

for (partition_name in names(loo_results)) {
  original_labels <- original_labels_list[[partition_name]]
  lat_lookup <- lat_lookups[[partition_name]]

  if (is.null(original_labels)) {
    warning(sprintf("No original_labels found for partition '%s' -- skipping", partition_name))
    next
  }
  if (is.null(lat_lookup)) {
    warning(sprintf("No lat_lookup found for partition '%s' -- skipping", partition_name))
    next
  }

  for (excluded in names(loo_results[[partition_name]])) {
    entry <- loo_results[[partition_name]][[excluded]]

    if (!isTRUE(entry$succeeded) || is.null(entry$final) || is.null(entry$final$labels)) {
      n_skipped_already_failed <- n_skipped_already_failed + 1
      next  # this run already failed originally -- nothing to fix
    }

    # Re-align the RAW labels (untouched, still stored from the original run)
    aligned_labels <- align_loo_labels_to_original(entry$final$labels, original_labels)

    if (is.null(aligned_labels)) {
      warning(sprintf("[%s] %s: alignment failed (no common countries) -- marking as failed",
                       partition_name, excluded))
      entry$succeeded <- FALSE
      entry$error_note <- "Post-hoc alignment failed (no common countries)"
      entry$aligned_labels <- NULL
      entry$proxy <- NULL
      loo_results[[partition_name]][[excluded]] <- entry
      n_newly_failed_alignment <- n_newly_failed_alignment + 1
      next
    }

    # Recompute the proxy using the now-aligned labels
    proxy_result <- assign_by_latitude_proxy(excluded, aligned_labels, lat_lookup)

    entry$aligned_labels <- aligned_labels
    entry$proxy <- proxy_result
    loo_results[[partition_name]][[excluded]] <- entry
    n_fixed <- n_fixed + 1
  }
}

cat(sprintf("\nFixed %d entries.\n", n_fixed))
cat(sprintf("Skipped %d entries that had already failed originally (nothing to fix).\n", n_skipped_already_failed))
if (n_newly_failed_alignment > 0) {
  cat(sprintf("WARNING: %d entries newly marked as failed because alignment could not be computed.\n", n_newly_failed_alignment))
}


# Rebuild loo_summary corrected loo_results
cat("\nRebuilding loo_summary from corrected results...\n")

summary_rows <- list()
for (partition_name in names(loo_results)) {
  for (excluded in names(loo_results[[partition_name]])) {
    entry <- loo_results[[partition_name]][[excluded]]

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      partition = partition_name,
      excluded_country = excluded,
      excluded_country_name = if (!is.null(country_lookups) && !is.null(country_lookups[[partition_name]])) {
        country_lookups[[partition_name]][excluded]
      } else NA,
      original_cluster = if (isTRUE(entry$succeeded)) unname(original_labels_list[[partition_name]][excluded]) else NA,
      loo_rerun_succeeded = isTRUE(entry$succeeded),
      proxy_cluster = if (isTRUE(entry$succeeded)) entry$proxy$proxy_cluster else NA,
      excluded_latitude = if (isTRUE(entry$succeeded)) entry$proxy$excluded_lat else NA,
      note = entry$error_note
    )
  }
}

loo_summary_corrected <- do.call(rbind, summary_rows)
rownames(loo_summary_corrected) <- NULL

cat("\nCorrected summary:\n")
print(loo_summary_corrected)

# Quick before/after comparison
cat("\n=== BEFORE vs AFTER comparison (proxy_cluster column only) ===\n")
comparison <- loo_summary %>%
  select(partition, excluded_country, original_cluster, proxy_cluster_OLD = proxy_cluster) %>%
  left_join(
    loo_summary_corrected %>% select(partition, excluded_country, proxy_cluster_NEW = proxy_cluster),
    by = c("partition", "excluded_country")
  ) %>%
  mutate(changed = !is.na(proxy_cluster_OLD) & !is.na(proxy_cluster_NEW) & proxy_cluster_OLD != proxy_cluster_NEW)

print(comparison)
cat(sprintf("\n%d / %d entries had their proxy_cluster value CHANGE after realignment.\n",
            sum(comparison$changed, na.rm = TRUE), nrow(comparison)))


## Save Corrected Output
loo_summary <- loo_summary_corrected

save(loo_results, loo_summary, loo_run_metadata, country_lookups, original_labels_list, lat_lookups,
     file = OUTPUT_FILE)

cat(sprintf("\n=== Corrected results saved to %s ===\n", OUTPUT_FILE))
cat("Your original file is untouched. Load the CORRECTED file in your Rmd going forward:\n")
cat(sprintf('  loo_data <- load_loo_results("%s")\n', OUTPUT_FILE))