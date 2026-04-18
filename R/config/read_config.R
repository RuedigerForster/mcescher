# read_config.R
# Read and validate a method YAML config; match detected peaks to named peaks.

library(yaml)

# ---------------------------------------------------------------------------
# read_method_config(path)
#
# Reads method.yaml and returns a validated, typed list.  Stops with a
# descriptive error if any required field is missing or has the wrong type.
# ---------------------------------------------------------------------------
read_method_config <- function(path) {
  if (!file.exists(path))
    stop(sprintf("Config file not found: %s", path))

  cfg <- yaml::read_yaml(path)

  # ---- required top-level keys -------------------------------------------
  required <- c("method", "inhibit", "peaks",
                "processing", "baseline", "detection", "integration")
  missing  <- setdiff(required, names(cfg))
  if (length(missing))
    stop(sprintf("method.yaml: missing required section(s): %s",
                 paste(missing, collapse = ", ")))

  # ---- type coercion / defaults ------------------------------------------

  # inhibit: list of length-2 numeric vectors
  cfg$inhibit <- lapply(cfg$inhibit, function(x) as.numeric(unlist(x)))

  # peaks: list of named lists with numeric RT_ref / RT_window
  cfg$peaks <- lapply(cfg$peaks, function(p) {
    p$RT_ref    <- as.numeric(p$RT_ref)
    p$RT_window <- as.numeric(p$RT_window)
    p
  })

  # processing
  cfg$processing$x_factor <- as.integer(cfg$processing$x_factor)
  cfg$processing$denoise  <- isTRUE(cfg$processing$denoise)

  # denoising (optional section; supply defaults if absent)
  if (is.null(cfg$denoising))
    cfg$denoising <- list(d = 1L, fc = 0.011, K = 1L, lam = 0.2)
  cfg$denoising$d   <- as.integer(cfg$denoising$d)
  cfg$denoising$K   <- as.integer(cfg$denoising$K)
  cfg$denoising$fc  <- as.numeric(cfg$denoising$fc)
  cfg$denoising$lam <- as.numeric(cfg$denoising$lam)

  # area_scale: convert integration output (counts·RT_unit) to counts·seconds
  # Empower and most CDS report areas in counts·s; RT is stored in minutes here.
  cfg$area_scale <- switch(tolower(as.character(cfg$method$RT_unit)),
    "min" = 60,
    "s"   = 1,
    "sec" = 1,
    { warning(sprintf("Unknown RT_unit '%s'; area_scale set to 1.", cfg$method$RT_unit)); 1 }
  )

  # baseline
  cfg$baseline$lambda <- as.numeric(cfg$baseline$lambda)
  cfg$baseline$ratio  <- as.numeric(cfg$baseline$ratio)

  # detection
  cfg$detection$amp_thresh          <- as.numeric(cfg$detection$amp_thresh)
  cfg$detection$smooth_width_factor <- as.integer(cfg$detection$smooth_width_factor)

  # integration
  cfg$integration$methods         <- as.character(unlist(cfg$integration$methods))
  cfg$integration$ensemble_method <- as.character(cfg$integration$ensemble_method)
  cfg$integration$run_ensemble    <- isTRUE(cfg$integration$run_ensemble)

  cfg
}


# ---------------------------------------------------------------------------
# match_peak_names(peaks_df, config)
#
# Adds a `name` column to peaks_df by matching each detected peak to the
# nearest configured peak (by RT) within that peak's RT_window.
#
# Matching is performed independently per chromatogram so that RT drift
# between runs does not cause cross-assignment.
#
# Warns when a configured peak has no detected match in a chromatogram.
# Unmatched detected peaks receive name = NA_character_.
# ---------------------------------------------------------------------------
match_peak_names <- function(peaks_df, config) {
  peaks_df$name <- NA_character_

  if (nrow(peaks_df) == 0L || length(config$peaks) == 0L)
    return(peaks_df)

  cfg_name <- vapply(config$peaks, `[[`, "", "name")
  cfg_RT   <- vapply(config$peaks, `[[`, 0.0, "RT_ref")
  cfg_win  <- vapply(config$peaks, `[[`, 0.0, "RT_window")

  for (chrom in unique(peaks_df$chrom)) {
    rows <- which(peaks_df$chrom == chrom)
    used <- integer(0)

    for (k in seq_along(cfg_name)) {
      cands <- setdiff(rows, used)
      if (length(cands) == 0L) {
        warning(sprintf("[%s] Peak '%s' (RT_ref=%.3f): no unassigned peaks remain.",
                        chrom, cfg_name[k], cfg_RT[k]))
        next
      }
      dists   <- abs(peaks_df$RT[cands] - cfg_RT[k])
      best_i  <- which.min(dists)
      if (dists[best_i] <= cfg_win[k]) {
        peaks_df$name[cands[best_i]] <- cfg_name[k]
        used <- c(used, cands[best_i])
      } else {
        warning(sprintf(
          "[%s] Peak '%s' (RT_ref=%.3f) not found — closest detected RT=%.3f (delta=%.3f > window=%.3f)",
          chrom, cfg_name[k], cfg_RT[k],
          peaks_df$RT[cands[best_i]], dists[best_i], cfg_win[k]))
      }
    }
  }

  # Re-order: name column right after peak_no / RT for readability
  nms <- names(peaks_df)
  anchor <- intersect(c("chrom", "peak_no", "RT"), nms)
  rest   <- setdiff(nms, c(anchor, "name"))
  peaks_df[, c(anchor, "name", rest)]
}
