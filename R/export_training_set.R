# export_training_set.R
# Export a quality-gated mcescher peak table paired with raw signal windows
# to a Parquet file for machine-learning use.
#
# Each row in the output represents one (chromatogram × compound) pair and
# contains:
#   - the normalized, resampled signal window centred on the peak apex
#   - the mcescher consensus area + GUM uncertainty
#   - compound metadata
#
# Requires: nanoparquet, chromConverter (for Agilent .D), ncdf4 (for CDF)

#' Export a training dataset of paired (signal window, area label) tuples.
#'
#' Reads raw chromatogram files, extracts a fixed-width signal window around
#' each peak listed in a quality-gated mcescher peak table, resamples to
#' \code{n_samples} points, and writes the result to a Parquet file.  The
#' output is immediately usable for training neural-network integrators in
#' Python, Julia, or R without any knowledge of the underlying GC data format.
#'
#' @param peaks_df       Quality-gated peak table, typically the output of
#'                       \code{\link{loq_filter}(}\code{\link{replicate_quality_gate}(...))}
#'                       after removing rows with \code{quality_ok = FALSE}.
#'                       Must contain columns \code{chrom}, \code{name},
#'                       \code{RT}, \code{area_cons}, \code{u_cons},
#'                       \code{area_cv_pct}.
#' @param data_root      Root directory that contains the chromatogram files.
#'                       The function searches recursively for files whose
#'                       \code{basename} matches the \code{chrom} column.
#'                       Supports Agilent \code{.D} directories
#'                       (requires \code{chromConverter}) and AIA/CDF
#'                       \code{.cdf} files (requires \code{ncdf4}).
#' @param output_parquet Path to the output \code{.parquet} file.
#' @param window_half_min  Half-width of the extracted signal window in
#'                       minutes.  The window spans
#'                       \code{[RT - window_half_min, RT + window_half_min]}.
#'                       Default 0.6 min.
#' @param n_samples      Number of evenly spaced samples after resampling.
#'                       Default 512.  Must match the CNN input width.
#' @param normalize      If \code{TRUE} (default), apply min-max normalisation
#'                       to \code{[0, 1]} within each window before storing.
#'
#' @return Invisibly returns the output data frame.  Called for its side
#'   effect of writing \code{output_parquet}.
#'
#' @details
#' \strong{Output schema:}
#' \describe{
#'   \item{chrom}{Chromatogram filename (e.g. \code{"0428_1_A.D"}).}
#'   \item{name}{Compound name as labelled by mcescher.}
#'   \item{rt_apex}{Peak apex retention time (minutes).}
#'   \item{rt_start}{Window start (minutes).}
#'   \item{rt_end}{Window end (minutes).}
#'   \item{signal}{List of \code{n_samples} float32 values: the normalised,
#'     resampled signal window.}
#'   \item{area_cons}{Consensus peak area (mV·s).}
#'   \item{u_cons}{Standard uncertainty of \code{area_cons} (mV·s, k=1).}
#'   \item{area_cv_pct}{Relative standard uncertainty (\%).}
#'   \item{aic_winner}{Integration method selected by AICc.}
#' }
#'
#' \strong{Reading in Python:}
#' \preformatted{
#' import pandas as pd, numpy as np
#' df = pd.read_parquet("training_set.parquet")
#' X  = np.stack(df["signal"])          # (N, 512) float32
#' y  = df["area_cons"].values           # (N,)
#' }
#'
#' @seealso \code{\link{replicate_quality_gate}}, \code{\link{loq_filter}}
#'
#' @export
export_training_set <- function(
  peaks_df,
  data_root,
  output_parquet,
  window_half_min = 0.6,
  n_samples       = 512L,
  normalize       = TRUE
) {
  if (!requireNamespace("nanoparquet", quietly = TRUE))
    stop("Package 'nanoparquet' is required.  Install with: install.packages('nanoparquet')")

  data_root <- normalizePath(data_root, mustWork = TRUE)

  # Build chrom → full path index (search once, reuse)
  cat("Indexing chromatogram files under", data_root, "...\n")
  all_d   <- list.dirs(data_root, recursive = TRUE, full.names = TRUE)
  all_d   <- all_d[grepl("\\.D$", all_d, ignore.case = TRUE)]
  all_cdf <- list.files(data_root, pattern = "\\.cdf$",
                        recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

  path_index <- c(
    stats::setNames(all_d,   basename(all_d)),
    stats::setNames(all_cdf, basename(all_cdf))
  )
  cat(sprintf("  %d .D dirs + %d .cdf files indexed.\n",
              length(all_d), length(all_cdf)))

  chroms_needed <- unique(peaks_df$chrom)
  missing <- setdiff(chroms_needed, names(path_index))
  if (length(missing) > 0L)
    warning(sprintf("%d chromatogram(s) not found under data_root and will be skipped:\n  %s",
                    length(missing), paste(head(missing, 10), collapse = "\n  ")))

  # Process chromatograms
  rows <- vector("list", nrow(peaks_df))
  row_idx <- 0L
  n_chrom_ok <- 0L; n_chrom_fail <- 0L

  for (chrom in intersect(chroms_needed, names(path_index))) {
    path     <- path_index[[chrom]]
    is_d_dir <- grepl("\\.D$", path, ignore.case = TRUE)

    rt_sig <- tryCatch(
      .load_signal(path, is_d_dir),
      error = function(e) {
        message(sprintf("  skip %s: %s", chrom, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(rt_sig)) { n_chrom_fail <- n_chrom_fail + 1L; next }

    rt  <- rt_sig$rt
    sig <- rt_sig$signal
    n_chrom_ok <- n_chrom_ok + 1L

    peak_rows <- peaks_df[peaks_df$chrom == chrom, ]
    for (i in seq_len(nrow(peak_rows))) {
      p         <- peak_rows[i, ]
      rt_start  <- p$RT - window_half_min
      rt_end    <- p$RT + window_half_min
      mask      <- rt >= rt_start & rt <= rt_end
      if (sum(mask) < 8L) next

      seg <- sig[mask]
      x_old <- seq(0, 1, length.out = sum(mask))
      x_new <- seq(0, 1, length.out = n_samples)
      window <- stats::approx(x_old, seg, xout = x_new)$y

      if (normalize) {
        lo  <- min(window); hi <- max(window)
        rng <- hi - lo
        window <- if (rng > 1e-9) (window - lo) / rng else rep(0, n_samples)
      }

      row_idx <- row_idx + 1L
      rows[[row_idx]] <- list(
        chrom       = chrom,
        name        = p$name,
        rt_apex     = p$RT,
        rt_start    = rt_start,
        rt_end      = rt_end,
        signal      = list(as.numeric(window)),
        area_cons   = p$area_cons,
        u_cons      = p$u_cons,
        area_cv_pct = p$area_cv_pct,
        aic_winner  = if ("aic_winner" %in% names(p)) p$aic_winner else NA_character_
      )
    }
  }

  cat(sprintf("Chromatograms: %d loaded, %d failed.\n", n_chrom_ok, n_chrom_fail))

  if (row_idx == 0L) stop("No rows produced; check data_root and peaks_df.")
  rows <- rows[seq_len(row_idx)]

  out <- data.frame(
    chrom       = vapply(rows, `[[`, "", "chrom"),
    name        = vapply(rows, `[[`, "", "name"),
    rt_apex     = vapply(rows, `[[`, 0.0, "rt_apex"),
    rt_start    = vapply(rows, `[[`, 0.0, "rt_start"),
    rt_end      = vapply(rows, `[[`, 0.0, "rt_end"),
    area_cons   = vapply(rows, `[[`, 0.0, "area_cons"),
    u_cons      = vapply(rows, `[[`, 0.0, "u_cons"),
    area_cv_pct = vapply(rows, `[[`, 0.0, "area_cv_pct"),
    aic_winner  = vapply(rows, `[[`, "", "aic_winner"),
    stringsAsFactors = FALSE
  )
  out$signal <- lapply(rows, function(r) r$signal[[1L]])

  nanoparquet::write_parquet(out, output_parquet)
  cat(sprintf("Written: %s  (%d rows, %d compounds, signal[%d])\n",
              output_parquet, nrow(out),
              length(unique(out$name)), n_samples))
  invisible(out)
}


# Internal: load (rt, signal) from a .D directory or .cdf file
.load_signal <- function(path, is_d_dir) {
  if (is_d_dir) {
    if (!requireNamespace("chromConverter", quietly = TRUE))
      stop("'chromConverter' required for .D files: install.packages('chromConverter')")
    x  <- chromConverter::read_agilent_d(path, format_out = "data.frame",
                                         data_format = "long",
                                         read_metadata = FALSE)
    list(rt = x$rt, signal = x$intensity)
  } else {
    if (!requireNamespace("ncdf4", quietly = TRUE))
      stop("'ncdf4' required for .cdf files: install.packages('ncdf4')")
    ds  <- ncdf4::nc_open(path)
    on.exit(ncdf4::nc_close(ds), add = TRUE)
    delay    <- ncdf4::ncvar_get(ds, "actual_delay_time")
    interval <- ncdf4::ncvar_get(ds, "actual_sampling_interval")
    sig      <- ncdf4::ncvar_get(ds, "ordinate_values")
    n   <- length(sig)
    rt  <- (delay + seq(0, n - 1) * interval) / 60.0
    list(rt = rt, signal = sig)
  }
}
