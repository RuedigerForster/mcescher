# chrom_result.R
# S3 class ChromResult — unified container for chromatography pipeline output.
#
# Fields:
#   $meta          named list  — run metadata and method info
#   $peaks         data.frame  — one row per peak, all columns merged:
#                                chrom, peak_no, RT, name, height, width,
#                                area_PD/TS/gauss/EGH, fit_rmse_*,
#                                left/right_valley_RT,
#                                area_mean/sd/cv_pct/min/max/range (ensemble)
#   $lod_loq       data.frame  — one row per chromatogram
#   $alignment     data.frame | NULL
#   $baseline_rmse data.frame | NULL
#   $report        character   — full text report (cat() output from pipeline_summary)
#
# Methods:
#   print()         compact overview
#   summary()       full text report
#   as.data.frame() returns $peaks
#   [i]             subset to chromatogram(s) i (index or name)


# ---------------------------------------------------------------------------
# Constructor (internal)
# ---------------------------------------------------------------------------
.new_ChromResult <- function(meta,
                              peaks,
                              lod_loq,
                              alignment     = NULL,
                              baseline_rmse = NULL,
                              report        = NULL) {
  structure(
    list(
      meta          = meta,
      peaks         = peaks,
      lod_loq       = lod_loq,
      alignment     = alignment,
      baseline_rmse = baseline_rmse,
      report        = report
    ),
    class = "ChromResult"
  )
}


# ---------------------------------------------------------------------------
# print.ChromResult — compact single-screen summary
# ---------------------------------------------------------------------------
print.ChromResult <- function(x, ...) {
  cat("ChromResult\n")

  if (!is.null(x$meta$method_name))
    cat(sprintf("  Method        : %s\n", x$meta$method_name))

  cat(sprintf("  Chromatograms : %d\n", x$meta$n_chrom))
  cat(sprintf("  RT range      : %.3f \u2013 %.3f min\n",
              x$meta$RT_range[1], x$meta$RT_range[2]))
  cat(sprintf("  Samples/chrom : %d  (dt = %.4f min)\n",
              x$meta$n_samp, x$meta$dt_min))

  if (!is.null(x$meta$integrate_methods))
    cat(sprintf("  Integration   : %s\n",
                paste(x$meta$integrate_methods, collapse = ", ")))

  cat("\n")

  if (!is.null(x$peaks) && nrow(x$peaks) > 0L) {
    has_names <- "name" %in% names(x$peaks)

    for (chrom in unique(x$peaks$chrom)) {
      pk <- x$peaks[x$peaks$chrom == chrom, , drop = FALSE]
      cat(sprintf("  %s  (%d peak%s)\n",
                  chrom, nrow(pk), if (nrow(pk) == 1L) "" else "s"))

      if (has_names) {
        lbl <- ifelse(is.na(pk$name),
                      sprintf("RT %.3f [unnamed]", pk$RT),
                      sprintf("%s (RT %.3f)", pk$name, pk$RT))
      } else {
        lbl <- sprintf("RT %.3f", pk$RT)
      }

      # Show area_mean if ensemble was run, otherwise area_TS
      if ("area_mean" %in% names(pk)) {
        area_vals <- sprintf("area=%.1f (CV %.2f%%)",
                             pk$area_mean, pk$area_cv_pct)
      } else if ("area_TS" %in% names(pk)) {
        area_vals <- sprintf("area_TS=%.4f", pk$area_TS)
      } else {
        area_vals <- rep("", nrow(pk))
      }

      for (i in seq_len(nrow(pk)))
        cat(sprintf("    [%d] %s  %s\n", pk$peak_no[i], lbl[i], area_vals[i]))
    }
  } else {
    cat("  No peaks detected.\n")
  }

  cat("\n")

  if (!is.null(x$lod_loq) && nrow(x$lod_loq) > 0L) {
    lod <- round(x$lod_loq$LOD[1], 4)
    loq <- round(x$lod_loq$LOQ[1], 4)
    cat(sprintf("  LOD = %s  |  LOQ = %s  (noise_sd = %s)\n",
                lod, loq, round(x$lod_loq$noise_sd[1], 5)))
  }

  cat("  Use summary() for the full pipeline report.\n")
  invisible(x)
}


# ---------------------------------------------------------------------------
# summary.ChromResult — print the full captured report
# ---------------------------------------------------------------------------
summary.ChromResult <- function(object, ...) {
  if (!is.null(object$report) && length(object$report) > 0L)
    cat(object$report, sep = "")
  else
    print(object)
  invisible(object)
}


# ---------------------------------------------------------------------------
# as.data.frame.ChromResult — return the unified peaks table
# ---------------------------------------------------------------------------
as.data.frame.ChromResult <- function(x, row.names = NULL,
                                       optional = FALSE, ...) {
  x$peaks
}


# ---------------------------------------------------------------------------
# [.ChromResult — subset by chromatogram index or name
# ---------------------------------------------------------------------------
`[.ChromResult` <- function(x, i) {
  nms <- x$meta$col_names
  if (is.character(i)) i <- match(i, nms)
  i <- i[!is.na(i) & i >= 1L & i <= length(nms)]
  if (length(i) == 0L)
    stop("No matching chromatograms.")

  sel <- nms[i]
  x$meta$col_names <- sel
  x$meta$n_chrom   <- length(i)

  if (!is.null(x$peaks) && nrow(x$peaks) > 0L)
    x$peaks <- x$peaks[x$peaks$chrom %in% sel, , drop = FALSE]

  if (!is.null(x$lod_loq))
    x$lod_loq <- x$lod_loq[x$lod_loq$chrom_name %in% sel, , drop = FALSE]

  if (!is.null(x$alignment))
    x$alignment <- x$alignment[x$alignment$chrom %in% sel, , drop = FALSE]

  if (!is.null(x$baseline_rmse))
    x$baseline_rmse <- x$baseline_rmse[
      x$baseline_rmse$chrom %in% sel, , drop = FALSE]

  x$report <- NULL   # report is invalidated after subsetting
  x
}
