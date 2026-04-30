# --------------------------------------------------------------------------- #
#  Accessor methods                                                            #
# --------------------------------------------------------------------------- #

setMethod("peaks",  "ChromResult", function(object) object@peaks)
setMethod("meta",   "ChromResult", function(object) object@meta)
setMethod("traces", "ChromResult", function(object) object@traces)


# --------------------------------------------------------------------------- #
#  show — compact single-screen summary (replaces print.ChromResult)          #
# --------------------------------------------------------------------------- #

setMethod("show", "ChromResult", function(object) {
  cat("ChromResult\n")

  if (!is.null(object@meta$method_name))
    cat(sprintf("  Method        : %s\n", object@meta$method_name))

  cat(sprintf("  Chromatograms : %d\n", object@meta$n_chrom))
  cat(sprintf("  RT range      : %.3f – %.3f min\n",
              object@meta$RT_range[1], object@meta$RT_range[2]))
  cat(sprintf("  Samples/chrom : %d  (dt = %.4f min)\n",
              object@meta$n_samp, object@meta$dt_min))

  if (!is.null(object@meta$integrate_methods))
    cat(sprintf("  Integration   : %s\n",
                paste(object@meta$integrate_methods, collapse = ", ")))

  cat("\n")

  pk <- object@peaks
  if (!is.null(pk) && nrow(pk) > 0L) {
    has_names <- "name" %in% names(pk)

    for (chrom in unique(pk$chrom)) {
      pkc <- pk[pk$chrom == chrom, , drop = FALSE]
      cat(sprintf("  %s  (%d peak%s)\n",
                  chrom, nrow(pkc), if (nrow(pkc) == 1L) "" else "s"))

      if (has_names) {
        lbl <- ifelse(is.na(pkc$name),
                      sprintf("RT %.3f [unnamed]", pkc$RT),
                      sprintf("%s (RT %.3f)", pkc$name, pkc$RT))
      } else {
        lbl <- sprintf("RT %.3f", pkc$RT)
      }

      if ("area_mean" %in% names(pkc)) {
        area_vals <- sprintf("area=%.1f (CV %.2f%%)", pkc$area_mean, pkc$area_cv_pct)
      } else if ("area_TS" %in% names(pkc)) {
        area_vals <- sprintf("area_TS=%.4f", pkc$area_TS)
      } else {
        area_vals <- rep("", nrow(pkc))
      }

      for (i in seq_len(nrow(pkc)))
        cat(sprintf("    [%d] %s  %s\n", pkc$peak_no[i], lbl[i], area_vals[i]))
    }
  } else {
    cat("  No peaks detected.\n")
  }

  cat("\n")

  ll <- object@lod_loq
  if (!is.null(ll) && nrow(ll) > 0L) {
    lod <- round(ll$LOD[1], 4)
    loq <- round(ll$LOQ[1], 4)
    cat(sprintf("  LOD = %s  |  LOQ = %s  (noise_sd = %s)\n",
                lod, loq, round(ll$noise_sd[1], 5)))
  }

  cat("  Use summary() for the full pipeline report.\n")
  invisible(object)
})


# --------------------------------------------------------------------------- #
#  summary — print the full captured pipeline report                          #
# --------------------------------------------------------------------------- #

setMethod("summary", "ChromResult", function(object, ...) {
  rpt <- object@report
  if (!is.null(rpt) && length(rpt) > 0L)
    cat(rpt, sep = "")
  else
    show(object)
  invisible(object)
})


# --------------------------------------------------------------------------- #
#  as.data.frame — return the peaks table                                     #
# --------------------------------------------------------------------------- #

setMethod("as.data.frame", "ChromResult",
  function(x, row.names = NULL, optional = FALSE, ...) {
    x@peaks
  }
)


# --------------------------------------------------------------------------- #
#  [ — subset by chromatogram index or name                                   #
# --------------------------------------------------------------------------- #

setMethod("[", "ChromResult", function(x, i, j, ..., drop = TRUE) {
  nms <- x@meta$col_names
  if (is.character(i)) i <- match(i, nms)
  i <- i[!is.na(i) & i >= 1L & i <= length(nms)]
  if (length(i) == 0L)
    stop("No matching chromatograms.")

  sel <- nms[i]
  x@meta$col_names <- sel
  x@meta$n_chrom   <- length(i)

  pk <- x@peaks
  if (!is.null(pk) && nrow(pk) > 0L)
    x@peaks <- pk[pk$chrom %in% sel, , drop = FALSE]

  ll <- x@lod_loq
  if (!is.null(ll))
    x@lod_loq <- ll[ll$chrom_name %in% sel, , drop = FALSE]

  al <- x@alignment
  if (!is.null(al) && is.data.frame(al))
    x@alignment <- al[al$chrom %in% sel, , drop = FALSE]

  br <- x@baseline_rmse
  if (!is.null(br) && is.data.frame(br))
    x@baseline_rmse <- br[br$chrom %in% sel, , drop = FALSE]

  x@report <- NULL
  tr <- x@traces
  if (!is.null(tr))
    x@traces <- tr[i]
  x
})


# --------------------------------------------------------------------------- #
#  plot — raw signal + consensus baseline + uncertainty band                  #
# --------------------------------------------------------------------------- #

setMethod("plot", "ChromResult", function(x, y, ...) {
  args <- list(...)
  chrom     <- if (!is.null(args$chrom))     args$chrom     else 1L
  col_signal <- if (!is.null(args$col_signal)) args$col_signal else .chrom_plot_palette$signal
  col_bl     <- if (!is.null(args$col_bl))    args$col_bl    else .chrom_plot_palette$baseline
  col_band   <- if (!is.null(args$col_band))  args$col_band  else .chrom_plot_palette$band
  main       <- args$main
  xlab       <- if (!is.null(args$xlab)) args$xlab else "Retention time (min)"
  ylab       <- if (!is.null(args$ylab)) args$ylab else "Signal"
  xlim       <- args$xlim
  ylim       <- args$ylim

  tr_list <- x@traces
  if (is.null(tr_list) || length(tr_list) < chrom || is.null(tr_list[[chrom]]))
    stop("No trace data in this ChromResult  re-run pipeline_summary with run_ensemble = TRUE.")

  tr   <- tr_list[[chrom]]
  name <- if (!is.null(x@meta$col_names)) x@meta$col_names[[chrom]] else paste0("chrom_", chrom)

  if (is.null(main))
    main <- sprintf("%s    %s",
                    if (!is.null(x@meta$method_name)) x@meta$method_name else "ChromResult",
                    name)

  pk <- x@peaks
  pk <- if (!is.null(pk) && nrow(pk) > 0L)
          pk[pk$chrom == name, , drop = FALSE]
        else
          data.frame()

  yr <- range(c(tr$signal,
                if (!is.null(tr$bl_hi)) tr$bl_hi else tr$bl_cons,
                if (!is.null(tr$bl_lo)) tr$bl_lo else tr$bl_cons),
              na.rm = TRUE)
  yr[1] <- min(yr[1], 0)

  plot(tr$RT, tr$signal, type = "n",
       xlim = if (is.null(xlim)) range(tr$RT) else xlim,
       ylim = if (is.null(ylim)) yr else ylim,
       xlab = xlab, ylab = ylab, main = main)

  if (!is.null(tr$bl_lo) && !is.null(tr$bl_hi)) {
    ok <- !is.na(tr$bl_lo) & !is.na(tr$bl_hi)
    polygon(c(tr$RT[ok], rev(tr$RT[ok])),
            c(tr$bl_hi[ok], rev(tr$bl_lo[ok])),
            col = col_band, border = NA)
  }

  lines(tr$RT, tr$signal,  col = col_signal, lwd = 1)
  lines(tr$RT, tr$bl_cons, col = col_bl,     lwd = 1.5, lty = 2)
  abline(h = 0, col = "grey70", lty = 3)

  if (nrow(pk) > 0L) {
    abline(v = pk$RT, col = .chrom_plot_palette$peak_ref, lty = 2, lwd = 0.8)
    lbl <- if ("name" %in% names(pk))
             ifelse(is.na(pk$name), sprintf("%.3f", pk$RT), pk$name)
           else
             sprintf("%.3f", pk$RT)
    y_top <- par("usr")[4]
    text(pk$RT, y_top, lbl, adj = c(1.1, 1), srt = 90,
         cex = 0.65, xpd = FALSE)
  }

  legend("topright", bty = "n",
         legend = c("Signal", "Consensus baseline", "±2σ band"),
         col    = c(col_signal, col_bl, adjustcolor(col_band, alpha.f = 1)),
         lty    = c(1, 2, NA),
         lwd    = c(1, 1.5, NA),
         pch    = c(NA, NA, 15),
         pt.cex = 1.8)

  invisible(x)
})
