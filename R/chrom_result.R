# chrom_result.R — Internal ChromResult constructor.
# S4 class definition : AllClasses.R
# S4 generics         : AllGenerics.R
# S4 methods          : AllMethods.R

.new_ChromResult <- function(meta,
                              peaks,
                              lod_loq,
                              alignment     = NULL,
                              baseline_rmse = NULL,
                              report        = NULL,
                              traces        = NULL) {
  new("ChromResult",
      meta          = meta,
      peaks         = peaks,
      lod_loq       = lod_loq,
      alignment     = alignment,
      baseline_rmse = baseline_rmse,
      report        = report,
      traces        = traces)
}
