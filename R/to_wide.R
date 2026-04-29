# to_wide.R
# Reshape a long-format peaks data frame into a wide numeric matrix.

#' Reshape a long-format peaks table to a wide matrix.
#'
#' @param df        Long-format data frame with at least \code{id_col},
#'                  \code{name_col}, and \code{area_col} columns.
#'                  Typically the \code{$peaks} slot of a \code{ChromResult}.
#' @param area_col  Name of the area column to pivot (e.g. \code{"area_cons"}).
#' @param id_col    Column identifying the chromatogram. Default \code{"chrom"}.
#' @param name_col  Column with peak names. Default \code{"name"}.
#' @return Numeric matrix: rows = chromatograms, columns = named peaks.
#'   NA where a peak was not detected in a chromatogram.
#' @export
to_wide <- function(df, area_col,
                    id_col   = "chrom",
                    name_col = "name") {
  stopifnot(all(c(id_col, name_col, area_col) %in% names(df)))
  tmp <- df[!is.na(df[[name_col]]), c(id_col, name_col, area_col)]
  mat <- reshape(tmp,
                 idvar     = id_col,
                 timevar   = name_col,
                 direction = "wide")
  rownames(mat) <- mat[[id_col]]
  keep <- grepl(paste0("^", area_col, "\\."), names(mat))
  mat  <- mat[, keep, drop = FALSE]
  colnames(mat) <- sub(paste0("^", area_col, "\\."), "", colnames(mat))
  as.matrix(mat)
}
