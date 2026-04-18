# inhibit.R
# Replaces SWArbled: manually define time windows where baseline fitting
# should be suppressed (e.g. valve-switching artifacts, dead time).
#
# Args:
#   RT              : numeric vector of retention times in minutes
#   inhibit_regions : list of c(t_start, t_end) pairs in minutes
#                     e.g. list(c(0, 1.5), c(12.3, 13.0))
#
# Returns: list(block, segments)
#   block    : logical vector, TRUE = inhibited / excluded
#   segments : list of named c(begin=, end=) index pairs for good regions
#
# Drop-in replacement for swarbled(); same output format consumed by the
# SASS loop in read_chroma.R.

inhibit <- function(RT, inhibit_regions = list()) {
  n <- length(RT)
  block <- rep(FALSE, n)

  for (r in inhibit_regions) {
    block[RT >= r[1] & RT <= r[2]] <- TRUE
  }

  # Build good-segment list from contiguous FALSE runs
  runs  <- rle(block)
  ends   <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)

  segments <- list()
  for (i in seq_along(runs$values)) {
    if (!runs$values[i]) {
      segments[[length(segments) + 1]] <- c(begin = starts[i], end = ends[i])
    }
  }

  list(block = block, segments = segments)
}
