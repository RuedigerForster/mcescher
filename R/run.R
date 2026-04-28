# run.R
# Full pipeline for a single chromatogram.
# Run from the R/ directory: source("./run.R")

source("./import/import_chromas.R")

raw <- import_chromas("../data/UNIS/", pattern = "^UNIS1160190\\.cdf$")
chroma_list <- list(raw)

source("./test_single.R")
