library(chromConverter)
# library(lubridate)

import_chromas <- function(path, pattern) {
  files <- list.files(path = path, pattern = pattern, ignore.case = TRUE)
  run_datetime <- sample_name <- source_file <- sample_type <- character()
  for (file in files) {
    if (file.exists(paste0(path, file))) {
      A <- chromConverter::read_cdf(path = paste0(path, file), format_out = "data.frame", metadata_format = "chromconverter", data_format = "long")
      RT = A[, 1]/60
      if (file == files[1]) {
        df <- data.frame(A[, 2])
      } else {
        df <- cbind(df, A[1:nrow(df), 2])
      }
      run_datetime <- c(run_datetime, as.character(attr(A,"run_datetime")))
      sample_name <- c(sample_name, attr(A, "sample_name"))
      source_file <-  c(source_file, gsub(path, '', attr(A, "source_file")))
      sample_type <- c(sample_type, as.character(attr(A, "sample_type")))
    }
  }
  df <- as.matrix(df)
  while (length(RT) != nrow(df)) {
    RT <- RT[-length(RT)]
  }
  
  attr(df, "dimnames") <- list(as.character(RT), run_datetime)
  attr(df, "retention_time") <- RT
  attr(df, "sample_name") <- sample_name
  attr(df, "run_datetime") <- run_datetime
  attr(df, "source_file") <- source_file
  attr(df, "sample_type") <- sample_type
  
  return(df)
}