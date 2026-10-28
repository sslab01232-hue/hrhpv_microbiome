# Shared settings for all R scripts (run from the repository root)

data_dir    <- "data"
results_dir <- "results"

seed     <- 2026
depth    <- 40000   # rarefaction depth
n_rarefy <- 100
n_perm   <- 999
n_cores  <- 2       # set to 1 on Windows

data_path <- function(...) file.path(data_dir, ...)

out_dir <- function(...) {
  p <- file.path(results_dir, ...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

# sample_metadata.csv: SampleID, Group, Age, HPV52, Dx, n_major_hrHPV, CST, CST_subtype
# (rows in the sample order used in the analysis)
read_metadata <- function() {
  md <- read.csv(data_path("sample_metadata.csv"))
  md$Group <- factor(md$Group, c("single", "multiple"))
  md$HPV52 <- factor(md$HPV52, c("N", "Y"))
  rownames(md) <- md$SampleID
  md
}

read_biom_tsv <- function(f) {
  skip <- startsWith(readLines(f, n = 1), "# Constructed")
  read.delim(f, skip = skip, row.names = 1, check.names = FALSE)
}
