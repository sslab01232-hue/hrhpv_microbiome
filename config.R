## =============================================================================
## config.R
## Paths and settings shared by all R scripts. Run every script from the
## repository root, e.g.  Rscript 03_diversity.R
## =============================================================================
DATA_DIR    <- "data"       # input files (see data/README.md)
RESULTS_DIR <- "results"    # all outputs are written here

SEED      <- 2026           # random seed used throughout
DEPTH     <- 40000          # rarefaction depth (reads per sample)
N_RAREFY  <- 100            # number of repeated rarefactions
N_PERM    <- 999            # permutations for PERMANOVA / PERMDISP
N_CORES   <- 2              # parallel::mclapply cores (use 1 on Windows)

data_path <- function(...) file.path(DATA_DIR, ...)
result_path <- function(...) {
  p <- file.path(RESULTS_DIR, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

## Sample metadata: one row per analysed sample (n = 37)
## columns: SampleID, Group (single / multiple), Age, HPV52 (N / Y), Dx (CIN2 / CIN3 / SCC),
##          n_major_hrHPV, CST (I-V), CST_subtype (e.g. III-A, IV-B)
## The row order must be the sample order used in the original analysis
## (permutation P values and SparCC estimates depend on it).
read_metadata <- function() {
  md <- read.csv(data_path("sample_metadata.csv"), stringsAsFactors = FALSE)
  md$Group <- factor(md$Group, c("single", "multiple"))
  md$HPV52 <- factor(md$HPV52, c("N", "Y"))
  md$Age   <- as.numeric(md$Age)
  rownames(md) <- md$SampleID
  md
}

## BIOM-converted TSV tables may start with "# Constructed from biom file"
read_biom_tsv <- function(f, ...) {
  skip <- if (startsWith(readLines(f, n = 1), "# Constructed")) 1 else 0
  read.delim(f, skip = skip, row.names = 1, check.names = FALSE, ...)
}
