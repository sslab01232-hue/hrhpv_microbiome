#!/usr/bin/env bash
###############################################################################
# 06_picrust2.sh
# PICRUSt2 functional prediction
#   1. all 1,288 ASVs of the final feature table (primary analysis)
#   2. the 230 ASVs assigned to the 20 most abundant species (comparison)
#
# Run inside the PICRUSt2 conda environment (conda activate picrust2).
# Default settings: --max_nsti 2, EPA-ng placement + gappa, HSP method mp,
# MinPath pathway inference, default reference files.
# Inputs (from 01_qiime2_processing.sh):
#   data/feature_table_1288.biom        raw counts, not rarefied
#   data/rep_seqs/dna-sequences.fasta   representative sequences
#   data/asv_table_1288.csv             ASV table with taxonomy (Species column)
#   data/network_species.csv            the 20 most abundant species (column Taxa)
###############################################################################
set -euo pipefail
THREADS=8

# -----------------------------------------------------------------------------
# 1. All 1,288 ASVs
# -----------------------------------------------------------------------------
picrust2_pipeline.py \
  -s data/rep_seqs/dna-sequences.fasta \
  -i data/feature_table_1288.biom \
  -o picrust2_out_1288 \
  -p ${THREADS} \
  --stratified \
  --verbose
add_descriptions.py -i picrust2_out_1288/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o picrust2_out_1288/pathways_out/path_abun_unstrat_descrip.tsv.gz

# -----------------------------------------------------------------------------
# 2. 230 ASVs assigned to the 20 most abundant species
# -----------------------------------------------------------------------------
mkdir -p picrust2_input_230
Rscript - <<'REOF'
ds   <- read.csv("data/asv_table_1288.csv", check.names = FALSE)
md   <- read.csv("data/sample_metadata.csv")
sp20 <- read.csv("data/network_species.csv", check.names = FALSE)$Taxa
keep <- ds[ds$Species %in% sp20, ]; stopifnot(nrow(keep) == 230)
fa   <- readLines("data/rep_seqs/dna-sequences.fasta")
hdr  <- grep("^>", fa); ids <- sub("^>", "", fa[hdr])
seqs <- setNames(sapply(seq_along(hdr), function(k)
          paste(fa[(hdr[k] + 1):(if (k < length(hdr)) hdr[k + 1] - 1 else length(fa))], collapse = "")), ids)
stopifnot(all(keep$ASV %in% names(seqs)))
tab <- keep[, c("ASV", md$SampleID)]; names(tab)[1] <- "#OTU ID"
write.table(tab, "picrust2_input_230/feature-table-230.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(as.vector(rbind(paste0(">", keep$ASV), seqs[keep$ASV])), "picrust2_input_230/dna-sequences-230.fasta")
cat("ASVs:", nrow(tab), "| reads:", sum(tab[, -1]), "\n")
REOF

biom convert -i picrust2_input_230/feature-table-230.tsv -o picrust2_input_230/feature-table-230.biom \
             --table-type="OTU table" --to-hdf5
picrust2_pipeline.py \
  -s picrust2_input_230/dna-sequences-230.fasta \
  -i picrust2_input_230/feature-table-230.biom \
  -o picrust2_out_230 \
  -p ${THREADS} \
  --stratified \
  --verbose
add_descriptions.py -i picrust2_out_230/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o picrust2_out_230/pathways_out/path_abun_unstrat_descrip.tsv.gz
