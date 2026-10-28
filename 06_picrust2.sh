#!/usr/bin/env bash
# PICRUSt2 (default settings, --max_nsti 2): all 1,288 ASVs, and the 230 ASVs
# assigned to the 20 most abundant species
set -euo pipefail

picrust2_pipeline.py -s data/rep_seqs/dna-sequences.fasta -i data/feature_table_1288.biom \
  -o picrust2_out_1288 -p 8 --stratified
add_descriptions.py -i picrust2_out_1288/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
  -o picrust2_out_1288/pathways_out/path_abun_unstrat_descrip.tsv.gz

mkdir -p picrust2_in_230
Rscript -e '
ds   <- read.csv("data/asv_table_1288.csv", check.names = FALSE)
md   <- read.csv("data/sample_metadata.csv")
keep <- ds[ds$Species %in% read.csv("data/network_species.csv")$Taxa, c("ASV", md$SampleID)]
fa   <- readLines("data/rep_seqs/dna-sequences.fasta")
seqs <- setNames(fa[c(FALSE, TRUE)], sub(">", "", fa[c(TRUE, FALSE)]))
writeLines(paste0(">", keep$ASV, "\n", seqs[keep$ASV]), "picrust2_in_230/seqs.fasta")
names(keep)[1] <- "#OTU ID"
write.table(keep, "picrust2_in_230/table.tsv", sep = "\t", quote = FALSE, row.names = FALSE)'
biom convert -i picrust2_in_230/table.tsv -o picrust2_in_230/table.biom --table-type="OTU table" --to-hdf5

picrust2_pipeline.py -s picrust2_in_230/seqs.fasta -i picrust2_in_230/table.biom \
  -o picrust2_out_230 -p 8 --stratified
