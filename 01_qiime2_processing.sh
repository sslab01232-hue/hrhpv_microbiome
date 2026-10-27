#!/usr/bin/env bash
###############################################################################
# 01_qiime2_processing.sh
#
# Read QC -> primer removal -> quality trimming -> DADA2 -> taxonomy
# -> sample / ASV filtering (1,288 ASVs, 37 samples) -> phylogenetic tree
# -> core-metrics-phylogenetic (40,000 reads) -> export for R
#
# Software: FastQC 0.12.1, Cutadapt 5.0, Trimmomatic 0.39, QIIME 2 2024.02
#           (q2-cutadapt, q2-dada2, q2-feature-table, q2-phylogeny,
#            q2-diversity), SILVA 138
###############################################################################
set -euo pipefail

THREADS=4
RAW_DIR=raw_fastq            # demultiplexed paired-end FASTQ
OUT=qiime2_out
mkdir -p ${OUT} fastqc_out trimmed_fastq

# -----------------------------------------------------------------------------
# 1. Read quality assessment - FastQC 0.12.1, default settings
# -----------------------------------------------------------------------------
fastqc -t ${THREADS} -o fastqc_out ${RAW_DIR}/*.fastq.gz

# -----------------------------------------------------------------------------
# 2. Quality trimming - Trimmomatic 0.39 (SLIDINGWINDOW:4:15, MINLEN:36)
# -----------------------------------------------------------------------------
for R1 in ${RAW_DIR}/*_R1*.fastq.gz; do
  R2=${R1/_R1/_R2}
  s=$(basename ${R1} | sed 's/_R1.*//')
  trimmomatic PE -threads ${THREADS} -phred33 \
    ${R1} ${R2} \
    trimmed_fastq/${s}_R1.fastq.gz trimmed_fastq/${s}_R1.unpaired.fastq.gz \
    trimmed_fastq/${s}_R2.fastq.gz trimmed_fastq/${s}_R2.unpaired.fastq.gz \
    SLIDINGWINDOW:4:15 MINLEN:36
done

# -----------------------------------------------------------------------------
# 3. Import (manifest of the sequenced samples)
# -----------------------------------------------------------------------------
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path pe-64-manifest \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path ${OUT}/demux-paired-end.qza

# -----------------------------------------------------------------------------
# 4. Primer removal - q2-cutadapt trim-paired
#    341F CCTACGGGNGGCWGCAG / 805R GACTACHVGGGTATCTAATCC
# -----------------------------------------------------------------------------
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences ${OUT}/demux-paired-end.qza \
  --p-cores 4 \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-error-rate 0.1 \
  --p-indels \
  --p-times 1 \
  --p-overlap 3 \
  --p-match-adapter-wildcards \
  --p-no-match-read-wildcards \
  --p-minimum-length 1 \
  --p-no-discard-untrimmed \
  --p-quality-cutoff-5end 0 \
  --p-quality-cutoff-3end 0 \
  --p-quality-base 33 \
  --o-trimmed-sequences ${OUT}/demux-trimmed.qza

qiime demux summarize \
  --i-data ${OUT}/demux-trimmed.qza \
  --o-visualization ${OUT}/demux-trimmed.qzv

# -----------------------------------------------------------------------------
# 5. Denoising - DADA2 denoise-paired
# -----------------------------------------------------------------------------
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ${OUT}/demux-trimmed.qza \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 220 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-max-ee-f 5 \
  --p-max-ee-r 5 \
  --p-trunc-q 2 \
  --p-min-overlap 12 \
  --p-pooling-method independent \
  --p-chimera-method consensus \
  --p-min-fold-parent-over-abundance 1.0 \
  --p-n-reads-learn 1000000 \
  --p-hashed-feature-ids \
  --p-n-threads 1 \
  --o-table ${OUT}/table.qza \
  --o-representative-sequences ${OUT}/rep-seqs.qza \
  --o-denoising-stats ${OUT}/denoising-stats.qza

# -----------------------------------------------------------------------------
# 6. Taxonomic assignment - naive Bayes classifier trained on SILVA 138
#    (species-level names were subsequently curated; see step 9)
# -----------------------------------------------------------------------------
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads ${OUT}/rep-seqs.qza \
  --o-classification ${OUT}/taxonomy.qza

# -----------------------------------------------------------------------------
# 7. Sample filtering -> 37 analysed samples
#    metadata.tsv lists the 40 eligible samples (Filter = Keep 37 / Remove 3;
#    the three samples failed DNA extraction or library preparation)
# -----------------------------------------------------------------------------
qiime feature-table filter-samples \
  --i-table ${OUT}/table.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[Filter]='Keep'" \
  --p-filter-empty-features \
  --o-filtered-table ${OUT}/table.37cases.qza

qiime feature-table filter-seqs \
  --i-data ${OUT}/rep-seqs.qza \
  --i-table ${OUT}/table.37cases.qza \
  --o-filtered-data ${OUT}/rep-seqs.37cases.qza

# -----------------------------------------------------------------------------
# 8. ASV filtering - ASVs with a total of < 10 reads or a relative abundance
#    < 0.01% in every sample are excluded (1,293 ASVs retained; 5 further ASVs
#    are removed with the curated taxonomy in step 9 -> 1,288 ASVs)
# -----------------------------------------------------------------------------
qiime tools export --input-path ${OUT}/table.37cases.qza --output-path ${OUT}/export_37
biom convert -i ${OUT}/export_37/feature-table.biom \
             -o ${OUT}/export_37/feature-table.tsv --to-tsv

Rscript - <<'REOF'
x  <- read.delim("qiime2_out/export_37/feature-table.tsv", skip = 1,
                 check.names = FALSE, row.names = 1)
ra <- sweep(x, 2, colSums(x), "/") * 100            # % relative abundance
keep <- rowSums(x) >= 10 & apply(ra, 1, max) >= 0.01
out <- data.frame(`feature-id` = rownames(x)[keep], check.names = FALSE)
write.table(out, "ASV_keep_list.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("ASVs retained:", sum(keep), "\n")              # 1,293
REOF

qiime feature-table filter-features \
  --i-table ${OUT}/table.37cases.qza \
  --m-metadata-file ASV_keep_list.tsv \
  --o-filtered-table ${OUT}/table.37cases.filtered.qza

qiime feature-table filter-seqs \
  --i-data ${OUT}/rep-seqs.37cases.qza \
  --i-table ${OUT}/table.37cases.filtered.qza \
  --o-filtered-data ${OUT}/rep-seqs.37cases.filtered.qza

# -----------------------------------------------------------------------------
# 9. Curated taxonomy (species names confirmed by BLAST; current SILVA/LPSN
#    nomenclature) -> import, then restrict table and sequences to it
# -----------------------------------------------------------------------------
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format TSVTaxonomyFormat \
  --input-path taxonomy.filtered.final.no2.tsv \
  --output-path ${OUT}/taxonomy.filtered.final.no2.qza

qiime feature-table filter-features \
  --i-table ${OUT}/table.37cases.filtered.qza \
  --m-metadata-file ${OUT}/taxonomy.filtered.final.no2.qza \
  --o-filtered-table ${OUT}/table.37cases.filtered.no2.qza

qiime feature-table filter-seqs \
  --i-data ${OUT}/rep-seqs.37cases.filtered.qza \
  --i-table ${OUT}/table.37cases.filtered.no2.qza \
  --o-filtered-data ${OUT}/rep-seqs.37cases.filtered.no2.qza

qiime feature-table summarize \
  --i-table ${OUT}/table.37cases.filtered.no2.qza \
  --m-sample-metadata-file metadata.tsv \
  --o-visualization ${OUT}/table.37cases.filtered.no2.qzv   # 1,288 ASVs x 37 samples

qiime feature-table relative-frequency \
  --i-table ${OUT}/table.37cases.filtered.no2.qza \
  --o-relative-frequency-table ${OUT}/table.37cases.filtered.no2.RA.qza

# -----------------------------------------------------------------------------
# 10. Phylogenetic tree - align-to-tree-mafft-fasttree
#     MAFFT -> mask (max-gap-frequency 1.0, min-conservation 0.4)
#     -> FastTree -> midpoint root
# -----------------------------------------------------------------------------
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${OUT}/rep-seqs.37cases.filtered.no2.qza \
  --p-mask-max-gap-frequency 1.0 \
  --p-mask-min-conservation 0.4 \
  --p-n-threads 1 \
  --o-alignment ${OUT}/aligned-rep-seqs.no2.qza \
  --o-masked-alignment ${OUT}/masked-aligned-rep-seqs.no2.qza \
  --o-tree ${OUT}/unrooted-tree.no2.qza \
  --o-rooted-tree ${OUT}/rooted-tree.no2.qza

# -----------------------------------------------------------------------------
# 11. Rarefaction (40,000 reads, without replacement) + alpha / beta metrics
#     core-metrics-phylogenetic
#     -> observed features, Shannon, Pielou, Faith's PD, Jaccard, Bray-Curtis,
#        unweighted / weighted UniFrac, PCoA
# -----------------------------------------------------------------------------
qiime diversity core-metrics-phylogenetic \
  --i-table ${OUT}/table.37cases.filtered.no2.qza \
  --i-phylogeny ${OUT}/rooted-tree.no2.qza \
  --p-sampling-depth 40000 \
  --m-metadata-file metadata.tsv \
  --p-n-jobs-or-threads 1 \
  --output-dir core-metrics-result-HM-40000-no2

# PERMDISP (unweighted UniFrac, 999 permutations)
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-result-HM-40000-no2/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column HPV_gr1 \
  --p-method permdisp \
  --p-permutations 999 \
  --p-no-pairwise \
  --o-visualization core-metrics-result-HM-40000-no2/unweighted_unifrac_distance_matrix_permdisp.qzv

# -----------------------------------------------------------------------------
# 12. Export for R (02-07) and PICRUSt2 (06_picrust2.sh) -> data/
# -----------------------------------------------------------------------------
mkdir -p data/core_metrics data/rep_seqs
qiime tools export --input-path ${OUT}/table.37cases.filtered.no2.qza    --output-path data/table
biom convert -i data/table/feature-table.biom -o data/feature_table_1288.tsv --to-tsv
cp data/table/feature-table.biom data/feature_table_1288.biom
qiime tools export --input-path ${OUT}/rep-seqs.37cases.filtered.no2.qza --output-path data/rep_seqs
qiime tools export --input-path ${OUT}/taxonomy.filtered.final.no2.qza   --output-path data/taxonomy
qiime tools export --input-path ${OUT}/rooted-tree.no2.qza               --output-path data/tree_export
cp data/tree_export/tree.nwk data/tree.nwk
for m in observed_features_vector shannon_vector evenness_vector faith_pd_vector \
         jaccard_distance_matrix bray_curtis_distance_matrix \
         unweighted_unifrac_distance_matrix weighted_unifrac_distance_matrix rarefied_table; do
  qiime tools export --input-path core-metrics-result-HM-40000-no2/${m}.qza \
                     --output-path data/core_metrics/${m}
done
biom convert -i data/core_metrics/rarefied_table/feature-table.biom \
             -o data/rarefied_table_40000.tsv --to-tsv
