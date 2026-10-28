#!/usr/bin/env bash
# Read processing, ASV filtering, phylogeny, and diversity metrics
# FastQC 0.12.1, Trimmomatic 0.39, QIIME 2 2024.02, SILVA 138
set -euo pipefail
OUT=qiime2_out
mkdir -p $OUT fastqc_out trimmed_fastq data/core_metrics

# Quality check and trimming
fastqc -t 4 -o fastqc_out raw_fastq/*.fastq.gz
for r1 in raw_fastq/*_R1*.fastq.gz; do
  s=$(basename $r1 | sed 's/_R1.*//')
  trimmomatic PE -threads 4 -phred33 $r1 ${r1/_R1/_R2} \
    trimmed_fastq/${s}_R1.fastq.gz /dev/null trimmed_fastq/${s}_R2.fastq.gz /dev/null \
    SLIDINGWINDOW:4:15 MINLEN:36
done

# Import, primer removal (341F/805R), and DADA2
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path pe-64-manifest --input-format PairedEndFastqManifestPhred33V2 \
  --output-path $OUT/demux.qza
qiime cutadapt trim-paired --i-demultiplexed-sequences $OUT/demux.qza \
  --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACHVGGGTATCTAATCC \
  --p-error-rate 0.1 --p-overlap 3 --p-cores 4 \
  --o-trimmed-sequences $OUT/demux-trimmed.qza
qiime dada2 denoise-paired --i-demultiplexed-seqs $OUT/demux-trimmed.qza \
  --p-trunc-len-f 280 --p-trunc-len-r 220 --p-max-ee-f 5 --p-max-ee-r 5 \
  --p-chimera-method consensus --p-pooling-method independent \
  --o-table $OUT/table.qza --o-representative-sequences $OUT/rep-seqs.qza \
  --o-denoising-stats $OUT/denoising-stats.qza

# Taxonomy (naive Bayes, SILVA 138; species names curated afterwards)
qiime feature-classifier classify-sklearn --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads $OUT/rep-seqs.qza --o-classification $OUT/taxonomy.qza

# Samples: 37 of the 40 eligible (Filter column in metadata.tsv)
qiime feature-table filter-samples --i-table $OUT/table.qza --m-metadata-file metadata.tsv \
  --p-where "[Filter]='Keep'" --p-filter-empty-features --o-filtered-table $OUT/table37.qza

# ASVs: >= 10 reads in total and >= 0.01% relative abundance in at least one sample (1,288 ASVs)
qiime tools export --input-path $OUT/table37.qza --output-path $OUT/table37
biom convert -i $OUT/table37/feature-table.biom -o $OUT/table37/table.tsv --to-tsv
Rscript -e '
x <- read.delim("qiime2_out/table37/table.tsv", skip = 1, row.names = 1, check.names = FALSE)
keep <- rowSums(x) >= 10 & apply(t(t(x) / colSums(x)), 1, max) >= 1e-4
write.table(data.frame("feature-id" = rownames(x)[keep], check.names = FALSE), "asv_keep.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)'
qiime feature-table filter-features --i-table $OUT/table37.qza --m-metadata-file asv_keep.tsv \
  --o-filtered-table $OUT/table-final.qza
qiime feature-table filter-seqs --i-data $OUT/rep-seqs.qza --i-table $OUT/table-final.qza \
  --o-filtered-data $OUT/rep-seqs-final.qza

# Curated taxonomy (species names checked by BLAST)
qiime tools import --type 'FeatureData[Taxonomy]' --input-format TSVTaxonomyFormat \
  --input-path taxonomy_curated.tsv --output-path $OUT/taxonomy-final.qza

# Phylogeny (MAFFT, FastTree, midpoint root; tree pruned to the final ASVs)
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $OUT/rep-seqs.qza \
  --o-alignment $OUT/aligned.qza --o-masked-alignment $OUT/masked.qza \
  --o-tree $OUT/unrooted-tree-all.qza --o-rooted-tree $OUT/rooted-tree-all.qza
qiime phylogeny filter-tree --i-tree $OUT/rooted-tree-all.qza --i-table $OUT/table-final.qza \
  --o-filtered-tree $OUT/rooted-tree.qza

# Core metrics (rarefied to 40,000 reads)
qiime diversity core-metrics-phylogenetic --i-table $OUT/table-final.qza \
  --i-phylogeny $OUT/rooted-tree.qza --p-sampling-depth 40000 \
  --m-metadata-file metadata.tsv --output-dir core-metrics

# Export for the R scripts
qiime tools export --input-path $OUT/table-final.qza    --output-path data
qiime tools export --input-path $OUT/rep-seqs-final.qza --output-path data/rep_seqs
qiime tools export --input-path $OUT/rooted-tree.qza    --output-path data
qiime tools export --input-path $OUT/taxonomy-final.qza --output-path data/taxonomy
biom convert -i data/feature-table.biom -o data/feature_table_1288.tsv --to-tsv
mv data/feature-table.biom data/feature_table_1288.biom
ls core-metrics/*_vector.qza core-metrics/*_distance_matrix.qza core-metrics/rarefied_table.qza |
  xargs -I{} sh -c 'qiime tools export --input-path {} --output-path data/core_metrics/$(basename {} .qza)'
biom convert -i data/core_metrics/rarefied_table/feature-table.biom -o data/rarefied_table_40000.tsv --to-tsv
