Analysis code for:
Vaginal Microbiome Alterations Associated with Multiple High-Risk Human Papillomavirus Infections in Women with High-Grade Cervical Lesions. 
Raw 16S rRNA gene sequencing reads: NCBI SRA BioProject **PRJNA1537372**.

## Scripts (run in order from the repository root)

| Script | Content |
|---|---|
| `01_qiime2_processing.sh` | Trimming, DADA2, ASV filtering, taxonomy, phylogeny, core metrics (QIIME 2 2024.02) |
| `02_differential_abundance.R` | Taxon-level MaAsLin2 and sensitivity analyses |
| `03_diversity.R` | Alpha and beta diversity, 100 rarefactions, rare-ASV filters, rarefaction curves |
| `03b_diversity_sensitivity.R` | Regression models, PERMANOVA robustness, covariate adjustment, exclusion analyses |
| `04_network.R` | SparCC networks (NetCoMi 1.3.0), thresholds, bootstrap, subsampling, HPV52 exclusion |
| `05_network_functions.R` | Functions used by `04_network.R` for resampling |
| `06_picrust2.sh` | PICRUSt2 2.6.3 for all ASVs and for the 230 ASVs of the 20 most abundant species |
| `07_pathway_analysis.R` | Pathway models, taxon contributions, NSTI, pathway detection |

`config.R` holds paths and settings (seed 2026); `R/maaslin_lm.R` fits linear models with the
MaAsLin2 settings used for the sensitivity analyses. The scripts print a comparison of these
implementations with the MaAsLin2, NetCoMi, and QIIME 2 output.

## License

Code: [MIT License]. Please cite the article when using this code.
