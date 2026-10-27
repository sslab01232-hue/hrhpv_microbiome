Analysis code for:
Vaginal Microbiome Alterations Associated with Multiple High-Risk Human Papillomavirus Infections in Women with High-Grade Cervical Lesions. 
Raw 16S rRNA gene sequencing reads: NCBI SRA BioProject **PRJNA1537372**.

## Scripts (run in order from the repository root)

| Script | Content | Software |
|---|---|---|
| `01_qiime2_processing.sh` | Read QC, primer removal, DADA2, taxonomy, sample and ASV filtering (1,288 ASVs, 37 samples), phylogenetic tree, core-metrics (40,000 reads), export to `data/` | FastQC 0.12.1, Trimmomatic 0.39, QIIME 2 2024.02 (q2-cutadapt, q2-dada2, q2-phylogeny, q2-diversity), SILVA 138 |
| `02_differential_abundance.R` | Taxon-level differential abundance and sensitivity analyses | MaAsLin2 1.18.0 |
| `03_diversity.R` | Read depth, rarefaction curves, alpha and beta diversity, 100 rarefactions, rare-ASV filters | vegan, ape |
| `03b_diversity_tests_and_sensitivity.R` | Alpha regression models, PERMANOVA robustness (subsampling range, PERMDISP, ANOSIM, Mantel, MRM), covariate-adjusted PERMANOVA, three-group analysis, HPV52 and SCC exclusion | vegan, ecodist |
| `04_network.R` | Co-occurrence networks, network comparison, edge lists, sensitivity analyses (thresholds, bootstrap, subsampling, HPV52 exclusion) | NetCoMi 1.3.0, SpiecEasi |
| `05_network_functions.R` | Functions used by `04_network.R` (sourced) | igraph, VGAM |
| `06_picrust2.sh` | PICRUSt2 for all 1,288 ASVs and for the 230 ASVs of the 20 most abundant species | PICRUSt2 2.6.3 |
| `07_pathway_analysis.R` | Pathway models (age and richness adjustment), taxon contributions, NSTI, pathway detection | — |
| `config.R`, `R/maaslin_lm.R` | Shared paths and settings; MaAsLin2-equivalent linear models | — |

All outputs are written to `results/`. The random seed is 2026 throughout.

## Notes

- `R/maaslin_lm.R` reproduces the MaAsLin2 settings used in this study and was used for the sensitivity models; its agreement with MaAsLin2 is checked in `02_differential_abundance.R` (section 2).
- `05_network_functions.R` transcribes the NetCoMi/SpiecEasi functions so that networks can be rebuilt under resampling; its agreement with NetCoMi is checked in `04_network.R` (section 4).
- Permutation P values and SparCC estimates depend on the sample order; keep the sample order of `sample_metadata.csv` and of the QIIME 2 artifacts.
- Descriptive comparisons of clinical characteristics (Table 1), CST and *Lactobacillus*-dominance distributions, and *Lactobacillus* relative abundance used standard tests (`wilcox.test`, `fisher.test`) as described in the Materials and Methods and can be reproduced from `sample_metadata.csv` and the ASV table.
- PERMANOVA and PERMDISP use 999 permutations, except for the analysis excluding SCC cases (9,999 permutations).
- `parallel::mclapply` is used for the resampling analyses; set `N_CORES <- 1` in `config.R` on Windows.

## License

Code: [MIT License]. Please cite the article when using this code.
