## =============================================================================
## 02_differential_abundance.R
## Taxon-level differential abundance (Supplementary Tables S3 and S4)
##
##  1. MaAsLin2 (1.18.0): unadjusted, age-adjusted, and HPV52-adjusted models
##  2. Validation of maaslin_lm() against the MaAsLin2 output
##  3. Sensitivity analyses: minimum abundance 0.1%, exclusion of SCC cases,
##     joint adjustment for age and HPV52, exclusion of HPV52-positive cases
##
## Inputs : data/asv_relative_abundance_1288.tsv  (ASV, Taxon, samples; %)
##          data/top20_species.csv                (species analysed)
##          data/sample_metadata.csv
## =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(Maaslin2) })
source("config.R"); source("R/maaslin_lm.R")
OUT <- file.path(RESULTS_DIR, "taxa"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

meta  <- read_metadata()
abund <- species_table(data_path("asv_relative_abundance_1288.tsv"), data_path("top20_species.csv"))
abund <- as.data.frame(abund[rownames(meta), ])
stopifnot(nrow(abund) == 37, !any(is.na(meta$Age)))

## -----------------------------------------------------------------------------
## 1. MaAsLin2: TSS + LOG + LM, min_prevalence 0.1, reference = single / HPV52 N
## -----------------------------------------------------------------------------
run_maaslin <- function(fixed, ref, tag) {
  fit <- Maaslin2(input_data = abund, input_metadata = meta,
                  output = file.path(OUT, paste0("maaslin2_", tag)),
                  fixed_effects = fixed, reference = ref,
                  normalization = "TSS", transform = "LOG", analysis_method = "LM",
                  standardize = FALSE, min_prevalence = 0.1, min_abundance = 0,
                  max_significance = 0.1, plot_heatmap = FALSE, plot_scatter = FALSE, cores = 1)
  fit$results %>% mutate(model = tag)
}
res <- bind_rows(
  run_maaslin("Group",             "Group,single",             "unadjusted"),
  run_maaslin(c("Group", "Age"),   "Group,single",             "age_adjusted"),
  run_maaslin(c("Group", "HPV52"), c("Group,single", "HPV52,N"), "HPV52_adjusted"))
write.csv(res, result_path("taxa", "maaslin2_taxa_all_models.csv"), row.names = FALSE)

## Supplementary Table S3: group term and covariate terms, taxa with P < 0.05 in any model
grp <- res %>% filter(metadata == "Group") %>% select(feature, model, coef, stderr, pval, qval)
sig_any <- grp %>% filter(pval < 0.05) %>% distinct(feature) %>% pull(feature)
fmt <- function(x, d = 4) formatC(x, digits = d, format = "f")
wide <- function(models, cov_term) {
  g <- grp %>% filter(model %in% models, feature %in% sig_any) %>%
    mutate(`Coef (SE)` = paste0(fmt(coef, 3), " (", fmt(stderr, 3), ")"), `P (Q)` = paste0(fmt(pval), " (", fmt(qval), ")")) %>%
    select(feature, model, `Coef (SE)`, `P (Q)`) %>%
    pivot_wider(names_from = model, values_from = c(`Coef (SE)`, `P (Q)`), names_glue = "{model}: {.value}")
  cv <- res %>% filter(metadata == cov_term, feature %in% sig_any) %>%
    transmute(feature, `Covariate: Coef (SE)` = paste0(fmt(coef, 3), " (", fmt(stderr, 3), ")"),
              `Covariate: P (Q)` = paste0(fmt(pval), " (", fmt(qval), ")"))
  g %>% left_join(cv, by = "feature") %>% arrange(feature)
}
write.csv(wide(c("unadjusted", "age_adjusted"), "Age"),     result_path("taxa", "taxa_unadjusted_vs_age.csv"),   row.names = FALSE)
write.csv(wide(c("unadjusted", "HPV52_adjusted"), "HPV52"), result_path("taxa", "taxa_unadjusted_vs_HPV52.csv"), row.names = FALSE)

## -----------------------------------------------------------------------------
## 2. Validation: maaslin_lm() vs MaAsLin2 (differences should be ~0)
## -----------------------------------------------------------------------------
orig <- res; orig$feature <- gsub("\\.", " ", orig$feature)   # MaAsLin2 replaces spaces in feature names by "."
for (m in list(list("unadjusted", "Group"), list("age_adjusted", c("Group", "Age")), list("HPV52_adjusted", c("Group", "HPV52")))) {
  r <- maaslin_lm(abund, meta, m[[2]]); o <- orig[orig$model == m[[1]], ]
  k <- merge(r, o, by = c("feature", "name"), suffixes = c(".lm", ".maaslin2"))
  cat(sprintf("%-15s tests %d/%d | max|coef diff| %.2e | max|P diff| %.2e | max|q diff| %.2e\n", m[[1]], nrow(r), nrow(o),
              max(abs(k$coef.lm - k$coef.maaslin2)), max(abs(k$pval.lm - k$pval.maaslin2)), max(abs(k$qval.lm - k$qval.maaslin2))))
}

## -----------------------------------------------------------------------------
## 3. Sensitivity analyses (Supplementary Table S4)
## -----------------------------------------------------------------------------
noscc  <- rownames(meta)[meta$Dx != "SCC"]       # three SCC cases excluded (n = 34)
keep52 <- rownames(meta)[meta$HPV52 == "N"]      # five HPV52-positive cases excluded (n = 32)
sens <- rbind(
  cbind(Model = "Minimum abundance 0.1%, unadjusted (n = 37)",          maaslin_lm(abund, meta, "Group", min_abund = 0.1)),
  cbind(Model = "Minimum abundance 0.1%, age + HPV52 (n = 37)",         maaslin_lm(abund, meta, c("Group", "Age", "HPV52"), min_abund = 0.1)),
  cbind(Model = "SCC excluded, unadjusted (n = 34)",                    maaslin_lm(abund[noscc, ], meta[noscc, ], "Group")),
  cbind(Model = "SCC excluded, age + HPV52 (n = 34)",                   maaslin_lm(abund[noscc, ], meta[noscc, ], c("Group", "Age", "HPV52"))),
  cbind(Model = "Adjusted for age and HPV52 (n = 37)",                  maaslin_lm(abund, meta, c("Group", "Age", "HPV52"))),
  cbind(Model = "HPV52-negative only, unadjusted (n = 32)",             maaslin_lm(abund[keep52, ], meta[keep52, ], "Group")),
  cbind(Model = "HPV52-negative only, age-adjusted (n = 32)",           maaslin_lm(abund[keep52, ], meta[keep52, ], c("Group", "Age"))))
write.csv(sens, result_path("taxa", "taxa_sensitivity_analyses.csv"), row.names = FALSE)
print(sens[sens$name == "Groupmultiple" & sens$feature == "Segatella copri", ], digits = 3)
