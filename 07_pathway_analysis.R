## =============================================================================
## 07_pathway_analysis.R
## Predicted MetaCyc pathways (Supplementary Tables S27-S29; Figure 4 data)
##
##  1. Differential abundance of pathways, PICRUSt2 on all 1,288 ASVs:
##     unadjusted, age-adjusted (primary), and age + log(observed features)
##  2. Comparison with the analysis restricted to the 230 ASVs of the 20 most
##     abundant species
##  3. Taxon-level contributions to the enriched pathways
##  4. NSTI
##  5. Richness diagnostics and pathway detection
##
## Models use maaslin_lm() (R/maaslin_lm.R; MaAsLin2 settings: TSS, LOG, LM,
## min_prevalence 0.1, standardize = FALSE, BH over all feature x term
## coefficients of each model; significance q < 0.1).
##
## Inputs : picrust2_out_1288/ and picrust2_out_230/   (06_picrust2.sh)
##          data/asv_relative_abundance_1288.tsv, data/sample_metadata.csv,
##          data/core_metrics/observed_features_vector/alpha-diversity.tsv
## =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
source("config.R"); source("R/maaslin_lm.R")
OUT <- file.path(RESULTS_DIR, "pathways"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
P1288 <- "picrust2_out_1288"; P230 <- "picrust2_out_230"

md  <- read_metadata()
obs <- read.delim(data_path("core_metrics", "observed_features_vector", "alpha-diversity.tsv"), row.names = 1)
md$observed_features <- obs[rownames(md), 1]                      # 40,000-read rarefied table
md$logRich <- log(md$observed_features)

read_pw <- function(dir) { u <- read.delim(file.path(dir, "pathways_out", "path_abun_unstrat.tsv.gz"), check.names = FALSE, row.names = 1)
  stopifnot(setequal(colnames(u), md$SampleID)); t(as.matrix(u[, rownames(md)])) }
ab <- read_pw(P1288)
de <- read.delim(file.path(P1288, "pathways_out", "path_abun_unstrat_descrip.tsv.gz"), check.names = FALSE)
dm <- setNames(de$description, de[[1]])
ent <- c("&alpha;" = "α", "&beta;" = "β", "&gamma;" = "γ", "&delta;" = "δ", "&epsilon;" = "ε", "&omega;" = "ω", "&prime;" = "′",
         "&rarr;" = "→", "&harr;" = "↔", "&amp;" = "&", "&lt;" = "<", "&gt;" = ">", "&quot;" = "\"", "&apos;" = "'")
for (e in names(ent)) dm <- gsub(e, ent[[e]], dm, fixed = TRUE)

## =============================================================================
## 1. Pathway models (1,288 ASVs)
## =============================================================================
models <- list(unadjusted = "Group", age_adjusted = c("Group", "Age"), age_logrichness_adjusted = c("Group", "Age", "logRich"))
res <- bind_rows(lapply(names(models), function(m) maaslin_lm(ab, md, models[[m]]) %>% mutate(model = m))) %>%
  mutate(description = unname(dm[feature]))
write.csv(res, file.path(OUT, "pathway_all_models_1288.csv"), row.names = FALSE)

grp <- res %>% filter(name == "Groupmultiple")
det <- data.frame(feature = colnames(ab), det_single = colSums(ab[md$Group == "single", ] > 0), det_multiple = colSums(ab[md$Group == "multiple", ] > 0))
summ <- grp %>% group_by(model) %>%
  summarise(n_tested = n(), pct_coef_positive = 100 * mean(coef > 0), n_q01_multiple = sum(qval < 0.1 & coef > 0),
            n_q01_single = sum(qval < 0.1 & coef < 0), n_q005 = sum(qval < 0.05), .groups = "drop") %>% slice(match(names(models), model))
rich <- res %>% filter(model == "age_logrichness_adjusted", name == "logRich")
summ$logrich_pct_positive <- c(NA, NA, 100 * mean(rich$coef > 0)); summ$logrich_n_q01 <- c(NA, NA, sum(rich$qval < 0.1))
write.csv(summ, file.path(OUT, "summary_by_model.csv"), row.names = FALSE); print(as.data.frame(summ), digits = 3)

wide <- grp %>% select(model, feature, description, coef, stderr, pval, qval) %>%
  pivot_wider(names_from = model, values_from = c(coef, stderr, pval, qval)) %>% left_join(det, by = "feature") %>% arrange(qval_age_adjusted)
write.csv(wide, file.path(OUT, "pathway_group_term_wide_1288.csv"), row.names = FALSE)
sigA <- grp %>% filter(model == "age_adjusted", qval < 0.1) %>% arrange(qval)                       # primary model
sigR <- grp %>% filter(model == "age_logrichness_adjusted", qval < 0.1, coef > 0)
six  <- intersect(sigA$feature[sigA$coef > 0], sigR$feature)                                         # enriched in both models
write.csv(wide %>% filter(feature %in% sigA$feature), file.path(OUT, "pathways_q01_age_adjusted.csv"), row.names = FALSE)
cat("age-adjusted q < 0.1:", nrow(sigA), "| also q < 0.1 after richness adjustment:", length(six), "\n")

## =============================================================================
## 2. Analysis restricted to the 230 ASVs of the 20 most abundant species
## =============================================================================
ab230  <- read_pw(P230)
res230 <- bind_rows(lapply(names(models)[1:2], function(m) maaslin_lm(ab230, md, models[[m]]) %>% mutate(model = m)))
write.csv(res230, file.path(OUT, "pathway_all_models_230.csv"), row.names = FALSE)
s230 <- res230 %>% filter(model == "age_adjusted", name == "Groupmultiple", qval < 0.1)
cmp  <- s230 %>% transmute(pathway = feature, coef230 = coef, q230 = qval) %>%
  left_join(wide %>% transmute(pathway = feature, coef1288 = coef_age_adjusted, q1288 = qval_age_adjusted,
                               q1288_richness = qval_age_logrichness_adjusted), by = "pathway")
write.csv(cmp, file.path(OUT, "comparison_230_vs_1288.csv"), row.names = FALSE); print(cmp, digits = 3)

## =============================================================================
## 3. Taxon-level contributions to the age-adjusted significant pathways
##    (contributions of ASVs with the same taxonomic label summed; mean % of
##     pathway abundance in samples in which the pathway was detected)
## =============================================================================
asv <- read.delim(data_path("asv_relative_abundance_1288.tsv"), check.names = FALSE); names(asv)[1:2] <- c("ASV", "Taxon")
lab <- setNames(taxon_label(asv$Taxon), asv$ASV)
con <- read.delim(file.path(P1288, "pathways_out", "path_abun_contrib.tsv.gz"), check.names = FALSE) %>% filter(`function` %in% sigA$feature)
cs_obs <- con %>% group_by(sample, `function`) %>% mutate(norm100 = taxon_function_abun / sum(taxon_function_abun) * 100) %>% ungroup() %>%
  mutate(taxon_label = unname(lab[taxon])) %>% group_by(sample, `function`, taxon_label) %>% summarise(norm100 = sum(norm100), .groups = "drop")
stopifnot(!any(is.na(cs_obs$taxon_label)))
pres <- as.data.frame(as.table(t(ab[, sigA$feature, drop = FALSE] > 0))) %>% setNames(c("function", "sample", "det")) %>%
  filter(det) %>% mutate(across(c(`function`, sample), as.character))
cs <- pres %>% select(-det) %>%
  inner_join(expand.grid(`function` = sigA$feature, taxon_label = unique(cs_obs$taxon_label), stringsAsFactors = FALSE),
             by = "function", relationship = "many-to-many") %>%
  left_join(cs_obs, by = c("sample", "function", "taxon_label")) %>% mutate(norm100 = coalesce(norm100, 0)) %>%
  left_join(md %>% select(SampleID, Group), by = c("sample" = "SampleID")) %>%
  group_by(`function`, taxon_label, Group) %>% summarise(mean_contrib = mean(norm100), n_contrib = sum(norm100 > 0), n_detected = n(), .groups = "drop")
write.csv(cs, file.path(OUT, "taxon_contributions_sig_pathways.csv"), row.names = FALSE)

## =============================================================================
## 4. NSTI
## =============================================================================
ns <- read.delim(Sys.glob(file.path(P1288, "*marker_predicted_and_nsti.tsv.gz"))[1])   # file name differs between PICRUSt2 versions
write.csv(data.frame(n_ASV_input = nrow(asv),
                     n_ASV_placed = nrow(ns), median = median(ns$metadata_NSTI), q1 = quantile(ns$metadata_NSTI, .25),
                     q3 = quantile(ns$metadata_NSTI, .75), max = max(ns$metadata_NSTI), n_NSTI_gt_2 = sum(ns$metadata_NSTI > 2)),
          file.path(OUT, "NSTI_per_ASV_summary.csv"), row.names = FALSE)
wn <- read.delim(file.path(P1288, "KO_metagenome_out", "weighted_nsti.tsv.gz")); names(wn) <- c("SampleID", "weighted_NSTI")
wn <- merge(wn, md[, c("SampleID", "Group")])
write.csv(data.frame(wn %>% group_by(Group) %>% summarise(n = n(), median = median(weighted_NSTI), min = min(weighted_NSTI), max = max(weighted_NSTI)),
                     Wilcoxon_P = wilcox.test(weighted_NSTI ~ Group, wn, exact = FALSE)$p.value),
          file.path(OUT, "NSTI_weighted_by_group.csv"), row.names = FALSE)

## =============================================================================
## 5. Richness diagnostics and pathway detection
## =============================================================================
X <- model.matrix(~ Group + Age + logRich, md)[, -1]
vif <- sapply(1:3, function(j) 1 / (1 - summary(lm(X[, j] ~ X[, -j]))$r.squared)); names(vif) <- colnames(X); print(round(vif, 2))
md$n_pathways <- rowSums(ab > 0)                                            # pathways detected per sample
ct <- cor.test(md$n_pathways, md$observed_features, method = "spearman", exact = FALSE)
cat("pathways detected per sample vs observed features: rho =", round(ct$estimate, 2), "P =", signif(ct$p.value, 2), "\n")

## detection of the six pathways by group (Fisher's exact test, BH across the six pathways)
fis <- do.call(rbind, lapply(six, function(p) { d <- ab[, p] > 0
  tb <- table(factor(d, c(TRUE, FALSE)), md$Group)
  data.frame(pathway = p, detected_single = sum(d & md$Group == "single"), detected_multiple = sum(d & md$Group == "multiple"),
             P = fisher.test(tb)$p.value) }))
fis$q_BH <- p.adjust(fis$P, "BH")
write.csv(fis, file.path(OUT, "pathway_detection_fisher.csv"), row.names = FALSE); print(fis, digits = 3)

## detection by presence of Segatella copri
sp <- rowsum(as.matrix(asv[, -(1:2)]), taxon_label(asv$Taxon)); scopri <- sp["Segatella copri", rownames(md)] > 0
sc <- do.call(rbind, lapply(six, function(p) { d <- ab[, p] > 0
  data.frame(pathway = p, detected_scopri_pos = sum(d & scopri), n_scopri_pos = sum(scopri),
             detected_scopri_neg = sum(d & !scopri), n_scopri_neg = sum(!scopri),
             P = fisher.test(table(factor(d, c(TRUE, FALSE)), factor(scopri, c(TRUE, FALSE))))$p.value) }))
write.csv(sc, file.path(OUT, "pathway_detection_by_Scopri.csv"), row.names = FALSE); print(sc, digits = 3)

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
