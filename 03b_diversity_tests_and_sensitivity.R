## =============================================================================
## 03b_diversity_tests_and_sensitivity.R
## Diversity tests, covariate adjustment, and sensitivity analyses
## (Figure 2A statistics; Supplementary Tables S8, S9, S11, S13, S14; Figure S2)
## Run after 03_diversity.R.
##
##  1. Alpha diversity: Wilcoxon (BH across the four indices), Kruskal-Wallis
##     across single / multiple-1 / multiple-2, and regression models
##     (observed features: quasi-Poisson, IRR; other indices: linear, beta)
##     unadjusted, adjusted for age, and adjusted for HPV52
##  2. Beta diversity (unweighted UniFrac, Jaccard): PERMANOVA with a 95%
##     subsampling range of R2, PERMDISP, ANOSIM, Mantel test, MRM
##  3. PERMANOVA adjusted for age or HPV52 (distance ~ covariate + Group)
##  4. Three-group PERMANOVA and pairwise distances (descriptive)
##  5. Joint adjustment for age and HPV52, and exclusion of HPV52-positive
##     participants (alpha; beta for five distances)
##  6. Exclusion of participants with squamous cell carcinoma (alpha, CST, beta)
##
## Inputs : data/core_metrics/ (QIIME 2 core-metrics, 40,000 reads)
##          results/diversity/primary_and_rarefaction/ (03_diversity.R)
##          data/sample_metadata.csv  (additionally: n_major_hrHPV, CST, CST_subtype)
## =============================================================================
suppressPackageStartupMessages({ library(vegan); library(dplyr); library(ecodist) })
source("config.R")
OUT <- file.path(RESULTS_DIR, "diversity", "tests_and_sensitivity"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
P03 <- file.path(RESULTS_DIR, "diversity", "primary_and_rarefaction")

meta <- read_metadata()
## multiple-1: one major hrHPV genotype; multiple-2: two or more major genotypes
meta$Group3 <- factor(ifelse(meta$Group == "single", "single", ifelse(meta$n_major_hrHPV >= 2, "multiple-2", "multiple-1")),
                      levels = c("single", "multiple-1", "multiple-2"))
print(table(meta$Group)); print(table(meta$Group3)); print(table(meta$HPV52, meta$Group))

CM <- data_path("core_metrics")
read_alpha <- function(n) { v <- read.delim(file.path(CM, paste0(n, "_vector"), "alpha-diversity.tsv"), row.names = 1); setNames(v[, 1], rownames(v)) }
df <- meta
df$Observed <- read_alpha("observed_features")[df$SampleID]; df$FaithPD <- read_alpha("faith_pd")[df$SampleID]
df$Shannon  <- read_alpha("shannon")[df$SampleID];           df$Evenness <- read_alpha("evenness")[df$SampleID]
stopifnot(!any(is.na(df[, c("Observed", "FaithPD", "Shannon", "Evenness")])))
idx <- c(Observed = "Observed features", FaithPD = "Faith's PD", Shannon = "Shannon", Evenness = "Pielou's evenness")

read_dm <- function(n) { m <- as.matrix(read.delim(file.path(CM, paste0(n, "_distance_matrix"), "distance-matrix.tsv"), row.names = 1, check.names = FALSE))
  stopifnot(all(rownames(m) == colnames(m))); as.dist(m) }
dms <- list("Unweighted UniFrac" = read_dm("unweighted_unifrac"), "Jaccard" = read_dm("jaccard"))
align_meta <- function(d) { m <- meta[match(labels(d), meta$SampleID), ]; stopifnot(all(labels(d) == m$SampleID)); m }

## =============================================================================
## 1. Alpha diversity
## =============================================================================
wilcox_tab <- bind_rows(lapply(names(idx), function(v) { w <- wilcox.test(df[[v]] ~ df$Group)
  data.frame(Index = idx[[v]], Median_single = median(df[[v]][df$Group == "single"]),
             Median_multiple = median(df[[v]][df$Group == "multiple"]), W = unname(w$statistic), P = w$p.value) }))
wilcox_tab$q_BH <- p.adjust(wilcox_tab$P, "BH")                                         # BH across the four indices
write.csv(wilcox_tab, file.path(OUT, "alpha_wilcoxon.csv"), row.names = FALSE); print(wilcox_tab, digits = 4)

kw_tab <- bind_rows(lapply(names(idx), function(v) { k <- kruskal.test(df[[v]] ~ df$Group3)
  data.frame(Index = idx[[v]], KW_chi2 = unname(k$statistic), df = unname(k$parameter), P = k$p.value) }))
write.csv(kw_tab, file.path(OUT, "alpha_3group_kruskal.csv"), row.names = FALSE)

get_term <- function(model, term, exp_it = FALSE) {
  co <- summary(model)$coefficients; ci <- suppressMessages(confint(model))
  e <- c(co[term, 1], ci[term, ]); if (exp_it) e <- exp(e)
  data.frame(Estimate = e[1], CI_low = e[2], CI_high = e[3], P = co[term, 4]) }
fit <- function(v, rhs, d) if (v == "Observed") glm(as.formula(paste(v, rhs)), data = d, family = quasipoisson()) else lm(as.formula(paste(v, rhs)), data = d)
terms <- c(Groupmultiple = "Group (multiple vs single)", Age = "Age (per 1-year increase)", HPV52Y = "HPV52 (positive vs negative)")
alpha_models <- function(model_set, data_list) {
  out <- list()
  for (v in names(idx)) { ex <- v == "Observed"
    for (mm in model_set) { m <- fit(v, mm$rhs, data_list[[mm$data]])
      for (t in intersect(names(terms), rownames(summary(m)$coefficients)))
        out[[length(out) + 1]] <- cbind(Index = idx[[v]], Model = mm$label, Predictor = terms[[t]],
                                        Scale = if (ex) "IRR" else "beta", get_term(m, t, ex)) } }
  r <- do.call(rbind, out); rownames(r) <- NULL; r }
alpha_adj <- alpha_models(list(list(rhs = "~ Group", label = "Unadjusted", data = "all"),
                               list(rhs = "~ Group + Age", label = "Adjusted for age", data = "all"),
                               list(rhs = "~ Group + HPV52", label = "Adjusted for HPV52", data = "all")), list(all = df))
write.csv(alpha_adj, file.path(OUT, "alpha_regression_single_covariate.csv"), row.names = FALSE)

## =============================================================================
## 2. Beta diversity: PERMANOVA, 95% subsampling range of R2, PERMDISP,
##    ANOSIM, Mantel test, MRM
##    R2 range: stratified subsampling without replacement (80% of each group,
##    999 replicates), 2.5th-97.5th percentiles. Subsampling is used instead of
##    bootstrapping because duplicated samples have a distance of 0.
## =============================================================================
N_SUB <- 999; SUB_FRAC <- 0.8
subsample_r2 <- function(d, mdf, rhs, term_name) {
  m <- as.matrix(d); idx_by <- split(seq_len(nrow(mdf)), mdf$Group)
  set.seed(SEED)
  vals <- vapply(seq_len(N_SUB), function(b) {
    bi <- unlist(lapply(idx_by, function(ix) sample(ix, round(SUB_FRAC * length(ix)), replace = FALSE)))
    a <- tryCatch(adonis2(as.formula(paste("as.dist(m[bi, bi]) ~", rhs)), data = mdf[bi, , drop = FALSE], permutations = 99, by = "terms"),
                  error = function(e) NULL)
    if (is.null(a)) NA_real_ else a$R2[rownames(a) == term_name] }, numeric(1))
  quantile(vals, c(0.025, 0.975), na.rm = TRUE) }
beta_rob <- bind_rows(lapply(names(dms), function(nm) {
  d <- dms[[nm]]; mdf <- align_meta(d)
  set.seed(SEED); ad <- adonis2(d ~ Group, data = mdf, permutations = N_PERM, by = "terms")
  rg <- subsample_r2(d, mdf, "Group", "Group")
  disp <- betadisper(d, mdf$Group, type = "median"); set.seed(SEED); pt <- permutest(disp, permutations = N_PERM)
  set.seed(SEED); an <- anosim(d, mdf$Group, permutations = N_PERM)
  design <- outer(mdf$Group, mdf$Group, function(a, b) ifelse(a == b, 0, 1)); diag(design) <- 0
  set.seed(SEED); mt <- vegan::mantel(d, as.dist(design), method = "pearson", permutations = N_PERM)
  set.seed(SEED); mr <- ecodist::MRM(d ~ as.dist(design), nperm = N_PERM)
  data.frame(Metric = nm, PERMANOVA_R2 = ad$R2[1], R2_subsampling_2.5 = rg[1], R2_subsampling_97.5 = rg[2],
             PERMANOVA_F = ad$F[1], PERMANOVA_P = ad$`Pr(>F)`[1], PERMDISP_F = pt$tab$F[1], PERMDISP_P = pt$tab$`Pr(>F)`[1],
             ANOSIM_R = an$statistic, ANOSIM_P = an$signif, Mantel_r = mt$statistic, Mantel_P = mt$signif,
             MRM_R2 = mr$r.squared[1], MRM_P = mr$r.squared[2]) }))
write.csv(beta_rob, file.path(OUT, "beta_permanova_robustness.csv"), row.names = FALSE); print(t(beta_rob))

## =============================================================================
## 3. PERMANOVA adjusted for age or HPV52 (sequential terms: distance ~ covariate + Group)
## =============================================================================
beta_adj <- bind_rows(lapply(names(dms), function(nm) { d <- dms[[nm]]; mdf <- align_meta(d)
  bind_rows(lapply(c("Age", "HPV52"), function(cv) {
    set.seed(SEED); ad <- adonis2(as.formula(paste("d ~", cv, "+ Group")), data = mdf, permutations = N_PERM, by = "terms")
    rg <- subsample_r2(d, mdf, paste(cv, "+ Group"), "Group")
    data.frame(Metric = nm, Covariate = cv, Group_R2 = ad$R2[2], Group_R2_subsampling_2.5 = rg[1], Group_R2_subsampling_97.5 = rg[2],
               Group_F = ad$F[2], Group_P = ad$`Pr(>F)`[2], Cov_R2 = ad$R2[1], Cov_F = ad$F[1], Cov_P = ad$`Pr(>F)`[1]) })) }))
write.csv(beta_adj, file.path(OUT, "beta_permanova_single_covariate.csv"), row.names = FALSE); print(beta_adj, digits = 4)

## =============================================================================
## 4. Three-group PERMANOVA (unweighted UniFrac) and pairwise distances
## =============================================================================
d <- dms[["Unweighted UniFrac"]]; mdf <- align_meta(d)
set.seed(SEED); ad3 <- adonis2(d ~ Group3, data = mdf, permutations = N_PERM, by = "terms")
write.csv(data.frame(Metric = "Unweighted UniFrac", R2 = ad3$R2[1], F = ad3$F[1], P = ad3$`Pr(>F)`[1]),
          file.path(OUT, "beta_3group_permanova.csv"), row.names = FALSE)
m <- as.matrix(d); pw <- which(upper.tri(m), arr.ind = TRUE)
pair <- data.frame(Distance = m[pw], g1 = as.character(mdf$Group[pw[, 1]]), g2 = as.character(mdf$Group[pw[, 2]])) %>%
  mutate(Comparison = ifelse(g1 == g2, paste0(g1, "-", g1), "multiple-single"))
write.csv(pair %>% group_by(Comparison) %>% summarise(n_pairs = n(), Median = median(Distance), Q1 = quantile(Distance, .25),
                                                      Q3 = quantile(Distance, .75), .groups = "drop"),
          file.path(OUT, "pairwise_unweighted_unifrac_summary.csv"), row.names = FALSE)

## =============================================================================
## 5. Joint adjustment for age and HPV52; exclusion of HPV52-positive participants
##    (alpha and distance matrices from 03_diversity.R)
## =============================================================================
keep52 <- rownames(meta)[meta$HPV52 == "N"]; print(table(meta[keep52, "Group"]))
al <- read.csv(file.path(P03, "primary_alpha.csv"), row.names = 1); da <- cbind(meta, al[rownames(meta), ])
alpha_joint <- alpha_models(list(list(rhs = "~ Group + Age + HPV52", label = "Adjusted for age and HPV52 (n = 37)", data = "all"),
                                 list(rhs = "~ Group", label = "HPV52-negative only, unadjusted (n = 32)", data = "neg"),
                                 list(rhs = "~ Group + Age", label = "HPV52-negative only, adjusted for age (n = 32)", data = "neg")),
                            list(all = da, neg = da[keep52, ]))
w52 <- do.call(rbind, lapply(names(idx), function(v) { d <- da[keep52, ]; w <- suppressWarnings(wilcox.test(d[[v]] ~ d$Group))
  data.frame(Index = idx[[v]], Median_single = median(d[[v]][d$Group == "single"]), Median_multiple = median(d[[v]][d$Group == "multiple"]),
             W = unname(w$statistic), P = w$p.value) }))
write.csv(alpha_joint, file.path(OUT, "alpha_joint_and_HPV52_excluded_regression.csv"), row.names = FALSE)
write.csv(w52, file.path(OUT, "alpha_wilcoxon_HPV52_excluded.csv"), row.names = FALSE); print(w52, digits = 3)

dists <- c(jaccard = "Jaccard", unweighted_unifrac = "Unweighted UniFrac", braycurtis = "Bray-Curtis",
           weighted_unifrac = "Weighted UniFrac", aitchison = "Aitchison")
Bt <- list()
for (k in names(dists)) {
  mm <- as.matrix(read.csv(file.path(P03, sprintf("primary_dm_%s.csv", k)), row.names = 1, check.names = FALSE))[rownames(meta), rownames(meta)]
  dd <- as.dist(mm)
  set.seed(SEED); a <- adonis2(dd ~ Age + HPV52 + Group, data = meta, permutations = N_PERM, by = "terms")
  for (t in c("Age", "HPV52", "Group"))
    Bt[[length(Bt) + 1]] <- data.frame(Distance = dists[[k]], Analysis = "Adjusted for age and HPV52 (n = 37)", Term = t,
                                       R2 = a[t, "R2"], F = a[t, "F"], P = a[t, "Pr(>F)"])
  mn <- meta[keep52, ]; dn <- as.dist(mm[keep52, keep52])
  set.seed(SEED); a0 <- adonis2(dn ~ Group, data = mn, permutations = N_PERM)
  set.seed(SEED); a1 <- adonis2(dn ~ Age + Group, data = mn, permutations = N_PERM, by = "terms")
  bd <- betadisper(dn, mn$Group, type = "median"); set.seed(SEED); pt <- permutest(bd, permutations = N_PERM)
  Bt[[length(Bt) + 1]] <- data.frame(Distance = dists[[k]], Analysis = "HPV52-negative only (n = 32)", Term = "Group",
                                     R2 = a0$R2[1], F = a0$F[1], P = a0$`Pr(>F)`[1])
  Bt[[length(Bt) + 1]] <- data.frame(Distance = dists[[k]], Analysis = "HPV52-negative only, adjusted for age (n = 32)", Term = "Group",
                                     R2 = a1["Group", "R2"], F = a1["Group", "F"], P = a1["Group", "Pr(>F)"])
  Bt[[length(Bt) + 1]] <- data.frame(Distance = dists[[k]], Analysis = "HPV52-negative only (n = 32)", Term = "PERMDISP",
                                     R2 = NA, F = pt$tab$F[1], P = pt$tab$`Pr(>F)`[1])
}
Bt <- do.call(rbind, Bt); write.csv(Bt, file.path(OUT, "beta_joint_and_HPV52_excluded.csv"), row.names = FALSE); print(Bt, digits = 3)

## =============================================================================
## 6. Exclusion of participants with squamous cell carcinoma (n = 34)
## =============================================================================
keep <- meta[meta$Dx != "SCC", ]
cat("SCC excluded:", paste(setdiff(meta$SampleID, keep$SampleID), collapse = ", "), "| n =", nrow(keep), "\n")
## alpha diversity (Wilcoxon; same test as the primary analysis)
dk <- df[keep$SampleID, ]
w_scc <- do.call(rbind, lapply(names(idx), function(v) { w <- suppressWarnings(wilcox.test(dk[[v]] ~ dk$Group))
  data.frame(Index = idx[[v]], Median_single = median(dk[[v]][dk$Group == "single"]), Median_multiple = median(dk[[v]][dk$Group == "multiple"]),
             W = unname(w$statistic), P = w$p.value) }))
write.csv(w_scc, file.path(OUT, "alpha_wilcoxon_SCC_excluded.csv"), row.names = FALSE); print(w_scc, digits = 3)
## CST distributions (Fisher's exact test)
cst <- table(keep$Group, keep$CST)
c3  <- table(keep$Group[keep$CST == "III"], keep$CST_subtype[keep$CST == "III"])
c4  <- table(keep$Group[keep$CST == "IV"],  keep$CST_subtype[keep$CST == "IV"])
write.csv(data.frame(Test = c("Overall CST (I-V)", "CST III subtypes", "CST IV subtypes"),
                     P = c(fisher.test(cst)$p.value, fisher.test(c3)$p.value, fisher.test(c4)$p.value)),
          file.path(OUT, "CST_fisher_SCC_excluded.csv"), row.names = FALSE)
## beta diversity (PERMANOVA and PERMDISP, 9,999 permutations)
beta_scc <- do.call(rbind, lapply(names(dms), function(nm) {
  m <- as.matrix(dms[[nm]])[keep$SampleID, keep$SampleID]; d <- as.dist(m); g <- keep$Group[match(labels(d), keep$SampleID)]
  set.seed(SEED); ad <- adonis2(d ~ g, permutations = 9999, by = "terms")
  disp <- betadisper(d, g, type = "median"); set.seed(SEED); pt <- permutest(disp, permutations = 9999)
  data.frame(Metric = nm, R2 = ad$R2[1], F = ad$F[1], P = ad$`Pr(>F)`[1], PERMDISP_F = pt$tab$F[1], PERMDISP_P = pt$tab$`Pr(>F)`[1]) }))
write.csv(beta_scc, file.path(OUT, "beta_SCC_excluded.csv"), row.names = FALSE); print(beta_scc, digits = 4)

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
