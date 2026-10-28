# Diversity: regression models, PERMANOVA robustness, covariate adjustment,
# three-group comparison, and exclusion of HPV52-positive or SCC cases.
# Run after 03_diversity.R.
library(vegan)
library(dplyr)
library(ecodist)
source("config.R")

out  <- out_dir("diversity", "sensitivity")
prim <- out_dir("diversity")
meta <- read_metadata()
meta$Group3 <- factor(ifelse(meta$Group == "single", "single",
                             ifelse(meta$n_major_hrHPV >= 2, "multiple-2", "multiple-1")),
                      levels = c("single", "multiple-1", "multiple-2"))

# QIIME 2 core-metrics output (40,000 reads)
cm <- function(...) file.path(data_path("core_metrics"), ...)
idx <- c(Observed = "observed_features", FaithPD = "faith_pd", Shannon = "shannon", Evenness = "evenness")
alpha <- cbind(meta, sapply(idx, function(n)
  read.delim(cm(paste0(n, "_vector"), "alpha-diversity.tsv"), row.names = 1)[meta$SampleID, 1]))
read_dm <- function(n) as.dist(as.matrix(read.delim(cm(paste0(n, "_distance_matrix"), "distance-matrix.tsv"),
                                                    row.names = 1, check.names = FALSE)))
dms <- list("Unweighted UniFrac" = read_dm("unweighted_unifrac"), "Jaccard" = read_dm("jaccard"))
meta_of <- function(d) meta[labels(d), ]

wilcox_tab <- function(d) data.frame(Index = names(idx), P = sapply(names(idx), function(v)
  suppressWarnings(wilcox.test(d[[v]] ~ d$Group))$p.value))

# Alpha diversity ------------------------------------------------------------
w <- wilcox_tab(alpha); w$q <- p.adjust(w$P, "BH")
write.csv(w, file.path(out, "alpha_wilcoxon.csv"), row.names = FALSE)
write.csv(data.frame(Index = names(idx), P = sapply(names(idx), function(v) kruskal.test(alpha[[v]] ~ alpha$Group3)$p.value)),
          file.path(out, "alpha_kruskal_3groups.csv"), row.names = FALSE)

# regression: quasi-Poisson for observed features (IRR), linear for the others
alpha_reg <- function(d, rhs, label) {
  do.call(rbind, lapply(names(idx), function(v) {
    f <- reformulate(rhs, v)
    m <- if (v == "Observed") glm(f, quasipoisson, d) else lm(f, d)
    co <- summary(m)$coefficients[-1, , drop = FALSE]
    ci <- suppressMessages(confint(m))[rownames(co), , drop = FALSE]
    est <- cbind(co[, 1], ci)
    if (v == "Observed") est <- exp(est)
    data.frame(Index = v, Model = label, Term = rownames(co), Estimate = est[, 1],
               CI_low = est[, 2], CI_high = est[, 3], P = co[, 4], row.names = NULL)
  }))
}
neg52 <- meta$HPV52 == "N"
reg <- rbind(alpha_reg(alpha, "Group", "Unadjusted"),
             alpha_reg(alpha, c("Group", "Age"), "Age"),
             alpha_reg(alpha, c("Group", "HPV52"), "HPV52"),
             alpha_reg(alpha, c("Group", "Age", "HPV52"), "Age + HPV52"),
             alpha_reg(alpha[neg52, ], "Group", "HPV52-negative"),
             alpha_reg(alpha[neg52, ], c("Group", "Age"), "HPV52-negative, age"))
write.csv(reg, file.path(out, "alpha_regression.csv"), row.names = FALSE)

# Beta diversity: PERMANOVA, 95% subsampling range of R2, PERMDISP, ANOSIM, Mantel, MRM
# (R2 range: 80% of each group sampled without replacement, 999 times)
subsample_r2 <- function(d, mdf, rhs) {
  m <- as.matrix(d)
  by_grp <- split(seq_len(nrow(mdf)), mdf$Group)
  set.seed(seed)
  r2 <- vapply(seq_len(999), function(b) {
    i <- unlist(lapply(by_grp, function(x) sample(x, round(0.8 * length(x)))))
    f <- as.formula(paste("as.dist(m[i, i]) ~", paste(rhs, collapse = " + ")))
    a <- tryCatch(adonis2(f, data = mdf[i, , drop = FALSE], permutations = 99, by = "terms"),
                  error = function(e) NULL)
    if (is.null(a)) NA_real_ else a["Group", "R2"]
  }, numeric(1))
  quantile(r2, c(0.025, 0.975), na.rm = TRUE)
}
robust <- do.call(rbind, lapply(names(dms), function(k) {
  d <- dms[[k]]; mdf <- meta_of(d)
  design <- dist(as.numeric(mdf$Group) - 1, "manhattan")              # 0 = same group, 1 = different
  set.seed(seed); ad <- adonis2(d ~ Group, mdf, permutations = n_perm)
  set.seed(seed); pd <- permutest(betadisper(d, mdf$Group, type = "median"), permutations = n_perm)
  set.seed(seed); an <- anosim(d, mdf$Group, permutations = n_perm)
  set.seed(seed); mt <- vegan::mantel(d, design, permutations = n_perm)
  set.seed(seed); mr <- MRM(d ~ design, nperm = n_perm)
  r <- subsample_r2(d, mdf, "Group")
  data.frame(Distance = k, R2 = ad$R2[1], R2_range_low = r[1], R2_range_high = r[2], P = ad$`Pr(>F)`[1],
             PERMDISP_P = pd$tab$`Pr(>F)`[1], ANOSIM_R = an$statistic, ANOSIM_P = an$signif,
             Mantel_r = mt$statistic, Mantel_P = mt$signif, MRM_R2 = mr$r.squared[1], MRM_P = mr$r.squared[2])
}))
write.csv(robust, file.path(out, "beta_robustness.csv"), row.names = FALSE)

# PERMANOVA adjusted for age or HPV52 (distance ~ covariate + Group, sequential terms)
adj <- do.call(rbind, lapply(names(dms), function(k) do.call(rbind, lapply(c("Age", "HPV52"), function(cv) {
  d <- dms[[k]]; mdf <- meta_of(d)
  set.seed(seed); a <- adonis2(reformulate(c(cv, "Group"), "d"), mdf, permutations = n_perm, by = "terms")
  r <- subsample_r2(d, mdf, c(cv, "Group"))
  data.frame(Distance = k, Covariate = cv, Group_R2 = a["Group", "R2"], R2_range_low = r[1], R2_range_high = r[2],
             Group_P = a["Group", "Pr(>F)"], Covariate_P = a[cv, "Pr(>F)"])
}))))
write.csv(adj, file.path(out, "beta_covariate_adjusted.csv"), row.names = FALSE)

# Three groups (unweighted UniFrac) and pairwise distances
d <- dms[["Unweighted UniFrac"]]; mdf <- meta_of(d)
set.seed(seed); print(adonis2(d ~ Group3, mdf, permutations = n_perm))
m <- as.matrix(d); ut <- which(upper.tri(m), arr.ind = TRUE)
pairs <- data.frame(dist = m[ut], g1 = mdf$Group[ut[, 1]], g2 = mdf$Group[ut[, 2]])
pairs$cmp <- ifelse(pairs$g1 == pairs$g2, paste(pairs$g1, pairs$g1, sep = "-"), "single-multiple")
write.csv(aggregate(dist ~ cmp, pairs, function(x) c(median = median(x), quantile(x, c(.25, .75)))),
          file.path(out, "pairwise_unweighted_unifrac.csv"), row.names = FALSE)

# Joint adjustment and HPV52 exclusion for the five distances (03_diversity.R output)
dnames <- c("jaccard", "unweighted_unifrac", "braycurtis", "weighted_unifrac", "aitchison")
joint <- do.call(rbind, lapply(dnames, function(k) {
  m  <- as.matrix(read.csv(file.path(prim, paste0("primary_dm_", k, ".csv")), row.names = 1, check.names = FALSE))
  m  <- m[rownames(meta), rownames(meta)]
  mn <- meta[neg52, ]; dn <- as.dist(m[neg52, neg52])
  set.seed(seed); a  <- adonis2(as.dist(m) ~ Age + HPV52 + Group, meta, permutations = n_perm, by = "terms")
  set.seed(seed); a0 <- adonis2(dn ~ Group, mn, permutations = n_perm)
  set.seed(seed); a1 <- adonis2(dn ~ Age + Group, mn, permutations = n_perm, by = "terms")
  set.seed(seed); pd <- permutest(betadisper(dn, mn$Group, type = "median"), permutations = n_perm)
  data.frame(Distance = k,
             Analysis = c("Age + HPV52", "HPV52-negative", "HPV52-negative, age", "HPV52-negative, PERMDISP"),
             R2 = c(a["Group", "R2"], a0$R2[1], a1["Group", "R2"], NA),
             P  = c(a["Group", "Pr(>F)"], a0$`Pr(>F)`[1], a1["Group", "Pr(>F)"], pd$tab$`Pr(>F)`[1]))
}))
write.csv(joint, file.path(out, "beta_joint_and_HPV52_excluded.csv"), row.names = FALSE)
write.csv(wilcox_tab(alpha[neg52, ]), file.path(out, "alpha_wilcoxon_HPV52_excluded.csv"), row.names = FALSE)

# SCC cases excluded (n = 34; PERMANOVA and PERMDISP with 9,999 permutations)
keep <- meta$SampleID[meta$Dx != "SCC"]
write.csv(wilcox_tab(alpha[keep, ]), file.path(out, "alpha_wilcoxon_SCC_excluded.csv"), row.names = FALSE)
k <- meta[keep, ]
write.csv(data.frame(Test = c("CST", "CST III subtypes", "CST IV subtypes"),
                     P = c(fisher.test(table(k$Group, k$CST))$p.value,
                           fisher.test(with(k[k$CST == "III", ], table(Group, CST_subtype)))$p.value,
                           fisher.test(with(k[k$CST == "IV", ],  table(Group, CST_subtype)))$p.value)),
          file.path(out, "CST_SCC_excluded.csv"), row.names = FALSE)
scc <- do.call(rbind, lapply(names(dms), function(nm) {
  d <- as.dist(as.matrix(dms[[nm]])[keep, keep]); g <- meta[labels(d), "Group"]
  set.seed(seed); a  <- adonis2(d ~ g, permutations = 9999)
  set.seed(seed); pd <- permutest(betadisper(d, g, type = "median"), permutations = 9999)
  data.frame(Distance = nm, R2 = a$R2[1], P = a$`Pr(>F)`[1], PERMDISP_P = pd$tab$`Pr(>F)`[1])
}))
write.csv(scc, file.path(out, "beta_SCC_excluded.csv"), row.names = FALSE)
