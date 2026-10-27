## =============================================================================
## 03_diversity.R
## Alpha and beta diversity (Supplementary Tables S6, S7, S10-S12; Figure S1 data)
##
##  1. Read depth by group and rarefaction-curve data
##  2. Primary analysis on the QIIME 2 table rarefied to 40,000 reads
##     (alpha: Wilcoxon rank-sum, BH across the four indices;
##      beta : PERMANOVA + PERMDISP for five distances, BH across distances)
##  3. 100 independent rarefactions and rare-ASV filters
##
## Faith's PD and (un)weighted UniFrac are computed from the rooted tree
## (branch x ASV incidence matrix) following the QIIME 2 definitions (Faith's PD
## includes all branches to the root; weighted UniFrac is un-normalized). The
## values for the primary rarefied table are compared with the QIIME 2
## core-metrics output (section 2).
##
## Inputs : data/feature_table_1288.tsv     raw counts, 1,288 ASVs x 37 samples
##          data/rarefied_table_40000.tsv   QIIME 2 core-metrics rarefied table
##          data/tree.nwk                   rooted phylogenetic tree
##          data/core_metrics/              exported QIIME 2 core-metrics output
##          data/sample_metadata.csv
## =============================================================================
suppressPackageStartupMessages({ library(vegan); library(ape); library(Matrix); library(parallel) })
source("config.R")
O1 <- file.path(RESULTS_DIR, "diversity", "read_depth")
O2 <- file.path(RESULTS_DIR, "diversity", "primary_and_rarefaction")
O3 <- file.path(RESULTS_DIR, "diversity", "rare_ASV_filters")
for (d in c(O1, O2, O3)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------- data
md    <- read_metadata()
rar_q <- as.matrix(read_biom_tsv(data_path("rarefied_table_40000.tsv")))
S     <- colnames(rar_q)          # sample order of the QIIME 2 artifacts (permutation P values depend on it)
full  <- as.matrix(read_biom_tsv(data_path("feature_table_1288.tsv")))
full  <- t(full[, S])                                                     # samples x ASVs
grp   <- md[S, "Group"]; names(grp) <- S
rar_q <- t(rar_q[colnames(full), S]); stopifnot(all(rowSums(rar_q) == DEPTH))
tree  <- read.tree(data_path("tree.nwk")); stopifnot(setequal(tree$tip.label, colnames(full)))
cat("samples", nrow(full), "| ASVs", ncol(full), "| reads", sum(full), "\n")

## ---------------------------------------------------------------- phylogenetic metrics
ntip <- length(tree$tip.label); par_edge <- integer(max(tree$edge)); par_edge[tree$edge[, 2]] <- seq_len(nrow(tree$edge))
root <- ntip + 1; ii <- integer(0); jj <- integer(0)
for (t in seq_len(ntip)) { node <- t; while (node != root) { e <- par_edge[node]; ii <- c(ii, e); jj <- c(jj, t); node <- tree$edge[e, 1] } }
Mtip <- sparseMatrix(i = ii, j = jj, x = 1, dims = c(nrow(tree$edge), ntip), dimnames = list(NULL, tree$tip.label))
blen <- tree$edge.length
branch <- function(X) as.matrix(Mtip[, colnames(X)] %*% t(X))            # edges x samples
faith  <- function(X) { B <- branch(X > 0); colSums((B > 0) * blen) }
unifrac_uw <- function(X) { B <- branch(X > 0) > 0; n <- ncol(B); D <- matrix(0, n, n, dimnames = list(rownames(X), rownames(X)))
  for (a in 1:(n - 1)) for (b in (a + 1):n) { u <- B[, a] | B[, b]; D[a, b] <- D[b, a] <- sum(blen[xor(B[, a], B[, b])]) / sum(blen[u]) }
  as.dist(D) }
unifrac_w <- function(X) { P <- branch(X / rowSums(X)); n <- ncol(P); D <- matrix(0, n, n, dimnames = list(rownames(X), rownames(X)))
  for (a in 1:(n - 1)) for (b in (a + 1):n) D[a, b] <- D[b, a] <- sum(blen * abs(P[, a] - P[, b]))
  as.dist(D) }
alpha <- function(X) { sh <- diversity(X, index = "shannon", base = 2); ob <- rowSums(X > 0)
  data.frame(SampleID = rownames(X), Observed = ob, FaithPD = faith(X), Shannon = sh, Evenness = sh / log2(ob)) }
aitchison <- function(X) { L <- log(X + 1); dist(L - rowMeans(L)) }       # CLR with pseudocount 1, Euclidean
dms <- function(X, rel_bc = FALSE) list(
  jaccard = vegdist(X > 0, "jaccard", binary = TRUE),
  braycurtis = vegdist(if (rel_bc) X / rowSums(X) else X, "bray"),
  unweighted_unifrac = unifrac_uw(X), weighted_unifrac = unifrac_w(X), aitchison = aitchison(X))

## ---------------------------------------------------------------- tests
DL  <- c(jaccard = "Jaccard", unweighted_unifrac = "Unweighted UniFrac", braycurtis = "Bray-Curtis",
         weighted_unifrac = "Weighted UniFrac", aitchison = "Aitchison")
IDX <- c(Observed = "Observed features", FaithPD = "Faith's PD", Shannon = "Shannon", Evenness = "Pielou's evenness")
test_dm <- function(d) { g <- grp[labels(d)]
  set.seed(SEED); ad <- adonis2(d ~ g, permutations = N_PERM)
  bd <- betadisper(d, g, type = "median"); set.seed(SEED); pt <- permutest(bd, permutations = N_PERM)
  c(R2 = ad$R2[1], F = ad$F[1], P = ad$`Pr(>F)`[1], disp_F = pt$tab$F[1], disp_P = pt$tab$`Pr(>F)`[1],
    disp_single = mean(bd$distances[g == "single"]), disp_multiple = mean(bd$distances[g == "multiple"])) }
wil <- function(x) suppressWarnings(wilcox.test(x ~ grp[names(x)]))
alpha_tab <- function(a) { r <- do.call(rbind, lapply(names(IDX), function(v) { x <- setNames(a[[v]], a$SampleID); ww <- wil(x)
    qs <- function(z) sprintf("%.2f (%.2f-%.2f)", median(z), quantile(z, .25), quantile(z, .75))
    data.frame(Index = IDX[[v]], Single = qs(x[grp[names(x)] == "single"]), Multiple = qs(x[grp[names(x)] == "multiple"]),
               W = unname(ww$statistic), P = ww$p.value) }))
  r$q_BH = p.adjust(r$P, "BH"); r }                                       # BH across the four indices
beta_tab <- function(D) { r <- do.call(rbind, lapply(names(DL), function(k) data.frame(Distance = DL[[k]], t(test_dm(D[[k]])))))
  r$q_BH = p.adjust(r$P, "BH"); r }                                       # BH across the five distances

## =============================================================================
## 1. Read depth and rarefaction-curve data (Figure S1)
## =============================================================================
depth <- rowSums(full)
d1 <- data.frame(SampleID = S, group = grp, reads = depth, observed_unrarefied = rowSums(full > 0),
                 faith_unrarefied = faith(full), shannon_unrarefied = diversity(full, "shannon", base = 2))
write.csv(d1, file.path(O1, "read_depth_per_sample.csv"), row.names = FALSE)
sm <- function(x) c(n = length(x), total = sum(x), mean = mean(x), sd = sd(x), median = median(x),
                    q1 = unname(quantile(x, .25)), q3 = unname(quantile(x, .75)), min = min(x), max = max(x))
summ <- rbind(data.frame(group = "single", t(sm(depth[grp == "single"]))), data.frame(group = "multiple", t(sm(depth[grp == "multiple"]))),
              data.frame(group = "all", t(sm(depth))))
summ$wilcoxon_P <- c(wilcox.test(reads ~ group, data = d1)$p.value, NA, NA)
write.csv(summ, file.path(O1, "read_depth_by_group.csv"), row.names = FALSE); print(summ, digits = 6)
cr <- do.call(rbind, lapply(c("observed_unrarefied", "faith_unrarefied", "shannon_unrarefied"), function(v)
  do.call(rbind, lapply(c("all", "single", "multiple"), function(sc) {
    dd <- if (sc == "all") d1 else d1[d1$group == sc, ]; ct <- suppressWarnings(cor.test(dd$reads, dd[[v]], method = "spearman"))
    data.frame(metric = v, scope = sc, rho = unname(ct$estimate), P = ct$p.value) }))))
write.csv(cr, file.path(O1, "read_depth_vs_richness_correlation.csv"), row.names = FALSE)

## rarefaction curves: 10 subsamples per sample and depth
depths <- sort(unique(c(500, 1000, 2500, 5000, seq(7500, 92500, 2500), DEPTH, min(depth))))
set.seed(SEED)
rc <- do.call(rbind, lapply(S, function(s) { x <- full[s, , drop = FALSE]
  do.call(rbind, lapply(depths[depths <= sum(x)], function(dp) {
    R <- rrarefy(x[rep(1, 10), , drop = FALSE], dp); rownames(R) <- paste0(s, "_", 1:10)
    data.frame(SampleID = s, depth = dp, iter = 1:10, observed = rowSums(R > 0), faith_pd = faith(R),
               shannon = diversity(R, "shannon", base = 2), goods_coverage = 1 - rowSums(R == 1) / dp) })) }))
rm_ <- aggregate(cbind(observed, faith_pd, shannon, goods_coverage) ~ SampleID + depth, rc, mean); rm_$group <- grp[rm_$SampleID]
write.csv(rm_, file.path(O1, "rarefaction_curve_mean_per_sample.csv"), row.names = FALSE)
## depth diagnostics: Good's coverage, richness relative to the full table, and new ASVs per 1,000 reads
diag <- do.call(rbind, lapply(c(10000, 20000, 30000, DEPTH, min(depth)), function(dp) { prev <- max(depths[depths < dp])
  a <- rm_[rm_$depth == dp, ]; b <- rm_[rm_$depth == prev, ]; b <- b[match(a$SampleID, b$SampleID), ]
  frac <- a$observed / d1$observed_unrarefied[match(a$SampleID, d1$SampleID)]; slope <- (a$observed - b$observed) / ((dp - prev) / 1000)
  data.frame(depth = dp, n_retained = nrow(a), goods_coverage_median = median(a$goods_coverage), goods_coverage_min = min(a$goods_coverage),
             pct_richness_vs_full_median = 100 * median(frac), pct_richness_vs_full_min = 100 * min(frac),
             new_ASV_per_1000reads_median = median(slope), new_ASV_per_1000reads_max = max(slope)) }))
write.csv(diag, file.path(O1, "rarefaction_depth_diagnostics.csv"), row.names = FALSE); print(diag, digits = 4)

## =============================================================================
## 2. Primary analysis: QIIME 2 table rarefied to 40,000 reads
## =============================================================================
CM <- data_path("core_metrics")
qa <- function(n) { v <- read.delim(file.path(CM, paste0(n, "_vector"), "alpha-diversity.tsv"), row.names = 1); setNames(v[, 1], rownames(v)) }
qd <- function(n) as.matrix(read.delim(file.path(CM, paste0(n, "_distance_matrix"), "distance-matrix.tsv"), row.names = 1, check.names = FALSE))[S, S]
pa <- alpha(rar_q); pdm <- dms(rar_q)
## comparison with the QIIME 2 core-metrics output (differences should be ~0)
q_alpha <- data.frame(Observed = qa("observed_features")[S], FaithPD = qa("faith_pd")[S], Shannon = qa("shannon")[S], Evenness = qa("evenness")[S])
cat("alpha, R vs QIIME 2, max|diff|:", sapply(names(IDX), function(v) signif(max(abs(pa[[v]] - q_alpha[[v]])), 3)), "\n")
qn <- c(jaccard = "jaccard", braycurtis = "bray_curtis", unweighted_unifrac = "unweighted_unifrac", weighted_unifrac = "weighted_unifrac")
for (k in names(qn)) cat(k, "R vs QIIME 2, max|diff| =", signif(max(abs(as.matrix(pdm[[k]])[S, S] - qd(qn[[k]]))), 3), "\n")
write.csv(pa, file.path(O2, "primary_alpha.csv"), row.names = FALSE)
for (k in names(pdm)) write.csv(as.matrix(pdm[[k]])[S, S], file.path(O2, sprintf("primary_dm_%s.csv", k)))   # used by 03b
write.csv(alpha_tab(pa), file.path(O2, "primary_alpha_wilcoxon.csv"), row.names = FALSE)
write.csv(beta_tab(pdm), file.path(O2, "primary_beta_permanova_permdisp.csv"), row.names = FALSE)

## =============================================================================
## 3. 100 independent rarefactions, with and without rare-ASV filters
##    F1: ASVs present in >= 2 samples; F2: ASVs > 0.1% relative abundance in
##    at least one sample; F3: both
## =============================================================================
ra <- full / rowSums(full); prev <- colSums(full > 0); maxra <- apply(ra, 2, max)
FL <- list(F1_prev2 = prev >= 2, F2_ra0.1 = maxra > 0.001); FL$F3_both <- FL$F1_prev2 & FL$F2_ra0.1
fs <- do.call(rbind, lapply(names(FL), function(k) { m <- FL[[k]]; kept <- rowSums(full[, m]) / rowSums(full)
  data.frame(filter = k, n_ASV = sum(m), pct_reads_retained_median = 100 * median(kept), pct_reads_retained_min = 100 * min(kept)) }))
write.csv(fs, file.path(O3, "filter_summary.csv"), row.names = FALSE)

set.seed(SEED); RT <- lapply(1:N_RAREFY, function(i) rrarefy(full, DEPTH))      # rarefaction without replacement
sub <- function(X, m) { X <- X[, m, drop = FALSE]; X[, colSums(X) > 0, drop = FALSE] }
run_set <- function(tag, mask, rel_bc, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(mask)) {                                                  # filtered primary table
    Xp <- sub(rar_q, mask); ap <- alpha(Xp); Dp <- dms(Xp, rel_bc)
    write.csv(cbind(Filter = tag, alpha_tab(ap)), file.path(outdir, "primary_alpha_wilcoxon.csv"), row.names = FALSE)
    write.csv(cbind(Filter = tag, beta_tab(Dp)),  file.path(outdir, "primary_beta_permanova_permdisp.csv"), row.names = FALSE) }
  res <- mclapply(1:N_RAREFY, function(i) { X <- if (is.null(mask)) RT[[i]] else sub(RT[[i]], mask); a <- alpha(X); D <- dms(X, rel_bc)
    list(alpha = cbind(iter = i, a),
         beta = do.call(rbind, lapply(names(DL), function(k) data.frame(Distance = DL[[k]], iter = i, t(test_dm(D[[k]])))))) },
    mc.cores = N_CORES)
  A <- do.call(rbind, lapply(res, `[[`, "alpha")); B <- do.call(rbind, lapply(res, `[[`, "beta"))
  write.csv(A, file.path(outdir, "rarefaction100_alpha.csv"), row.names = FALSE)
  write.csv(cbind(Filter = tag, B), file.path(outdir, "rarefaction100_beta_tests.csv"), row.names = FALSE)
  as_ <- do.call(rbind, lapply(names(IDX), function(v) { p <- sapply(split(A, A$iter), function(d) wil(setNames(d[[v]], d$SampleID))$p.value)
    data.frame(Filter = tag, Index = IDX[[v]], P_median = median(p), P_min = min(p), P_max = max(p), n_P_lt_0.05 = sum(p < 0.05)) }))
  write.csv(as_, file.path(outdir, "rarefaction100_alpha_summary.csv"), row.names = FALSE)
  bs <- do.call(rbind, lapply(DL, function(k) { d <- B[B$Distance == k, ]
    data.frame(Filter = tag, Distance = k, R2_median = median(d$R2), R2_min = min(d$R2), R2_max = max(d$R2),
               P_median = median(d$P), P_min = min(d$P), P_max = max(d$P), pct_P_lt_0.05 = 100 * mean(d$P < 0.05),
               dispP_median = median(d$disp_P), dispP_min = min(d$disp_P), dispP_max = max(d$disp_P)) }))
  write.csv(bs, file.path(outdir, "rarefaction100_beta_summary.csv"), row.names = FALSE)
  cat("\n=====", tag, "=====\n"); print(as_[, -1], digits = 3); print(bs[, -1], digits = 3) }
run_set("ALL", NULL, FALSE, O2)
for (k in names(FL)) run_set(k, FL[[k]], TRUE, file.path(O3, k))

writeLines(capture.output(sessionInfo()), result_path("diversity", "sessionInfo.txt"))
