## =============================================================================
## 04_network.R
## Co-occurrence networks (Supplementary Tables S15-S26; edge lists)
##
##  1. Network input: read counts of the 20 most abundant species per group
##  2. Primary analysis with NetCoMi 1.3.0 (SparCC, |r| >= 0.3, signed
##     dissimilarity, fast-greedy clustering, eigenvector hubs; netCompare with
##     1,000 permutations and BH adjustment; seed 2026)
##  3. Tables from the netCompare summary
##  4. Equivalence check of 05_network_functions.R against the NetCoMi output
##  5. Sensitivity analyses (05_network_functions.R):
##     thresholds 0.2/0.3/0.4, bootstrap edge stability, subsampling of the
##     single group to n = 14, exclusion of HPV52-positive participants
##
## Inputs : data/asv_table_1288.csv     ASV table with taxonomy (Supplementary
##                                      Table S2; column Species = Genus_species)
##          data/network_species.csv    20 species (column Taxa)
##          data/sample_metadata.csv
## Note   : sections 4-5 are computationally intensive (1,000 SparCC runs per
##          analysis); parallel::mclapply uses N_CORES (set N_CORES <- 1 on Windows).
## =============================================================================
suppressPackageStartupMessages({ library(NetCoMi); library(dplyr); library(igraph); library(parallel) })
source("config.R"); source("05_network_functions.R")
OUT <- file.path(RESULTS_DIR, "network"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
md <- read_metadata()

## =============================================================================
## 1. Network input (species-level read counts, samples x 20 species)
## =============================================================================
ds   <- read.csv(data_path("asv_table_1288.csv"), check.names = FALSE)
sp20 <- read.csv(data_path("network_species.csv"), check.names = FALSE)$Taxa
agg  <- ds %>% filter(Species %in% sp20) %>% group_by(Species) %>%
  summarise(across(all_of(md$SampleID), sum), .groups = "drop")
stopifnot(setequal(agg$Species, sp20))
agg  <- agg[match(sp20, agg$Species), ]
counts <- function(g) { s <- md$SampleID[md$Group == g]; m <- t(as.matrix(agg[, s])); colnames(m) <- agg$Species; m }
single <- counts("single"); multiple <- counts("multiple")
write.csv(t(single),   file.path(OUT, "input_counts_single.csv"))
write.csv(t(multiple), file.path(OUT, "input_counts_multiple.csv"))
cat("single:", nrow(single), "samples | multiple:", nrow(multiple), "samples | species:", ncol(single), "\n")

## =============================================================================
## 2. Primary analysis (NetCoMi)
## =============================================================================
set.seed(SEED)
net <- netConstruct(data = single, data2 = multiple, dataType = "counts",
                    measure = "sparcc", zeroMethod = "none", normMethod = "none",
                    sparsMethod = "threshold", thresh = 0.3, dissFunc = "signed", seed = SEED, verbose = 2)
ana <- netAnalyze(net, centrLCC = TRUE, clustMethod = "cluster_fast_greedy", hubPar = "eigenvector",
                  weightDeg = FALSE, normDeg = TRUE, verbose = 0)
cmp <- netCompare(ana, permTest = TRUE, nPerm = 1000, adjust = "BH", testRand = TRUE, nPermRand = 1000, gcd = TRUE,
                  cores = 1, seed = SEED, returnPermProps = TRUE, returnPermCentr = TRUE, verbose = TRUE)
saveRDS(net, file.path(OUT, "network_construct.rds")); saveRDS(ana, file.path(OUT, "network_analysis.rds"))
saveRDS(cmp, file.path(OUT, "network_compare.rds"))
writeLines(capture.output(summary(ana, groupNames = c("Single hrHPV", "Multiple hrHPV"))), file.path(OUT, "netAnalyze_summary.txt"))
writeLines(capture.output(summary(cmp, groupNames = c("Single hrHPV", "Multiple hrHPV"),
                                  showCentr = c("degree", "betweenness", "closeness", "eigenvector"), numbNodes = 20)),
           file.path(OUT, "netCompare_summary.txt"))

## edge lists (|SparCC| >= 0.3)
edges <- function(a, d, lab) { ut <- which(upper.tri(d) & d != 0, arr.ind = TRUE)
  e <- data.frame(taxon1 = rownames(a)[ut[, 1]], taxon2 = colnames(a)[ut[, 2]], sparcc = a[ut]); e <- e[order(-abs(e$sparcc)), ]
  write.csv(e, file.path(OUT, paste0("edges_", lab, ".csv")), row.names = FALSE)
  cat(lab, ": ", nrow(e), " edges (", sum(e$sparcc > 0), " positive)\n", sep = "") }
edges(net$assoMat1, net$adjaMat1, "single"); edges(net$assoMat2, net$adjaMat2, "multiple")
cat("hubs: single =", paste(ana$hubs$hubs1, collapse = ", "), "| multiple =", paste(ana$hubs$hubs2, collapse = ", "), "\n")

## =============================================================================
## 3. Tables from the netCompare summary
## =============================================================================
L <- readLines(file.path(OUT, "netCompare_summary.txt"))
sec <- function(start_pat, end_pat) { s <- grep(start_pat, L)[1]; e <- grep(end_pat, L); e <- e[e > s][1]; L[(s + 1):(e - 1)] }
parse_rows <- function(lines, ncol_num) {
  lines <- lines[nzchar(trimws(lines)) & !grepl("^\\s*(Single|-----|\\*|p-values)", lines)]
  do.call(rbind, lapply(lines, function(l) {
    tok <- unlist(strsplit(trimws(l), "\\s+"))
    if (length(tok) <= ncol_num || is.na(suppressWarnings(as.numeric(tok[length(tok)]))) && !grepl("^[.*]+$", tok[length(tok)])) return(NULL)
    tok <- tok[!grepl("^[.*]+$", tok)]                                   # drop significance codes
    k <- length(tok)
    data.frame(Metric = paste(tok[1:(k - ncol_num)], collapse = " "), t(as.numeric(tok[(k - ncol_num + 1):k])), stringsAsFactors = FALSE) }))
}
w <- parse_rows(sec("^Whole network:", "^-----"), 4); names(w)[2:5] <- c("Single hrHPV", "Multiple hrHPV", "Abs. diff.", "P-value")
five <- c("Clustering coefficient", "Modularity", "Positive edge percentage", "Edge density", "Natural connectivity")
w$`Adj. P-value (BH)` <- NA; i5 <- w$Metric %in% five
w$`Adj. P-value (BH)`[i5] <- p.adjust(w$`P-value`[i5], "BH")          # BH across the five reported properties
write.csv(w, file.path(OUT, "global_properties.csv"), row.names = FALSE); print(w)
j <- parse_rows(sec("^Jaccard index \\(similarity", "^-----"), 3); j <- j[!grepl("^Jacc", j$Metric), ]
names(j)[2:4] <- c("Jaccard index", "P(<= Jacc)", "P(>= Jacc)")
write.csv(j, file.path(OUT, "jaccard_index.csv"), row.names = FALSE)
ia <- grep("^ARI ", L); ig <- grep("^GCD ", L); tk <- function(i) unlist(strsplit(trimws(L[i]), "\\s+"))
write.csv(data.frame(Measure = c("Adjusted Rand index", "Graphlet correlation distance"),
                     Value = as.numeric(c(tk(ia)[2], tk(ig)[2])), P = as.numeric(c(tk(ia + 1)[2], tk(ig + 1)[2]))),
          file.path(OUT, "ari_gcd.csv"), row.names = FALSE)
cen <- function(title, fname) {
  s <- grep(paste0("^", title), L)[1]; e <- s + 1; while (e <= length(L) && nzchar(trimws(L[e]))) e <- e + 1
  d <- parse_rows(L[(s + 2):(e - 1)], 5); names(d) <- c("Taxon", "Single hrHPV", "Multiple hrHPV", "Abs. diff.", "P-value", "Adj. P-value")
  d$Taxon <- gsub("_", " ", d$Taxon); write.csv(d, file.path(OUT, fname), row.names = FALSE); d }
cen("Degree \\(normalized\\)", "centrality_degree.csv"); cen("Betweenness centrality \\(normalized\\)", "centrality_betweenness.csv")
cen("Closeness centrality \\(normalized\\)", "centrality_closeness.csv"); cen("Eigenvector centrality \\(normalized\\)", "centrality_eigenvector.csv")

## =============================================================================
## 4. Equivalence of 05_network_functions.R with the NetCoMi output
## =============================================================================
for (g in 1:2) {                                                         # properties from the NetCoMi SparCC estimates
  pr <- calc_props(build_net(net[[paste0("assoEst", g)]], 0.3))
  cat("group", g, "global properties max|diff| =",
      max(abs(unlist(pr[GLOB]) - unlist(ana$globalProps[paste0(GLOB, g)]))), "\n") }
set.seed(SEED); a1 <- sparcc(net$countMat1); a2 <- sparcc(net$countMat2)  # same RNG order as netConstruct
cat("SparCC max|diff|: single =", max(abs(a1 - net$assoEst1[rownames(a1), colnames(a1)])),
    "| multiple =", max(abs(a2 - net$assoEst2[rownames(a2), colnames(a2)])), "\n")
set.seed(SEED); pm <- perm_group_mat(nrow(net$countMat1), nrow(net$countMat2), 1000); seeds <- sample.int(1e8, size = 1000)
xb <- rbind(net$countMat1, net$countMat2)
for (p in 1:5) { set.seed(seeds[p]); b1 <- sparcc(xb[pm[p, ] == 1, ]); b2 <- sparcc(xb[pm[p, ] == 2, ])   # first permutations
  d <- abs(unlist(calc_props(build_net(b1, .3))[GLOB]) - unlist(calc_props(build_net(b2, .3))[GLOB]))
  cat("permutation", p, "max|diff| =", signif(max(abs(d - cmp$permDiffGlobal[p, c("diffClustCoef", "diffModul", "diffPEP", "diffDensity", "diffNatConnect")])), 3), "\n") }

## =============================================================================
## 5. Sensitivity analyses
## =============================================================================
c1 <- net$countMat1; c2 <- net$countMat2; A1 <- net$assoEst1; A2 <- net$assoEst2
TH <- c(0.2, 0.3, 0.4); B <- 1000

## permutation comparison of two count matrices at the given thresholds
compare_groups <- function(x1, x2, th_set) {
  set.seed(SEED); pm <- perm_group_mat(nrow(x1), nrow(x2), B); seeds <- sample.int(1e8, size = B); xb <- rbind(x1, x2)
  mclapply(1:B, function(p) { set.seed(seeds[p])
    b1 <- sparcc(xb[pm[p, ] == 1, , drop = FALSE]); b2 <- sparcc(xb[pm[p, ] == 2, , drop = FALSE])
    lapply(th_set, function(th) { p1 <- calc_props(build_net(b1, th)); p2 <- calc_props(build_net(b2, th))
      list(glob = abs(unlist(p1[GLOB]) - unlist(p2[GLOB])), cent = sapply(names(CENT), function(k) abs(p1[[k]] - p2[[k]]))) }) },
    mc.cores = N_CORES) }
summarise_perm <- function(perm, o1, o2, j, label) {
  obs <- abs(unlist(o1[GLOB]) - unlist(o2[GLOB])); pg <- t(sapply(perm, function(x) x[[j]]$glob))
  pv <- (colSums(sweep(pg, 2, obs, ">=")) + 1) / (B + 1)
  G <- data.frame(Analysis = label, Metric = GLOB_LAB[GLOB], Single = unlist(o1[GLOB]), Multiple = unlist(o2[GLOB]), Abs_diff = obs,
                  P = pv, P_adj_BH = p.adjust(pv, "BH"), Edges_single = o1$nEdges, Edges_multiple = o2$nEdges)
  C <- do.call(rbind, lapply(names(CENT), function(k) { ob <- abs(o1[[k]] - o2[[k]]); pc <- t(sapply(perm, function(x) x[[j]]$cent[, k]))
    pv <- (colSums(sweep(pc, 2, ob, ">=")) + 1) / (B + 1)
    data.frame(Analysis = label, Centrality = CENT[[k]], Taxon = names(ob), Single = o1[[k]], Multiple = o2[[k]],
               Abs_diff = ob, P = pv, P_adj_BH = p.adjust(pv, "BH")) }))
  list(global = G, centrality = C) }

## 5a. correlation thresholds 0.2 / 0.3 / 0.4 (Supplementary Tables S24-S25)
perm <- compare_groups(c1, c2, TH)
th_res <- lapply(seq_along(TH), function(j) summarise_perm(perm, calc_props(build_net(A1, TH[j])), calc_props(build_net(A2, TH[j])), j,
                                                            paste0("threshold_", TH[j])))
write.csv(do.call(rbind, lapply(th_res, `[[`, "global")),     file.path(OUT, "threshold_global_properties.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(th_res, `[[`, "centrality")), file.path(OUT, "threshold_centralities.csv"),      row.names = FALSE)

## 5b. bootstrap edge stability, 1,000 resamples with replacement per group (Supplementary Tables S22-S23)
boot <- list()
for (g in 1:2) {
  cm <- if (g == 1) c1 else c2; A <- if (g == 1) A1 else A2
  set.seed(SEED + g); idx <- lapply(1:B, function(b) sample(nrow(cm), replace = TRUE)); seeds <- sample.int(1e8, size = B)
  boots <- mclapply(1:B, function(b) { set.seed(seeds[b]); sparcc(cm[idx[[b]], , drop = FALSE]) }, mc.cores = N_CORES)
  for (th in TH) {
    ut <- which(upper.tri(A) & is_edge(A, th), arr.ind = TRUE); s <- sign(A[ut])
    ret <- matrix(sapply(boots, function(bm) is_edge(bm[ut], th) & sign(bm[ut]) == s), nrow = nrow(ut))
    nb  <- sapply(boots, function(bm) sum(is_edge(bm, th)[upper.tri(bm)]))
    boot[[length(boot) + 1]] <- data.frame(Group = c("single", "multiple")[g], Threshold = th,
      Taxon1 = rownames(A)[ut[, 1]], Taxon2 = colnames(A)[ut[, 2]], SparCC = A[ut], Retention = rowMeans(ret),
      Boot_edges_median = median(nb), Boot_edges_q025 = unname(quantile(nb, .025)), Boot_edges_q975 = unname(quantile(nb, .975)))
  }
}
boot <- do.call(rbind, boot); write.csv(boot, file.path(OUT, "bootstrap_edge_retention.csv"), row.names = FALSE)
print(aggregate(Retention ~ Group + Threshold, boot, function(z) c(n = length(z), median = median(z), pct_ge50 = 100 * mean(z >= .5))))

## 5c. single group subsampled to n = 14 without replacement, 1,000 times (Supplementary Table S21)
n2 <- nrow(c2); set.seed(SEED); idx <- lapply(1:B, function(b) sort(sample(nrow(c1), n2, replace = FALSE))); seeds <- sample.int(1e8, size = B)
ss <- do.call(rbind, mclapply(1:B, function(b) { set.seed(seeds[b]); a <- sparcc(c1[idx[[b]], , drop = FALSE])
  do.call(rbind, lapply(TH, function(th) { p <- calc_props(build_net(a, th))
    data.frame(iter = b, Threshold = th, t(unlist(p[GLOB])), nEdges = p$nEdges,
               deg_Gardnerella_vaginalis = p$deg[["Gardnerella_vaginalis"]], deg_Segatella_copri = p$deg[["Segatella_copri"]]) })) },
  mc.cores = N_CORES))
write.csv(ss, file.path(OUT, "subsample_single_n14.csv"), row.names = FALSE)
met <- c(GLOB, "nEdges", "deg_Gardnerella_vaginalis", "deg_Segatella_copri")
vals <- function(o) c(unlist(o[GLOB]), nEdges = o$nEdges, deg_Gardnerella_vaginalis = o$deg[["Gardnerella_vaginalis"]],
                      deg_Segatella_copri = o$deg[["Segatella_copri"]])
sub_sum <- do.call(rbind, lapply(TH, function(th) { sv <- vals(calc_props(build_net(A1, th))); mv <- vals(calc_props(build_net(A2, th)))
  d <- ss[ss$Threshold == th, ]
  do.call(rbind, lapply(met, function(m) { z <- d[[m]]
    data.frame(Threshold = th, Metric = m, Single_full_n23 = sv[[m]], Sub_median = median(z), Sub_q025 = unname(quantile(z, .025)),
               Sub_q975 = unname(quantile(z, .975)), Multiple_n14 = mv[[m]],
               Percentile_of_multiple = 100 * mean(z < mv[[m]]) + 50 * mean(z == mv[[m]]),
               P_ge = (sum(z >= mv[[m]]) + 1) / (B + 1), P_le = (sum(z <= mv[[m]]) + 1) / (B + 1)) })) }))
write.csv(sub_sum, file.path(OUT, "subsample_single_n14_summary.csv"), row.names = FALSE)
print(sub_sum[sub_sum$Metric %in% c("density", "natConnect", "nEdges"), ], digits = 3)

## 5d. HPV52-positive participants excluded (single n = 22, multiple n = 10; Supplementary Table S26)
pos <- md$SampleID[md$HPV52 == "Y"]
x1 <- c1[!rownames(c1) %in% pos, ]; x2 <- c2[!rownames(c2) %in% pos, ]
set.seed(SEED); B1 <- sparcc(x1); B2 <- sparcc(x2)                       # same RNG order as netConstruct
r52 <- summarise_perm(compare_groups(x1, x2, 0.3), calc_props(build_net(B1, 0.3)), calc_props(build_net(B2, 0.3)), 1, "HPV52_excluded")
write.csv(r52$global,     file.path(OUT, "network_HPV52_excluded_global.csv"),       row.names = FALSE)
write.csv(r52$centrality, file.path(OUT, "network_HPV52_excluded_centralities.csv"), row.names = FALSE)
print(r52$global, digits = 3)

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
