# Co-occurrence networks (NetCoMi 1.3.0) and sensitivity analyses
library(NetCoMi)
library(igraph)
library(parallel)
source("config.R")
source("05_network_functions.R")

out <- out_dir("network")
md  <- read_metadata()

# Read counts of the 20 most abundant species ---------------------------------
asv  <- read.csv(data_path("asv_table_1288.csv"), check.names = FALSE)
sp20 <- read.csv(data_path("network_species.csv"))$Taxa
in20 <- asv$Species %in% sp20
counts   <- rowsum(as.matrix(asv[in20, md$SampleID]), asv$Species[in20])[sp20, ]
single   <- t(counts[, md$Group == "single"])
multiple <- t(counts[, md$Group == "multiple"])

# Primary networks --------------------------------------------------------------
set.seed(seed)
net <- netConstruct(data = single, data2 = multiple, dataType = "counts", measure = "sparcc",
                    zeroMethod = "none", normMethod = "none", sparsMethod = "threshold", thresh = 0.3,
                    dissFunc = "signed", seed = seed)
ana <- netAnalyze(net, centrLCC = TRUE, clustMethod = "cluster_fast_greedy", hubPar = "eigenvector",
                  weightDeg = FALSE, normDeg = TRUE)
cmp <- netCompare(ana, permTest = TRUE, nPerm = 1000, adjust = "BH", testRand = TRUE, nPermRand = 1000,
                  gcd = TRUE, seed = seed, returnPermProps = TRUE, returnPermCentr = TRUE)
saveRDS(list(net = net, ana = ana, cmp = cmp), file.path(out, "netcomi.rds"))
groups <- c("Single hrHPV", "Multiple hrHPV")
writeLines(capture.output(summary(cmp, groupNames = groups, numbNodes = 20,
                                  showCentr = c("degree", "betweenness", "closeness", "eigenvector"))),
           file.path(out, "netCompare_summary.txt"))

# edge lists
edge_list <- function(a, adj) {
  ut <- which(upper.tri(adj) & adj != 0, arr.ind = TRUE)
  data.frame(taxon1 = rownames(a)[ut[, 1]], taxon2 = colnames(a)[ut[, 2]], sparcc = a[ut])
}
write.csv(edge_list(net$assoMat1, net$adjaMat1), file.path(out, "edges_single.csv"), row.names = FALSE)
write.csv(edge_list(net$assoMat2, net$adjaMat2), file.path(out, "edges_multiple.csv"), row.names = FALSE)

# Tables from the netCompare summary -------------------------------------------
txt <- readLines(file.path(out, "netCompare_summary.txt"))
block <- function(title, n_num, stop = "^\\s*$") {          # rows from a title to the next stop line
  from <- grep(title, txt)[1] + 1
  rest <- txt[from:length(txt)]
  rows <- rest[seq_len(grep(stop, rest)[1] - 1)]
  rows <- sub("\\s+[.*]+\\s*$", "", rows)                       # drop significance codes
  m <- regmatches(rows, regexec(paste0("^\\s*(.*?)", strrep("\\s+(\\S+)", n_num), "\\s*$"), rows))
  m <- do.call(rbind, m[lengths(m) == n_num + 2])
  v <- suppressWarnings(matrix(as.numeric(m[, -(1:2)]), nrow(m)))
  data.frame(name = m[, 2], v)[complete.cases(v), ]
}
glob <- setNames(block("^Whole network:", 4, stop = "^-{3,}"), c("Property", groups, "Abs_diff", "P"))
sel  <- glob$Property %in% GLOB
glob$q <- NA
glob$q[sel] <- p.adjust(glob$P[sel], "BH")                       # BH across the five reported properties
write.csv(glob, file.path(out, "global_properties.csv"), row.names = FALSE)

cent <- do.call(rbind, lapply(c("Degree", "Betweenness centrality", "Closeness centrality", "Eigenvector centrality"),
  function(k) cbind(Centrality = k, setNames(block(paste0("^", k, " \\(normalized\\)"), 5),
                                             c("Taxon", groups, "Abs_diff", "P", "q")))))
write.csv(cent, file.path(out, "centralities.csv"), row.names = FALSE)
i <- grep("^(ARI|GCD) ", txt)
print(txt[sort(c(i, i + 1))])

# The functions in 05_network_functions.R reproduce the NetCoMi estimates ------
set.seed(seed); a1 <- sparcc(net$countMat1); a2 <- sparcc(net$countMat2)
print(c(max(abs(a1 - net$assoEst1[rownames(a1), colnames(a1)])), max(abs(a2 - net$assoEst2[rownames(a2), colnames(a2)]))))
set.seed(seed); pm <- perm_group_mat(nrow(net$countMat1), nrow(net$countMat2), 1000); ps <- sample.int(1e8, 1000)
xb <- rbind(net$countMat1, net$countMat2)
print(sapply(1:5, function(p) {
  set.seed(ps[p])
  b1 <- calc_props(build_net(sparcc(xb[pm[p, ] == 1, ]), 0.3))
  b2 <- calc_props(build_net(sparcc(xb[pm[p, ] == 2, ]), 0.3))
  max(abs(abs(unlist(b1[names(GLOB)]) - unlist(b2[names(GLOB)])) - cmp$permDiffGlobal[p, 1:5]))
}))

# Sensitivity analyses ---------------------------------------------------------
c1 <- net$countMat1; c2 <- net$countMat2; A1 <- net$assoEst1; A2 <- net$assoEst2
TH <- c(0.2, 0.3, 0.4)
B  <- 1000

# permutation test (as in netCompare) at one or more thresholds
perm_compare <- function(x1, x2, th) {
  set.seed(seed)
  pm <- perm_group_mat(nrow(x1), nrow(x2), B)
  ps <- sample.int(1e8, B)
  xb <- rbind(x1, x2)
  mclapply(seq_len(B), function(p) {
    set.seed(ps[p])
    b1 <- sparcc(xb[pm[p, ] == 1, , drop = FALSE])
    b2 <- sparcc(xb[pm[p, ] == 2, , drop = FALSE])
    lapply(th, function(t) {
      p1 <- calc_props(build_net(b1, t)); p2 <- calc_props(build_net(b2, t))
      list(glob = abs(unlist(p1[names(GLOB)]) - unlist(p2[names(GLOB)])),
           cent = sapply(names(CENT), function(k) abs(p1[[k]] - p2[[k]])))
    })
  }, mc.cores = n_cores)
}
perm_summary <- function(perm, o1, o2, j, label) {
  obs <- abs(unlist(o1[names(GLOB)]) - unlist(o2[names(GLOB)]))
  pg  <- t(sapply(perm, function(x) x[[j]]$glob))
  p   <- (colSums(sweep(pg, 2, obs, ">=")) + 1) / (B + 1)
  g <- data.frame(Analysis = label, Property = GLOB, Single = unlist(o1[names(GLOB)]), Multiple = unlist(o2[names(GLOB)]),
                  P = p, q = p.adjust(p, "BH"), Edges_single = o1$nEdges, Edges_multiple = o2$nEdges)
  ce <- do.call(rbind, lapply(names(CENT), function(k) {
    ob <- abs(o1[[k]] - o2[[k]])
    pc <- t(sapply(perm, function(x) x[[j]]$cent[, k]))
    pv <- (colSums(sweep(pc, 2, ob, ">=")) + 1) / (B + 1)
    data.frame(Analysis = label, Centrality = CENT[[k]], Taxon = names(ob), P = pv, q = p.adjust(pv, "BH"))
  }))
  list(global = g, centrality = ce)
}

# correlation thresholds 0.2, 0.3, 0.4
perm <- perm_compare(c1, c2, TH)
thr  <- Map(function(t, j) perm_summary(perm, calc_props(build_net(A1, t)), calc_props(build_net(A2, t)), j,
                                        paste("threshold", t)), TH, seq_along(TH))
write.csv(do.call(rbind, lapply(thr, `[[`, "global")),     file.path(out, "threshold_global.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(thr, `[[`, "centrality")), file.path(out, "threshold_centrality.csv"), row.names = FALSE)

# bootstrap edge stability (1,000 resamples per group)
boots <- lapply(1:2, function(g) {
  cm <- list(c1, c2)[[g]]
  set.seed(seed + g)
  idx <- lapply(seq_len(B), function(b) sample(nrow(cm), replace = TRUE))
  bs  <- sample.int(1e8, B)
  mclapply(seq_len(B), function(b) { set.seed(bs[b]); sparcc(cm[idx[[b]], , drop = FALSE]) }, mc.cores = n_cores)
})
boot <- do.call(rbind, Map(function(g, t) {
  A  <- list(A1, A2)[[g]]
  ut <- which(upper.tri(A) & is_edge(A, t), arr.ind = TRUE)
  kept <- sapply(boots[[g]], function(bm) is_edge(bm[ut], t) & sign(bm[ut]) == sign(A[ut]))
  data.frame(Group = c("single", "multiple")[g], Threshold = t, Taxon1 = rownames(A)[ut[, 1]],
             Taxon2 = colnames(A)[ut[, 2]], SparCC = A[ut], Retention = rowMeans(matrix(kept, nrow(ut))))
}, rep(1:2, each = length(TH)), rep(TH, 2)))
write.csv(boot, file.path(out, "bootstrap_edge_retention.csv"), row.names = FALSE)
aggregate(Retention ~ Group + Threshold, boot, function(r) 100 * mean(r >= 0.5))

# single group subsampled to the size of the multiple group (1,000 times)
set.seed(seed)
sub_idx <- lapply(seq_len(B), function(b) sort(sample(nrow(c1), nrow(c2))))
sub_seeds <- sample.int(1e8, B)
sub <- do.call(rbind, mclapply(seq_len(B), function(b) {
  set.seed(sub_seeds[b])
  a <- sparcc(c1[sub_idx[[b]], , drop = FALSE])
  do.call(rbind, lapply(TH, function(t) {
    p <- calc_props(build_net(a, t))
    data.frame(iter = b, Threshold = t, t(unlist(p[names(GLOB)])), nEdges = p$nEdges)
  }))
}, mc.cores = n_cores))
write.csv(sub, file.path(out, "subsample_single_n14.csv"), row.names = FALSE)
sub_summary <- do.call(rbind, lapply(TH, function(t) {
  obs <- unlist(calc_props(build_net(A2, t))[c("density", "natConnect", "nEdges")])
  s <- as.matrix(sub[sub$Threshold == t, c("density", "natConnect", "nEdges")])
  data.frame(Threshold = t, Metric = names(obs), Multiple = obs,
             Sub_median = apply(s, 2, median), Sub_q025 = apply(s, 2, quantile, 0.025), Sub_q975 = apply(s, 2, quantile, 0.975),
             P = (colSums(sweep(s, 2, obs, ">=")) + 1) / (B + 1))
}))
write.csv(sub_summary, file.path(out, "subsample_summary.csv"), row.names = FALSE)

# HPV52-positive participants excluded
pos <- md$SampleID[md$HPV52 == "Y"]
x1  <- c1[!rownames(c1) %in% pos, ]
x2  <- c2[!rownames(c2) %in% pos, ]
set.seed(seed); B1 <- sparcc(x1); B2 <- sparcc(x2)
hpv52 <- perm_summary(perm_compare(x1, x2, 0.3), calc_props(build_net(B1, 0.3)), calc_props(build_net(B2, 0.3)),
                      1, "HPV52 excluded")
write.csv(hpv52$global, file.path(out, "HPV52_excluded_global.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
