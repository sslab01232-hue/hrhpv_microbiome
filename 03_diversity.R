# Alpha and beta diversity: primary analysis (40,000 reads), 100 rarefactions,
# rare-ASV filters, and rarefaction curves
library(vegan)
library(ape)
library(Matrix)
library(parallel)
source("config.R")

out  <- out_dir("diversity")
md   <- read_metadata()
rar  <- as.matrix(read_biom_tsv(data_path("rarefied_table_40000.tsv")))   # QIIME 2 core-metrics table
S    <- colnames(rar)                                                    # sample order of the QIIME 2 output
full <- t(as.matrix(read_biom_tsv(data_path("feature_table_1288.tsv")))[, S])
rar  <- t(rar[colnames(full), S])
grp  <- setNames(md[S, "Group"], S)
tree <- read.tree(data_path("tree.nwk"))

# Phylogenetic metrics -------------------------------------------------------
# branch x tip incidence: each tip is linked to the branches on its path to the root
paths  <- nodepath(tree)
edge_i <- match(seq_len(max(tree$edge)), tree$edge[, 2])
inc <- sparseMatrix(i = unlist(lapply(paths, function(p) edge_i[p[-1]])),
                    j = rep(seq_along(paths), lengths(paths) - 1),
                    dims = c(nrow(tree$edge), Ntip(tree)), dimnames = list(NULL, tree$tip.label))
blen <- tree$edge.length

branch <- function(X) as.matrix(inc[, colnames(X)] %*% t(X))   # branches x samples
faith  <- function(X) colSums((branch(X > 0) > 0) * blen)

unifrac_uw <- function(X) {
  B <- (branch(X > 0) > 0) * 1
  shared <- crossprod(B, B * blen)
  total  <- colSums(B * blen)
  as.dist(1 - shared / (outer(total, total, "+") - shared))
}
unifrac_w <- function(X) dist(t(branch(X / rowSums(X)) * blen), "manhattan")   # un-normalized, as in QIIME 2

alpha <- function(X) {
  sh <- diversity(X, "shannon", base = 2)
  ob <- rowSums(X > 0)
  data.frame(SampleID = rownames(X), Observed = ob, FaithPD = faith(X), Shannon = sh, Evenness = sh / log2(ob))
}
distances <- function(X, rel_bc = FALSE) {
  clr <- log(X + 1) - rowMeans(log(X + 1))
  list(jaccard            = vegdist(X > 0, "jaccard", binary = TRUE),
       braycurtis         = vegdist(if (rel_bc) X / rowSums(X) else X, "bray"),
       unweighted_unifrac = unifrac_uw(X),
       weighted_unifrac   = unifrac_w(X),
       aitchison          = dist(clr))
}

# Tests ----------------------------------------------------------------------
idx <- c(Observed = "Observed features", FaithPD = "Faith's PD", Shannon = "Shannon", Evenness = "Pielou's evenness")

alpha_test <- function(a) {
  r <- data.frame(Index = idx, P = sapply(names(idx), function(v)
    suppressWarnings(wilcox.test(a[[v]] ~ grp[a$SampleID]))$p.value))
  r$q <- p.adjust(r$P, "BH")
  r
}
beta_test <- function(D) {
  r <- do.call(rbind, lapply(names(D), function(k) {
    d <- D[[k]]; g <- grp[labels(d)]
    set.seed(seed); ad <- adonis2(d ~ g, permutations = n_perm)
    set.seed(seed); pd <- permutest(betadisper(d, g, type = "median"), permutations = n_perm)
    data.frame(Distance = k, R2 = ad$R2[1], F = ad$F[1], P = ad$`Pr(>F)`[1], PERMDISP_P = pd$tab$`Pr(>F)`[1])
  }))
  r$q <- p.adjust(r$P, "BH")
  r
}

# Read depth -----------------------------------------------------------------
reads <- rowSums(full)
write.csv(data.frame(SampleID = S, Group = grp, reads), file.path(out, "read_depth.csv"), row.names = FALSE)
print(tapply(reads, grp, median))
print(wilcox.test(reads ~ grp)$p.value)

# Rarefaction curves (Supplementary Figure S1): 10 subsamples per depth
steps <- sort(unique(c(500, 1000, 2500, 5000, seq(7500, 92500, 2500), depth, min(reads))))
set.seed(seed)
curve <- do.call(rbind, lapply(S, function(s) {
  x <- full[s, , drop = FALSE]
  do.call(rbind, lapply(steps[steps <= sum(x)], function(dp) {
    R <- rrarefy(x[rep(1, 10), , drop = FALSE], dp)
    rownames(R) <- paste0(s, "_", 1:10)
    data.frame(SampleID = s, depth = dp, observed = mean(rowSums(R > 0)), faith_pd = mean(faith(R)),
               goods_coverage = mean(1 - rowSums(R == 1) / dp))
  }))
}))
write.csv(curve, file.path(out, "rarefaction_curves.csv"), row.names = FALSE)

# Primary analysis on the QIIME 2 rarefied table ----------------------------
pa  <- alpha(rar)
pdm <- distances(rar)

# agreement with QIIME 2 core-metrics output
cm <- function(...) file.path(data_path("core_metrics"), ...)
q_alpha <- sapply(c(Observed = "observed_features", FaithPD = "faith_pd", Shannon = "shannon", Evenness = "evenness"),
                  function(n) read.delim(cm(paste0(n, "_vector"), "alpha-diversity.tsv"), row.names = 1)[S, 1])
q_dist  <- c(jaccard = "jaccard", braycurtis = "bray_curtis",
             unweighted_unifrac = "unweighted_unifrac", weighted_unifrac = "weighted_unifrac")
print(colSums(abs(as.matrix(pa[, colnames(q_alpha)]) - q_alpha)))
print(sapply(names(q_dist), function(k) {
  q <- as.matrix(read.delim(cm(paste0(q_dist[[k]], "_distance_matrix"), "distance-matrix.tsv"),
                            row.names = 1, check.names = FALSE))
  max(abs(as.matrix(pdm[[k]])[S, S] - q[S, S]))
}))

write.csv(pa, file.path(out, "primary_alpha.csv"), row.names = FALSE)
invisible(lapply(names(pdm), function(k)
  write.csv(as.matrix(pdm[[k]]), file.path(out, paste0("primary_dm_", k, ".csv")))))
write.csv(alpha_test(pa), file.path(out, "primary_alpha_tests.csv"), row.names = FALSE)
write.csv(beta_test(pdm), file.path(out, "primary_beta_tests.csv"), row.names = FALSE)

# 100 rarefactions, with and without rare-ASV filters ------------------------
prev  <- colSums(full > 0)
maxra <- apply(full / rowSums(full), 2, max)
filters <- list(all = NULL, prev2 = prev >= 2, ra0.1 = maxra > 0.001, both = prev >= 2 & maxra > 0.001)

set.seed(seed)
rarefied <- lapply(seq_len(n_rarefy), function(i) rrarefy(full, depth))

keep_asv <- function(X, m) if (is.null(m)) X else { X <- X[, m, drop = FALSE]; X[, colSums(X) > 0, drop = FALSE] }

run_filter <- function(name) {
  m <- filters[[name]]
  rel_bc <- !is.null(m)
  if (rel_bc) {                                                  # filtered primary table
    X <- keep_asv(rar, m)
    write.csv(alpha_test(alpha(X)), file.path(out, paste0("filter_", name, "_alpha_tests.csv")), row.names = FALSE)
    write.csv(beta_test(distances(X, TRUE)), file.path(out, paste0("filter_", name, "_beta_tests.csv")), row.names = FALSE)
  }
  res <- mclapply(seq_len(n_rarefy), function(i) {
    X <- keep_asv(rarefied[[i]], m)
    list(alpha = alpha_test(alpha(X)), beta = beta_test(distances(X, rel_bc)))
  }, mc.cores = n_cores)
  a <- do.call(rbind, Map(cbind, iter = seq_len(n_rarefy), lapply(res, `[[`, "alpha")))
  b <- do.call(rbind, Map(cbind, iter = seq_len(n_rarefy), lapply(res, `[[`, "beta")))
  write.csv(a, file.path(out, paste0("rarefy100_", name, "_alpha.csv")), row.names = FALSE)
  write.csv(b, file.path(out, paste0("rarefy100_", name, "_beta.csv")), row.names = FALSE)
  aggregate(P ~ Index, a, function(p) c(median = median(p), min = min(p), max = max(p), n_sig = sum(p < 0.05)))
}
summaries <- lapply(setNames(names(filters), names(filters)), run_filter)
print(summaries)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
