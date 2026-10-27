## =============================================================================
## 05_network_functions.R
## Functions for the network sensitivity analyses (sourced by 04_network.R).
##
## The functions transcribe the NetCoMi 1.3.0 / SpiecEasi pipeline used for the
## primary network analysis, so that the same network construction and
## comparison can be repeated under resampling:
##  - SparCC (SpiecEasi::sparcc; iter = 100, inner_iter = 20, th = 0.1)
##  - threshold sparsification, "signed" dissimilarity, similarity = 1 - dissimilarity
##  - netAnalyze: centrLCC = TRUE, fast-greedy clustering, weighted clustering
##    coefficient (Barrat), normalized natural connectivity, normalized centralities
##  - netCompare: permutation test on absolute differences,
##    P = (#{perm >= obs} + 1) / (nPerm + 1), BH adjustment
## Equivalence with the NetCoMi output is checked in 04_network.R (section 4).
## =============================================================================
suppressPackageStartupMessages({ library(igraph); library(VGAM) })
## ---- SpiecEasi::sparcc (verbatim logic) ----
norm_to_total <- function(x) x / sum(x)
clr_vec <- function(x.f, base = exp(1), tol = .Machine$double.eps) {
  nzero <- (x.f >= tol); LOG <- log(ifelse(nzero, x.f, 1), base); ifelse(nzero, LOG - mean(LOG) / mean(nzero), 0.0) }
clr_mat <- function(x.f, mar = 2) apply(x.f, mar, clr_vec)
cor2cov <- function(cor, sds) cor * sds * rep(sds, each = nrow(cor))
norm_diric <- function(x, rep = 1) { dmat <- VGAM::rdiric(rep, x + 1); norm_to_total(colMeans(dmat)) }
av <- function(data) { cov.clr <- cov(clr_mat(data)); J <- matrix(1, ncol(data), ncol(data))
  (J %*% diag(diag(cov.clr))) + (diag(diag(cov.clr)) %*% J) - (2 * cov.clr) }
basis_var <- function(T_mat, CovMat = matrix(0, nrow(T_mat), ncol(T_mat)),
                      M = matrix(1, nrow(T_mat), ncol(T_mat)) + (diag(ncol(T_mat)) * (ncol(T_mat) - 2)),
                      excluded = NULL, Vmin = 1e-4) {
  if (!is.null(excluded)) T_mat[excluded] <- 0
  Ti <- matrix(rowSums(T_mat)); CovVec <- matrix(rowSums(CovMat - diag(diag(CovMat))))
  M.I <- tryCatch(solve(M), error = function(e) MASS::ginv(M))
  Vbase <- M.I %*% (Ti + 2 * CovVec); Vbase[Vbase < Vmin] <- Vmin; list(Vbase = Vbase, M = M) }
C_from_V <- function(T_mat, Vbase) {
  J <- matrix(1, nrow(T_mat), ncol(T_mat)); Vdiag <- diag(c(Vbase))
  CovMat <- .5 * ((J %*% Vdiag) + (Vdiag %*% J) - T_mat); CovMat <- (CovMat + t(CovMat)) / 2
  CorMat <- cov2cor(CovMat); CorMat[abs(CorMat) > 1] <- sign(CorMat[abs(CorMat) > 1])
  CovMat <- cor2cov(CorMat, sqrt(as.vector(Vbase))); list(Cov = CovMat, Cor = CorMat) }
exclude_pairs <- function(Cor, M, th = 0.1, excluded = NULL) {
  break_flag <- FALSE; C_temp <- abs(Cor - diag(diag(Cor)))
  if (!is.null(excluded)) C_temp[excluded] <- 0
  exclude <- which(abs(C_temp - max(C_temp)) < .Machine$double.eps * 100)[seq_len(2)]
  if (max(C_temp) > th) { i <- na.exclude(arrayInd(exclude, c(nrow(M), ncol(M)))[, 1]); M[i, i] <- M[i, i] - 1
    excluded_new <- c(excluded, exclude) } else { excluded_new <- excluded; break_flag <- TRUE }
  list(M = M, excluded = excluded_new, break_flag = break_flag) }
sparccinner <- function(data.f, T_mat = NULL, iter = 10, th = 0.1) {
  if (is.null(T_mat)) T_mat <- av(data.f)
  res.bv <- basis_var(T_mat); Vbase <- res.bv$Vbase; M <- res.bv$M
  cbase <- C_from_V(T_mat, Vbase); Cov <- cbase$Cov; Cor <- cbase$Cor; excluded <- NULL
  for (i in seq_len(iter)) {
    res.excl <- exclude_pairs(Cor, M, th, excluded); M <- res.excl$M; excluded <- res.excl$excluded
    if (res.excl$break_flag) break
    res.bv <- basis_var(T_mat, M = M, excluded = excluded); Vbase <- res.bv$Vbase; M <- res.bv$M
    cbase <- C_from_V(T_mat, Vbase); Cov <- cbase$Cov; Cor <- cbase$Cor }
  list(Cov = Cov, Cor = Cor) }
sparcc <- function(data, iter = 100, inner_iter = 20, th = .1) {   # NetCoMi defaults
  sparccs <- lapply(1:iter, function(i) sparccinner(t(apply(data, 1, norm_diric)), iter = inner_iter, th = th))
  cors <- array(unlist(lapply(sparccs, function(x) x$Cor)), c(ncol(data), ncol(data), iter))
  corMed <- apply(cors, 1:2, median); dimnames(corMed) <- list(colnames(data), colnames(data)); corMed }
## ---- netConstruct (threshold + signed) ----
build_net <- function(asso, thresh) {
  a <- asso; a[abs(a) < thresh] <- 0
  x <- a; xv <- x[lower.tri(x)]; dv <- sqrt(0.5 * (1 - xv)); d <- x; d[lower.tri(d)] <- dv
  d[upper.tri(d)] <- t(d)[upper.tri(t(d))]; d[x == 0] <- Inf; diag(d) <- 0
  s <- 1 - d; s[d == Inf] <- 0
  list(asso = a, diss = d, adja = s) }
natural_connectivity <- function(G, norm = TRUE) {   # pulsar::natural.connectivity
  eig <- eigen(G); estrind <- exp(eig$values); nc <- log(mean(estrind))
  if (norm) { n <- length(estrind); nc <- nc / (n - log(n)) }; nc }
## ---- netAnalyze (.calcProps subset; NetCoMi defaults + centrLCC = TRUE) ----
calc_props <- function(nw) {
  adjaMat <- nw$adja; dissMat <- nw$diss; assoMat <- nw$asso; nNodes <- ncol(adjaMat); nm <- colnames(adjaMat)
  zero <- setNames(numeric(nNodes), nm)
  if (all(adjaMat[lower.tri(adjaMat)] == 0))
    return(list(clustCoef = 0, modularity = 0, pep = 0, density = 0, natConnect = 0, nEdges = 0,
                deg = zero, betw = zero, close = zero, eigen = zero))
  net <- graph_from_adjacency_matrix(adjaMat, weighted = TRUE, mode = "undirected", diag = FALSE)
  dg <- decompose(net); net_lcc <- dg[[which.max(sapply(dg, vcount))]]; lccSize <- vcount(net_lcc)
  lccNames <- V(net_lcc)$name
  dm_lcc <- dissMat[lccNames, lccNames]; dni <- dissMat; dni[is.infinite(dni)] <- 0
  dni_lcc <- dm_lcc; dni_lcc[is.infinite(dni_lcc)] <- 0
  dissnet_lcc <- graph_from_adjacency_matrix(dni_lcc, weighted = TRUE, mode = "undirected", diag = FALSE)
  cr <- cluster_fast_greedy(net); clust <- setNames(cr$membership, cr$names)
  ct <- table(clust); clust[clust %in% which(ct == 1)] <- 0
  sPath_lcc <- distances(dissnet_lcc, algorithm = "dijkstra")
  dv_lcc <- dm_lcc[lower.tri(dm_lcc)]; dv_lcc[is.infinite(dv_lcc)] <- 1; avDiss_lcc <- mean(dv_lcc)
  sPath_lcc <- sPath_lcc / avDiss_lcc
  natC <- natural_connectivity(adjaMat, norm = TRUE)
  cc <- mean(transitivity(net, type = "barrat"), na.rm = TRUE); if (is.nan(cc)) cc <- 0
  modul <- modularity(net, clust + 1, weights = E(net)$weight)   # igraph >= 2 uses edge weights by default (as in the original run)
  dens <- edge_density(net)
  lt <- assoMat[lower.tri(assoMat)]; pep <- sum(lt > 0) / sum(lt != 0) * 100
  deg <- degree(net, normalized = TRUE); deg[!names(deg) %in% lccNames] <- 0
  betw <- zero; bt <- betweenness(dissnet_lcc, normalized = TRUE); betw[names(bt)] <- bt; betw[is.nan(betw)] <- 0
  close <- zero; sr <- 1 / sPath_lcc; diag(sr) <- NA; cl <- rowSums(sr, na.rm = TRUE); close[lccNames] <- cl / (lccSize - 1)
  eig <- zero; ev <- eigen_centrality(net_lcc)$vector; eig[names(ev)] <- ev
  list(clustCoef = cc, modularity = modul, pep = pep, density = dens, natConnect = natC, nEdges = sum(adjaMat[lower.tri(adjaMat)] != 0),   # edges of the graph (pairs with SparCC = -1 have similarity 0 and are not edges)
       deg = deg[nm], betw = betw[nm], close = close[nm], eigen = eig[nm]) }
GLOB <- c("clustCoef", "modularity", "pep", "density", "natConnect")
GLOB_LAB <- c(clustCoef = "Clustering coefficient", modularity = "Modularity", pep = "Positive edge percentage",
              density = "Edge density", natConnect = "Natural connectivity")
CENT <- c(deg = "Degree", betw = "Betweenness", close = "Closeness", eigen = "Eigenvector")
## ---- netCompare permutation (RNG sequence as in NetCoMi: set.seed(seed); permuted labels; per-permutation seeds) ----
perm_group_mat <- function(n1, n2, nPerm) {
  groupvec <- c(rep(1, n1), rep(2, n2))
  clean <- function(m) { m <- unique.matrix(m); o <- sapply(1:nrow(m), function(i) identical(m[i, ], groupvec)); if (any(o)) m <- m[!o, ]; m }
  m <- t(sapply(1:(nPerm + 10), function(x) sample(groupvec))); m <- clean(m)
  while (nrow(m) < nPerm) { m <- rbind(m, sample(groupvec)); m <- clean(m) }
  m[1:nPerm, ] }
perm_test <- function(count1, count2, threshs, nPerm = 1000, seed = 2026, asso1 = NULL, asso2 = NULL) {
  ## returns observed props and permutation absolute differences for every threshold (same permutations reused)
  set.seed(seed); pm <- perm_group_mat(nrow(count1), nrow(count2), nPerm); seeds <- sample.int(1e8, size = nPerm)
  xbind <- rbind(count1, count2)
  perm <- lapply(1:nPerm, function(p) {
    set.seed(seeds[p])
    a1 <- sparcc(xbind[pm[p, ] == 1, , drop = FALSE]); a2 <- sparcc(xbind[pm[p, ] == 2, , drop = FALSE])
    lapply(setNames(threshs, threshs), function(th) {
      p1 <- calc_props(build_net(a1, th)); p2 <- calc_props(build_net(a2, th))
      list(glob = abs(unlist(p1[GLOB]) - unlist(p2[GLOB])),
           cent = sapply(names(CENT), function(k) abs(p1[[k]] - p2[[k]]))) }) })
  perm }
## an edge requires |SparCC| >= threshold; a coefficient of exactly -1 has
## similarity 0 under the signed transformation and is therefore not an edge
is_edge <- function(m, th) abs(m) >= th & m > -1
