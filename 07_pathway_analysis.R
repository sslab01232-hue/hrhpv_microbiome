# Predicted MetaCyc pathways: models with age and richness adjustment, comparison
# with the 230-ASV analysis, taxon contributions, NSTI, and pathway detection
library(dplyr)
library(tidyr)
source("config.R")
source("R/maaslin_lm.R")

out  <- out_dir("pathways")
p1288 <- "picrust2_out_1288"
p230  <- "picrust2_out_230"

md <- read_metadata()
md$observed <- read.delim(data_path("core_metrics", "observed_features_vector", "alpha-diversity.tsv"),
                          row.names = 1)[rownames(md), 1]
md$logRich <- log(md$observed)

read_pw <- function(dir) {
  u <- read.delim(file.path(dir, "pathways_out", "path_abun_unstrat.tsv.gz"), row.names = 1, check.names = FALSE)
  t(as.matrix(u[, rownames(md)]))
}
ab <- read_pw(p1288)
desc <- read.delim(file.path(p1288, "pathways_out", "path_abun_unstrat_descrip.tsv.gz"), check.names = FALSE)
desc <- setNames(desc$description, desc[[1]])

# Pathway models (1,288 ASVs) -------------------------------------------------
models <- list(unadjusted = "Group", age = c("Group", "Age"), age_richness = c("Group", "Age", "logRich"))
res <- bind_rows(lapply(names(models), function(m) mutate(maaslin_lm(ab, md, models[[m]]), model = m)))
res$description <- desc[res$feature]
write.csv(res, file.path(out, "pathway_models_1288.csv"), row.names = FALSE)

grp <- filter(res, name == "Groupmultiple")
grp %>% group_by(model) %>%
  summarise(n = n(), pct_positive = 100 * mean(coef > 0), n_q0.1 = sum(qval < 0.1 & coef > 0), n_q0.05 = sum(qval < 0.05 & coef > 0))

sig_age  <- filter(grp, model == "age", qval < 0.1, coef > 0)$feature
sig_rich <- filter(grp, model == "age_richness", qval < 0.1, coef > 0)$feature
six <- intersect(sig_age, sig_rich)

# Pathways restricted to the 230 ASVs of the 20 most abundant species
ab230  <- read_pw(p230)
res230 <- bind_rows(lapply(names(models)[1:2], function(m) mutate(maaslin_lm(ab230, md, models[[m]]), model = m)))
write.csv(res230, file.path(out, "pathway_models_230.csv"), row.names = FALSE)
sig230 <- filter(res230, model == "age", name == "Groupmultiple", qval < 0.1)$feature
intersect(sig230, sig_age)

# Taxon contributions to the age-adjusted significant pathways -----------------
# ASVs with the same taxonomic label are summed; mean % of pathway abundance in
# samples in which the pathway was detected
asv <- read.delim(data_path("asv_relative_abundance_1288.tsv"), check.names = FALSE)
lab <- setNames(taxon_label(asv[[2]]), asv[[1]])
con <- read.delim(file.path(p1288, "pathways_out", "path_abun_contrib.tsv.gz"), check.names = FALSE) %>%
  filter(`function` %in% sig_age) %>%
  group_by(sample, `function`) %>%
  mutate(pct = 100 * taxon_function_abun / sum(taxon_function_abun), taxon = lab[taxon]) %>%
  group_by(sample, `function`, taxon) %>%
  summarise(pct = sum(pct), .groups = "drop") %>%
  complete(nesting(sample, `function`), taxon, fill = list(pct = 0)) %>%
  mutate(Group = md[sample, "Group"]) %>%
  group_by(`function`, taxon, Group) %>%
  summarise(mean_pct = mean(pct), n_samples = sum(pct > 0), .groups = "drop")
write.csv(con, file.path(out, "taxon_contributions.csv"), row.names = FALSE)

# NSTI -------------------------------------------------------------------------
nsti <- read.delim(Sys.glob(file.path(p1288, "*marker_predicted_and_nsti.tsv.gz"))[1])$metadata_NSTI
summary(nsti); sum(nsti > 2)
wn <- read.delim(file.path(p1288, "KO_metagenome_out", "weighted_nsti.tsv.gz"), row.names = 1)
wilcox.test(wn[rownames(md), 1] ~ md$Group, exact = FALSE)

# Detection --------------------------------------------------------------------
md$n_pathways <- rowSums(ab > 0)
cor.test(md$n_pathways, md$observed, method = "spearman", exact = FALSE)

detected <- ab[, six] > 0
fisher_p <- function(x, g) fisher.test(table(factor(x, c(TRUE, FALSE)), g))$p.value
det <- data.frame(pathway = six,
                  single = colSums(detected[md$Group == "single", ]),
                  multiple = colSums(detected[md$Group == "multiple", ]),
                  P = apply(detected, 2, fisher_p, g = md$Group))
det$q <- p.adjust(det$P, "BH")
write.csv(det, file.path(out, "pathway_detection.csv"), row.names = FALSE)

# detection by presence of Segatella copri
scopri <- rowsum(as.matrix(asv[, -(1:2)]), lab)["Segatella copri", rownames(md)] > 0
data.frame(pathway = six,
           detected_scopri_pos = colSums(detected[scopri, ]), detected_scopri_neg = colSums(detected[!scopri, ]),
           P = apply(detected, 2, fisher_p, g = factor(scopri, c(TRUE, FALSE))))

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
