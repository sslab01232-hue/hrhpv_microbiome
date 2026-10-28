# Taxon-level differential abundance: MaAsLin2 1.18.0 and sensitivity analyses
library(dplyr)
library(Maaslin2)
source("config.R")
source("R/maaslin_lm.R")

out   <- out_dir("taxa")
meta  <- read_metadata()
abund <- species_table(data_path("asv_relative_abundance_1288.tsv"), data_path("top20_species.csv"))
abund <- as.data.frame(abund[rownames(meta), ])

# MaAsLin2 (TSS, LOG, LM; reference = single)
models <- list(unadjusted     = "Group",
               age_adjusted   = c("Group", "Age"),
               HPV52_adjusted = c("Group", "HPV52"))
refs   <- list(unadjusted = "Group,single", age_adjusted = "Group,single",
               HPV52_adjusted = c("Group,single", "HPV52,N"))

res <- bind_rows(lapply(names(models), function(m) {
  fit <- Maaslin2(abund, meta, output = file.path(out, m),
                  fixed_effects = models[[m]], reference = refs[[m]],
                  normalization = "TSS", transform = "LOG", analysis_method = "LM",
                  standardize = FALSE, min_prevalence = 0.1, min_abundance = 0,
                  plot_heatmap = FALSE, plot_scatter = FALSE)
  mutate(fit$results, model = m)
}))
write.csv(res, file.path(out, "maaslin2_all_models.csv"), row.names = FALSE)

# maaslin_lm() reproduces the MaAsLin2 estimates
res$feature <- gsub("\\.", " ", res$feature)
check <- sapply(names(models), function(m) {
  k <- merge(maaslin_lm(abund, meta, models[[m]]), filter(res, model == m), by = c("feature", "name"))
  max(abs(k$qval.x - k$qval.y))
})
print(check)

# Sensitivity analyses (Supplementary Table S4)
no_scc <- rownames(meta)[meta$Dx != "SCC"]
neg52  <- rownames(meta)[meta$HPV52 == "N"]
all37  <- rownames(meta)
adj    <- c("Group", "Age", "HPV52")

sens <- list(
  "Minimum abundance 0.1%"             = list(all37,  "Group", 0.1),
  "Minimum abundance 0.1%, age + HPV52" = list(all37,  adj,     0.1),
  "SCC excluded"                       = list(no_scc, "Group", 0),
  "SCC excluded, age + HPV52"          = list(no_scc, adj,     0),
  "Adjusted for age and HPV52"         = list(all37,  adj,     0),
  "HPV52-negative only"                = list(neg52,  "Group", 0),
  "HPV52-negative only, age"           = list(neg52,  c("Group", "Age"), 0))

sens_res <- bind_rows(lapply(names(sens), function(s) {
  x <- sens[[s]]
  mutate(maaslin_lm(abund[x[[1]], ], meta[x[[1]], ], x[[2]], min_abund = x[[3]]), analysis = s)
}))
write.csv(sens_res, file.path(out, "sensitivity_analyses.csv"), row.names = FALSE)
filter(sens_res, feature == "Segatella copri", name == "Groupmultiple")
