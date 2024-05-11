# 04_mh_profiles.R — Unsupervised learning: latent mental-health profiles (Gaussian mixture models),
# cross-checked with k-means, then profiled against sociodemographic and COVID-19 exposure variables.
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(mclust); library(cluster); library(nnet) })

mh_all <- readRDS(P("data/processed/mh_clean.rds"))
# Templated records (identical response blocks, flagged in 03) form a zero-variance spike that a
# mixture model captures as its own degenerate class; profiles are therefore estimated without them.
mh <- filter(mh_all, !suspect)
# Pessimism (alpha = 0.35) is excluded: an unreliable indicator would add noise to the clustering.
ind <- c(Anxiety = "anxiety", Depression = "depression", Stress = "stress", Optimism = "optimism", Mindfulness = "mindfulness")
Z <- scale(mh[, ind]); colnames(Z) <- names(ind)

# ---- GMM model selection by BIC ------------------------------------------------
set.seed(SEED)
# Diagonal (conditional-independence) parameterisations = classical latent profile analysis.
LPA <- c("EEI", "VEI", "EVI", "VVI")
bic <- mclustBIC(Z, G = 1:8, modelNames = LPA)
gmm_bic <- Mclust(Z, x = bic)
bic_df <- as.data.frame(unclass(bic)) %>% mutate(G = 1:8) %>% pivot_longer(-G, names_to = "model", values_to = "BIC") %>% filter(!is.na(BIC))
top_models <- bic_df %>% group_by(model) %>% summarise(m = max(BIC)) %>% slice_max(m, n = 5) %>% pull(model)
fig_bic <- bic_df %>% filter(model %in% top_models) %>%
  ggplot(aes(G, BIC, colour = model)) + geom_line() + geom_point(size = 1.2) +
  geom_vline(xintercept = gmm_bic$G, linetype = 2, colour = "grey40") +
  scale_x_continuous(breaks = 1:8) + scale_colour_brewer(palette = "Dark2", name = "Covariance model") +
  labs(x = "Number of mixture components (G)", y = "BIC (higher is better)", title = "Gaussian mixture model selection")
save_fig(fig_bic, "fig_profiles_bic", width = 6, height = 3.4)

# Integrated Completed Likelihood as a second criterion that penalises overlapping classes
icl <- mclustICL(Z, G = 1:8, modelNames = LPA)
icl_best <- summary(icl)[1] %>% names()
best_family <- strsplit(icl_best, ",")[[1]][1]
G_sel <- as.integer(strsplit(icl_best, ",")[[1]][2])
gmm <- Mclust(Z, G = G_sel, modelNames = best_family)   # primary solution: ICL-optimal

# Full-sample check: with the templated records included, does a degenerate class appear?
Z_all <- scale(mh_all[, ind])
icl_all <- mclustICL(Z_all, G = 1:8, modelNames = LPA); best_all <- strsplit(names(summary(icl_all))[1], ",")[[1]]
gmm_all <- Mclust(Z_all, G = 4, modelNames = "VEI")
spike <- which.min(sapply(1:4, function(g) mean(apply(Z_all[gmm_all$classification == g, , drop = FALSE], 2, sd))))
spike_overlap <- mean(mh_all$suspect[gmm_all$classification == spike])
# How general is this? Across all diagonal models with G = 3..6, how many full-sample fits contain a
# small class (< 5% of respondents) made up mostly (>= 80%) of flagged records?
spike_grid <- expand_grid(model = LPA, G = 3:6) %>%
  mutate(has_spike = map2_lgl(model, G, function(m, g) {
    fit <- Mclust(Z_all, G = g, modelNames = m, verbose = FALSE)
    if (is.null(fit)) return(NA)
    any(sapply(seq_len(fit$G), function(k) { i <- fit$classification == k; mean(i) < .05 && mean(mh_all$suspect[i]) >= .8 }))
  }))

# ---- Order and name profiles by mean depression ----------------------------------
cl <- gmm$classification
ord <- order(tapply(Z[, "Depression"], cl, mean))
mh$profile_id <- match(cl, ord)
K <- gmm$G
means <- as_tibble(Z) %>% mutate(profile = mh$profile_id) %>% group_by(profile) %>% summarise(across(everything(), mean), n = n())
name_profile <- function(r) {
  if (r$Depression < -0.5) "Low band"
  else if (r$Depression < 0.5) "Middle band"
  else "High band"
}
means$label <- map_chr(seq_len(nrow(means)), ~ name_profile(means[.x, ]))
means$label <- make.unique(means$label, sep = " ")
means$label <- sprintf("P%d: %s", means$profile, means$label)
mh$profile <- factor(means$label[mh$profile_id], means$label)
ref_prof <- means$label[grepl("Low band", means$label)][1]
avepp <- tapply(apply(gmm$z, 1, max), cl, mean)[ord]            # average posterior probability per class
entropy <- 1 + sum(gmm$z * log(pmax(gmm$z, 1e-12))) / (nrow(Z) * log(K))   # relative entropy (1 = perfect separation)

fig_prof <- means %>% pivot_longer(all_of(names(ind)), names_to = "measure", values_to = "z") %>%
  mutate(measure = factor(measure, names(ind)), lab = sprintf("%s (n = %d)", label, n)) %>%
  ggplot(aes(measure, z, group = lab, colour = lab)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_line(linewidth = .9) + geom_point(size = 2) +
  scale_colour_manual(values = unname(PAL[c("green", "sky", "orange", "red", "purple", "blue")])[1:K], name = NULL) +
  guides(colour = guide_legend(ncol = 2)) +
  labs(x = NULL, y = "Mean standardised score (z)", title = sprintf("Latent mental-health profiles (GMM, %s, G = %d)", gmm$modelName, K))
save_fig(fig_prof, "fig_profiles_means", width = 6.3, height = 3.8)

# ---- Cross-check with k-means (same K) -------------------------------------------
set.seed(SEED)
km <- kmeans(Z, centers = K, nstart = 50, iter.max = 100)
ari <- adjustedRandIndex(km$cluster, cl)
sil_k <- map_dbl(2:8, function(k) { set.seed(SEED); kk <- kmeans(Z, k, nstart = 25); mean(silhouette(kk$cluster, dist(Z))[, 3]) })
sil_best <- (2:8)[which.max(sil_k)]
sil_gmm <- mean(silhouette(cl, dist(Z))[, 3])
# Profiles are partly defined by PHQ-9; refit without the depression indicator to see how much they depend on it.
gmm_nodep <- Mclust(Z[, setdiff(colnames(Z), "Depression")], G = K, modelNames = gmm$modelName, verbose = FALSE)
ari_nodep <- adjustedRandIndex(gmm_nodep$classification, cl)
# Dimensional check: how well does a single severity score (first principal component) order the profiles?
pc1 <- prcomp(Z)$x[, 1]; pc1_pve <- summary(prcomp(Z))$importance[2, 1]
r2_pc1 <- summary(lm(pc1 ~ factor(cl)))$r.squared

# ---- Profile-level descriptives -------------------------------------------------
desc <- mh %>% group_by(profile) %>%
  summarise(n = n(), `%` = pct(n() / nrow(mh)),
            `PHQ-9` = sprintf("%s (%s)", f1(mean(depression)), f1(sd(depression))),
            `BAI` = sprintf("%s (%s)", f1(mean(anxiety)), f1(sd(anxiety))),
            `PSS-4` = sprintf("%s (%s)", f1(mean(stress)), f1(sd(stress))),
            `Optimism` = sprintf("%s (%s)", f1(mean(optimism)), f1(sd(optimism))),
            `MAAS` = sprintf("%s (%s)", f1(mean(mindfulness)), f1(sd(mindfulness))),
            `Dep. %` = pct(mean(dep)), `Anx. %` = pct(mean(anx)))
avepp_v <- f2(avepp)
desc <- desc %>% mutate(AvePP = avepp_v) %>% rename(Profile = profile)
save_tab(desc, "tab_profiles", digits = 1, font = "footnotesize",
         caption = "Characteristics of the latent mental-health profiles: mean (SD) of each scale, share meeting clinical thresholds, and average posterior membership probability (AvePP).",
         notes = sprintf("Model %s with G = %d selected by ICL (N = %s after excluding templated records). Relative entropy = %s. AvePP = average posterior probability of membership in the assigned profile.", gmm$modelName, K, fint(nrow(mh)), f2(entropy)))

# ---- Who belongs to the distressed profiles? Multinomial logistic regression ------
mn <- multinom(relevel(profile, ref_prof) ~ sex + age + mh_history + profession + residence + ses + covid_self + covid_family + vaccine,
               data = mh, trace = FALSE)
s <- summary(mn); zst <- s$coefficients / s$standard.errors
rrr <- map_dfr(rownames(s$coefficients), function(r) tibble(
  profile = r, term = colnames(s$coefficients), RRR = exp(s$coefficients[r, ]),
  lo = exp(s$coefficients[r, ] - 1.96 * s$standard.errors[r, ]), hi = exp(s$coefficients[r, ] + 1.96 * s$standard.errors[r, ]),
  p = 2 * pnorm(-abs(zst[r, ])))) %>% filter(term != "(Intercept)") %>% mutate(q = p.adjust(p, "BH"))
lab <- function(t) t %>% str_replace("^sex", "Sex: ") %>% str_replace("^mh_history", "MH history: ") %>%
  str_replace("^profession", "Profession: ") %>% str_replace("^residence", "Residence: ") %>%
  str_replace("^covid_self", "Own COVID: ") %>% str_replace("^covid_family", "Family COVID: ") %>%
  str_replace("^vaccine", "Vaccine: ") %>% recode(age = "Age (per year)", ses = "SES ladder (per step)")
fig_rrr <- rrr %>% mutate(term = lab(term)) %>%
  ggplot(aes(RRR, fct_rev(factor(term)), colour = profile, alpha = q < .05)) +
  geom_vline(xintercept = 1, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(.6), size = .2) +
  scale_x_log10() + scale_alpha_manual(values = c(`FALSE` = .35, `TRUE` = 1), guide = "none") +
  scale_colour_manual(values = unname(PAL[c("sky", "orange", "red", "purple", "blue")])[1:(K - 1)], name = NULL) +
  labs(x = sprintf("Relative risk ratio vs. %s (log scale, 95%% CI)", ref_prof), y = NULL,
       title = "Predictors of profile membership (multinomial logistic regression)", caption = "Faded: BH-adjusted q ≥ 0.05.")
save_fig(fig_rrr, "fig_profiles_rrr", width = 6.5, height = 5)

# ---- Macros --------------------------------------------------------------------
macro("GmmBicK", gmm_bic$G); macro("GmmBicModel", gmm_bic$modelName); macro("GmmK", K); macro("GmmModel", gmm$modelName); macro("GmmBIC", f1(gmm$bic)); macro("GmmEntropy", f2(entropy))
macro("GmmICLbest", icl_best); macro("GmmFamily", best_family); macro("GmmN", fint(nrow(mh)));
macro("GmmAllICL", paste(best_all, collapse = ",")); macro("SpikeN", sum(gmm_all$classification == spike)); macro("SpikeOverlap", pct(spike_overlap, 0)); macro("GmmMinAvePP", f2(min(avepp))); macro("GmmSil", f2(sil_gmm))
macro("KmARI", f2(ari)); macro("KmSilBest", sil_best); macro("KmSilMax", f2(max(sil_k)))
macro("SpikeFits", sum(spike_grid$has_spike, na.rm = TRUE)); macro("SpikeFitsTotal", sum(!is.na(spike_grid$has_spike)))
macro("AriNoDep", f2(ari_nodep)); macro("ProfPcOnePve", pct(pc1_pve)); macro("ProfPcOneRsq", f2(r2_pc1))
for (i in seq_len(K)) {
  nm <- c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight")[i]
  macro(paste0("Prof", nm, "Label"), means$label[i]); macro(paste0("Prof", nm, "Pct"), pct(means$n[i] / nrow(mh)))
  macro(paste0("Prof", nm, "Dep"), pct(mean(mh$dep[mh$profile_id == i])))
  macro(paste0("Prof", nm, "Phq"), f1(mean(mh$depression[mh$profile_id == i]))); macro(paste0("Prof", nm, "Bai"), f1(mean(mh$anxiety[mh$profile_id == i])))
}
top <- levels(mh$profile)[K]
r_top <- rrr %>% filter(profile == top)
macro("ProfTopPct", pct(mean(mh$profile == top))); macro("ProfTopN", sum(mh$profile == top))
gr <- function(t) r_top %>% filter(term == t)
macro("ProfTopFamDeath", f2(gr("covid_familyInfected, a death")$RRR)); macro("ProfTopFamDeathq", fp(gr("covid_familyInfected, a death")$q))
macro("ProfTopHist", f2(gr("mh_historyYes")$RRR)); macro("ProfTopHistq", fp(gr("mh_historyYes")$q))
macro("ProfTopFemale", f2(gr("sexFemale")$RRR)); macro("ProfTopFemaleq", fp(gr("sexFemale")$q))
macro("ProfTopTwoDose", f2(gr("vaccineTwo doses")$RRR)); macro("ProfTopTwoDoseq", fp(gr("vaccineTwo doses")$q))
macro("ProfTopSES", f2(gr("ses")$RRR)); macro("ProfTopSESq", fp(gr("ses")$q))
macro("ProfTopStudent", f2(gr("professionStudent")$RRR)); macro("ProfTopStudentq", fp(gr("professionStudent")$q))
write_macros("mh_profiles")
saveRDS(mh %>% select(id, profile), P("data/processed/mh_profiles.rds"))
cat("04 done: G =", K, gmm$modelName, " ARI =", ari, " entropy =", entropy, " ICL best:", icl_best, "\n")
print(means); print(rrr %>% filter(q < .05))
