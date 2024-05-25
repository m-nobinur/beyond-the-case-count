# 05_mh_ml.R — Supervised machine learning: predicting probable depression (PHQ-9 >= 10).
# Four learners x three nested feature tiers, tuned by cross-validation on a training split,
# evaluated once on an untouched stratified hold-out test set (TRIPOD-style reporting).
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(glmnet); library(ranger); library(xgboost); library(pROC); library(parallel) })

mh <- readRDS(P("data/processed/mh_clean.rds"))

# Explicit feature whitelist — never the PHQ items, totals or authors' derived columns (leakage).
TIER <- list(
  A = c("sex", "age", "education", "mh_history", "division", "residence", "living", "profession",
        "frontliner", "ses", "married", "children", "covid_self", "covid_family", "vaccine"),
  B = c("optimism", "pessimism", "mindfulness"),
  C = c("anxiety", "stress"))
TIER$B <- c(TIER$A, TIER$B); TIER$C <- c(TIER$B, TIER$C)
TIER_LAB <- c(A = "Tier A: demographics + COVID", B = "Tier B: + optimism, pessimism, mindfulness", C = "Tier C: + anxiety, stress")
stopifnot(!any(c("depression", "dep", "phq2") %in% unlist(TIER)))

# ---- Stratified hold-out split -------------------------------------------------
set.seed(SEED)
test_idx <- unlist(lapply(split(seq_len(nrow(mh)), mh$dep), function(i) sample(i, round(.25 * length(i)))))
train <- mh[-test_idx, ]; test <- mh[test_idx, ]
stopifnot(abs(mean(train$dep) - mean(test$dep)) < .01)

mm <- function(df, vars, ref) {   # design matrix with dummy levels fixed by the training data
  for (v in vars) if (is.factor(ref[[v]])) df[[v]] <- factor(df[[v]], levels(ref[[v]]))
  model.matrix(reformulate(vars), df)[, -1, drop = FALSE]
}
K_FOLD <- 10
strat_folds <- function(y, k) { f <- integer(length(y)); for (i in split(seq_along(y), y)) f[i] <- sample(rep_len(1:k, length(i))); f }
set.seed(SEED)
foldid <- strat_folds(train$dep, K_FOLD)
stopifnot(all(foldid %in% 1:K_FOLD))

auc_fast <- function(y, p) as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<")))

# ---- Learners: each returns list(predict = function(newdf) prob, cv_pred = OOF probs, info) ----
# `light = TRUE` gives a cheaper tuning grid, single-threaded, for the repeated-split evaluation.
fit_logit <- function(vars, y = "dep", tr = train, fid = foldid, light = FALSE) {
  f <- reformulate(vars, y)
  oof <- numeric(nrow(tr))
  if (!light) for (k in unique(fid)) {
    m <- glm(f, data = tr[fid != k, ], family = binomial)
    oof[fid == k] <- predict(m, tr[fid == k, ], type = "response")
  }
  m <- glm(f, data = tr, family = binomial)
  list(predict = function(nd) predict(m, nd, type = "response"), oof = oof, info = "unpenalised", model = m)
}
fit_enet <- function(vars, y = "dep", tr = train, fid = foldid, light = FALSE) {
  X <- mm(tr, vars, tr); alphas <- if (light) c(0, .5, 1) else c(0, .25, .5, .75, 1)
  # tuned on binomial deviance: AUC ignores the scale of the predictions and over-shrinks
  cvs <- map(alphas, ~ cv.glmnet(X, tr[[y]], family = "binomial", alpha = .x,
                                 foldid = fid, type.measure = "deviance", keep = TRUE))
  best <- which.min(map_dbl(cvs, ~ min(.x$cvm)))
  cv <- cvs[[best]]; a <- alphas[best]
  oof <- plogis(cv$fit.preval[, which(cv$lambda == cv$lambda.min)])
  list(predict = function(nd) as.numeric(predict(cv, mm(nd, vars, tr), s = "lambda.min", type = "response")),
       oof = oof, info = sprintf("alpha = %s, lambda = %s", a, signif(cv$lambda.min, 3)), model = cv, alpha = a)
}
fit_rf <- function(vars, y = "dep", tr = train, fid = foldid, light = FALSE) {
  d <- tr[, c(vars, y)] %>% mutate(!!y := factor(.data[[y]]))
  nt <- if (light) 1 else 4
  grid <- if (light) expand_grid(mtry = unique(pmin(c(2, 4), length(vars))), min.node.size = c(15, 30))
          else expand_grid(mtry = unique(pmin(c(2, 3, 5, 8), length(vars))), min.node.size = c(5, 15, 30))
  grid$auc <- pmap_dbl(grid, function(mtry, min.node.size) {
    m <- ranger(reformulate(vars, y), d, num.trees = if (light) 300 else 500, mtry = mtry, min.node.size = min.node.size,
                probability = TRUE, seed = SEED, num.threads = nt)
    auc_fast(tr[[y]], m$predictions[, "1"])           # out-of-bag AUC
  })
  g <- grid %>% slice_max(auc, n = 1, with_ties = FALSE)
  m <- ranger(reformulate(vars, y), d, num.trees = if (light) 500 else 1000, mtry = g$mtry, min.node.size = g$min.node.size,
              probability = TRUE, seed = SEED, num.threads = nt)
  list(predict = function(nd) predict(m, nd)$predictions[, "1"], oof = m$predictions[, "1"],
       info = sprintf("mtry = %d, min.node = %d", g$mtry, g$min.node.size), model = m)
}
fit_xgb <- function(vars, y = "dep", tr = train, fid = foldid, light = FALSE) {
  X <- mm(tr, vars, tr); dtr <- xgb.DMatrix(X, label = tr[[y]])
  folds <- map(sort(unique(fid)), ~ which(fid == .x))
  nt <- if (light) 1 else 4
  grid <- if (light) expand_grid(max_depth = c(2, 3), min_child_weight = 5) else expand_grid(max_depth = c(2, 3, 4), min_child_weight = c(1, 5))
  res <- pmap_dfr(grid, function(max_depth, min_child_weight) {
    set.seed(SEED)
    cv <- xgb.cv(params = list(objective = "binary:logistic", eval_metric = "auc", eta = if (light) .06 else .03,
                               max_depth = max_depth, min_child_weight = min_child_weight,
                               subsample = .8, colsample_bytree = .8, nthread = nt),
                 data = dtr, nrounds = 1500, folds = folds, early_stopping_rounds = 50, verbose = 0, prediction = !light)
    ev <- as.data.frame(cv$evaluation_log)
    tibble(max_depth, min_child_weight, nrounds = which.max(ev$test_auc_mean), auc = max(ev$test_auc_mean),
           oof = list(if (light) numeric(0) else as.numeric(cv$cv_predict$pred)))
  })
  g <- res %>% slice_max(auc, n = 1, with_ties = FALSE)
  set.seed(SEED)
  m <- xgb.train(params = list(objective = "binary:logistic", eta = if (light) .06 else .03, max_depth = g$max_depth,
                               min_child_weight = g$min_child_weight, subsample = .8, colsample_bytree = .8, nthread = nt),
                 data = dtr, nrounds = g$nrounds, verbose = 0)
  list(predict = function(nd) predict(m, xgb.DMatrix(mm(nd, vars, tr))), oof = g$oof[[1]],
       info = sprintf("depth = %d, mcw = %d, rounds = %d", g$max_depth, g$min_child_weight, g$nrounds), model = m)
}
LEARNERS <- list(`Logistic regression` = fit_logit, `Elastic net` = fit_enet, `Random forest` = fit_rf, XGBoost = fit_xgb)

# ---- Evaluation helpers --------------------------------------------------------
cal_fit <- function(y, p) {   # calibration intercept (in the large) and slope
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  slope <- coef(glm(y ~ lp, family = binomial))[2]
  int <- coef(glm(y ~ offset(lp), family = binomial))[1]
  c(int = unname(int), slope = unname(slope))
}
evaluate <- function(y, p, thr) {
  r <- roc(y, p, quiet = TRUE, direction = "<"); ci <- ci.auc(r, method = "delong")
  pred <- as.integer(p >= thr); cf <- cal_fit(y, p)
  tibble(AUC = as.numeric(auc(r)), AUC_lo = ci[1], AUC_hi = ci[3], Brier = mean((p - y)^2),
         Sens = mean(pred[y == 1] == 1), Spec = mean(pred[y == 0] == 0),
         Acc = mean(pred == y), Cal_int = cf["int"], Cal_slope = cf["slope"])
}

# ---- Fit all learners x tiers -------------------------------------------------
fits <- list(); results <- list(); test_preds <- list()
for (t in names(TIER)) for (l in names(LEARNERS)) {
  cat("fitting", l, "tier", t, "\n")
  f <- LEARNERS[[l]](TIER[[t]])
  # threshold chosen on training OOF predictions (Youden), never on the test set
  thr <- coords(roc(train$dep, f$oof, quiet = TRUE, direction = "<"), "best", best.method = "youden", ret = "threshold")$threshold[1]
  p_test <- f$predict(test)
  fits[[paste(t, l)]] <- f
  test_preds[[paste(t, l)]] <- tibble(tier = t, learner = l, y = test$dep, p = p_test)
  results[[paste(t, l)]] <- evaluate(test$dep, p_test, thr) %>%
    mutate(tier = t, learner = l, cv_auc = auc_fast(train$dep, f$oof), thr = thr, tuning = f$info)
}
res <- bind_rows(results)
preds <- bind_rows(test_preds)
write_csv(res, P("output", "ml_results_r.csv"))

# Paired DeLong tests: best learner vs. logistic baseline within each tier; tier B vs. tier A for each learner
delong <- function(a, b) roc.test(roc(test$dep, test_preds[[a]]$p, quiet = TRUE, direction = "<"),
                                  roc(test$dep, test_preds[[b]]$p, quiet = TRUE, direction = "<"), method = "delong")$p.value
# "Best" learner per tier is chosen on training CV AUC; the test set is never used for selection.
best_by_tier <- res %>% group_by(tier) %>% slice_max(cv_auc, n = 1, with_ties = FALSE) %>% ungroup()
best_vs_logit <- best_by_tier %>% mutate(p_vs_logit = map2_dbl(tier, learner, ~ if (.y == "Logistic regression") NA
                                                               else delong(paste(.x, .y), paste(.x, "Logistic regression"))))
tierB_vs_A <- map_dbl(names(LEARNERS), ~ delong(paste("B", .x), paste("A", .x))) %>% set_names(names(LEARNERS))

# Sensitivity: tier-B test AUC after dropping the templated records (flagged in 03) from the test set
clean_auc <- map_dbl(names(LEARNERS), function(l) auc_fast(test$dep[!test$suspect], test_preds[[paste("B", l)]]$p[!test$suspect])) %>%
  set_names(names(LEARNERS))

# ---- Repeated stratified hold-out: is the learner ranking stable across splits? ----
# One 75/25 split gives test AUCs with SD about 0.02, as large as the gaps between learners. Each
# repetition re-splits, re-tunes (lighter grids) and re-tests all 12 models.
N_REP <- 50
one_split <- function(r) {
  set.seed(SEED + 1000 + r)
  ti <- unlist(lapply(split(seq_len(nrow(mh)), mh$dep), function(i) sample(i, round(.25 * length(i)))))
  tr <- mh[-ti, ]; te <- mh[ti, ]; fid <- strat_folds(tr$dep, 5)
  out <- list()
  for (t in names(TIER)) for (l in names(LEARNERS)) {
    f <- LEARNERS[[l]](TIER[[t]], tr = tr, fid = fid, light = TRUE)
    out[[length(out) + 1]] <- tibble(rep = r, tier = t, learner = l, auc = auc_fast(te$dep, f$predict(te)))
  }
  # feature ablation with logistic regression on the same split
  for (a in names(ABL)) {
    v <- ABL[[a]]; m <- glm(reformulate(v, "dep"), data = tr, family = binomial)
    out[[length(out) + 1]] <- tibble(rep = r, tier = "ablation", learner = a, auc = auc_fast(te$dep, predict(m, te)))
  }
  bind_rows(out)
}
ABL <- list(`Tier A (15 variables)` = TIER$A,
            `Tier A without division` = setdiff(TIER$A, "division"),
            `Tier A + optimism` = c(TIER$A, "optimism"),
            `Tier A + pessimism` = c(TIER$A, "pessimism"),
            `Tier A + mindfulness` = c(TIER$A, "mindfulness"),
            `Tier B (A + all three scales)` = TIER$B,
            `Mindfulness alone` = "mindfulness",
            `Optimism alone` = "optimism",
            `Reference: PHQ-2 score alone` = "phq2")
stopifnot(!"phq2" %in% unlist(TIER))
cat("repeated splits:", N_REP, "x 12 models\n")
rep_res <- bind_rows(mclapply(seq_len(N_REP), one_split, mc.cores = min(8, detectCores()), mc.set.seed = FALSE))
stopifnot(nrow(rep_res) == N_REP * (12 + length(ABL)))
rep_ml <- rep_res %>% filter(tier != "ablation")
rep_vs_logit <- rep_ml %>% group_by(rep, tier) %>% mutate(d = auc - auc[learner == "Logistic regression"]) %>% ungroup()
rep_sum <- rep_vs_logit %>% group_by(tier, learner) %>%
  summarise(mean = mean(auc), sd = sd(auc), lo = quantile(auc, .025), hi = quantile(auc, .975),
            d_lo = quantile(d, .025), d_hi = quantile(d, .975), win = mean(d > 0), d = mean(d), .groups = "drop")
abl_sum <- rep_res %>% filter(tier == "ablation") %>% group_by(learner) %>%
  summarise(mean = mean(auc), sd = sd(auc), lo = quantile(auc, .025), hi = quantile(auc, .975), .groups = "drop") %>%
  mutate(learner = factor(learner, names(ABL))) %>% arrange(learner)
save_tab(rep_sum %>% mutate(win = ifelse(learner == "Logistic regression", "--", pct(win, 0)),
                            dtxt = ifelse(learner == "Logistic regression", "--", sprintf("%+.3f [%+.3f, %+.3f]", d, d_lo, d_hi))) %>%
           transmute(Tier = TIER_LAB[tier], Learner = learner, `Mean AUC (SD)` = sprintf("%s (%s)", f3(mean), f3(sd)),
                     `95% range` = sprintf("%s–%s", f3(lo), f3(hi)), `Diff. vs. logistic [95% range]` = dtxt,
                     `Splits beating logistic (%)` = win),
         "tab_ml_repeated", font = "footnotesize", escape = FALSE,
         caption = sprintf("Test AUC over %d repeated stratified 75/25 splits. Every split re-tunes and re-tests all twelve models.", N_REP),
         notes = "Differences are paired within split. The 95\\\\% range is the 2.5th--97.5th percentile across splits. Tuning grids are smaller than in the fixed-split analysis.")
save_tab(abl_sum %>% transmute(`Feature set (logistic regression)` = as.character(learner),
                               `Mean AUC (SD)` = sprintf("%s (%s)", f3(mean), f3(sd)), `95% range` = sprintf("%s–%s", f3(lo), f3(hi))),
         "tab_ml_ablation", font = "footnotesize", escape = FALSE,
         caption = sprintf("Feature ablation: test AUC of logistic regression over the same %d splits.", N_REP),
         notes = "The PHQ-2 row uses two of the nine PHQ-9 items, so it overlaps with the outcome by construction. It is shown only as a ceiling for comparison and is never a model feature.")
LCOL0 <- setNames(unname(PAL[c("grey", "blue", "green", "orange")]), names(LEARNERS))
fig_rep <- rep_ml %>% mutate(tier = factor(TIER_LAB[tier], TIER_LAB), learner = factor(learner, rev(names(LEARNERS)))) %>%
  ggplot(aes(auc, learner, colour = learner)) +
  geom_jitter(height = .15, width = 0, alpha = .35, size = .9) +
  stat_summary(fun = mean, fun.min = function(x) quantile(x, .025), fun.max = function(x) quantile(x, .975),
               colour = "black", size = .3, linewidth = .5) +
  facet_wrap(~tier, nrow = 1, scales = "free_x") + scale_colour_manual(values = LCOL0, guide = "none") +
  labs(x = sprintf("Test AUC across %d random splits (dots); mean and 95%% range (black)", N_REP), y = NULL,
       title = "Model complexity does not buy discrimination")
save_fig(fig_rep, "fig_ml_repeated", width = 7, height = 2.9)

# ---- Results table -----------------------------------------------------------
tab_ml <- res %>% mutate(Tier = tier, Learner = learner,
                         `CV AUC` = f3(cv_auc),
                         `Test AUC [95% CI]` = sprintf("%s [%s, %s]", f3(AUC), f3(AUC_lo), f3(AUC_hi)),
                         Brier = f3(Brier), Sens. = f2(Sens), Spec. = f2(Spec),
                         `Cal. slope` = f2(Cal_slope)) %>%
  select(Tier, Learner, `CV AUC`, `Test AUC [95% CI]`, Brier, Sens., Spec., `Cal. slope`)
save_tab(tab_ml, "tab_ml_results", font = "footnotesize",
         caption = sprintf("Discrimination and calibration of four learners across three nested feature tiers for probable depression (training n = %d, test n = %d).", nrow(train), nrow(test)),
         notes = "CV AUC = 10-fold out-of-fold AUC on the training set (out-of-bag AUC for the random forest); test metrics computed once on the untouched hold-out set; 95% CI by DeLong's method. Sensitivity/specificity at the Youden-optimal threshold estimated from training out-of-fold predictions. Tier A: sociodemographic and COVID-19 exposure; B: + optimism, pessimism, mindfulness; C: + anxiety and stress scores.")
save_tab(res %>% distinct(tier, learner, tuning) %>% rename(Tier = tier, Learner = learner, `Selected hyper-parameters` = tuning),
         "tab_ml_tuning", font = "footnotesize",
         caption = "Hyper-parameters selected by cross-validation on the training set.")

# ---- Figures: ROC curves, AUC comparison, calibration ---------------------------
roc_df <- preds %>% group_by(tier, learner) %>% group_modify(function(d, k) {
  r <- roc(d$y, d$p, quiet = TRUE, direction = "<"); tibble(fpr = 1 - r$specificities, tpr = r$sensitivities)
}) %>% ungroup() %>% mutate(tier = TIER_LAB[tier]) %>% arrange(tier, learner, fpr, tpr)
LCOL <- setNames(unname(PAL[c("grey", "blue", "green", "orange")]), names(LEARNERS))
fig_roc <- ggplot(roc_df, aes(fpr, tpr, colour = learner)) + geom_abline(linetype = 3, colour = "grey60") +
  geom_step(linewidth = .55) + facet_wrap(~tier, ncol = 3) + coord_equal() +
  scale_colour_manual(values = LCOL, name = NULL) +
  scale_x_continuous(labels = percent, breaks = c(0, .5, 1)) + scale_y_continuous(labels = percent, breaks = c(0, .5, 1)) +
  labs(x = "False-positive rate (1 − specificity)", y = "True-positive rate", title = "ROC curves on the hold-out test set")
save_fig(fig_roc, "fig_ml_roc", width = 7, height = 3.4)

fig_auc <- res %>% mutate(tier = factor(TIER_LAB[tier], TIER_LAB)) %>%
  ggplot(aes(AUC, fct_rev(factor(learner, names(LEARNERS))), colour = tier)) +
  geom_pointrange(aes(xmin = AUC_lo, xmax = AUC_hi), position = position_dodge(.5), size = .3) +
  geom_point(aes(x = cv_auc), shape = 4, position = position_dodge(.5), size = 2) +
  scale_colour_manual(values = unname(PAL[c("sky", "blue", "purple")]), name = NULL) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = "AUC (point & 95% CI: test set; ×: training CV)", y = NULL, title = "Model comparison by feature tier")
save_fig(fig_auc, "fig_ml_auc", width = 6.5, height = 3.4)

best_B <- best_by_tier %>% filter(tier == "B")
cal_df <- preds %>% filter(tier == "B") %>% group_by(learner) %>%
  mutate(bin = ntile(p, 10)) %>% group_by(learner, bin) %>%
  summarise(pred = mean(p), obs = mean(y), n = n(), .groups = "drop") %>%
  mutate(se = sqrt(obs * (1 - obs) / n))
fig_cal <- ggplot(cal_df, aes(pred, obs, colour = learner)) + geom_abline(linetype = 2, colour = "grey50") +
  geom_line() + geom_point(size = 1.4) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = LCOL, name = NULL) + facet_wrap(~learner, nrow = 1) +
  scale_x_continuous(labels = c("0", "50", "100"), breaks = c(0, .5, 1)) + scale_y_continuous(labels = c("0", "50", "100"), breaks = c(0, .5, 1)) +
  labs(x = "Mean predicted risk, % (decile)", y = "Observed depressed, %", title = "Calibration on the test set (tier B)") +
  theme(legend.position = "none")
save_fig(fig_cal, "fig_ml_calibration", width = 7, height = 2.8)

# ---- Interpretation: grouped permutation importance on the test set -------------
perm_imp <- function(f, vars, n_rep = 30) {
  base <- auc_fast(test$dep, f$predict(test))
  map_dfr(vars, function(v) {
    drops <- map_dbl(1:n_rep, function(r) { set.seed(SEED + r); nd <- test; nd[[v]] <- sample(nd[[v]]); base - auc_fast(test$dep, f$predict(nd)) })
    tibble(var = v, imp = mean(drops), lo = quantile(drops, .05), hi = quantile(drops, .95))
  })
}
VLAB <- c(sex = "Sex", age = "Age", education = "Education", mh_history = "Mental-health history", division = "Division",
          residence = "Residence", living = "Living arrangement", profession = "Profession", frontliner = "Frontline worker",
          ses = "SES ladder", married = "Marital status", children = "Children", covid_self = "Own COVID-19 infection",
          covid_family = "Family COVID-19 experience", vaccine = "Vaccination status", optimism = "Optimism",
          pessimism = "Pessimism", mindfulness = "Mindfulness", anxiety = "Anxiety (BAI)", stress = "Stress (PSS-4)")
# Interpreted model: logistic regression, since no learner outperformed it (repeated-split analysis above).
imp_B <- perm_imp(fits[["B Logistic regression"]], TIER$B) %>% mutate(model = "Tier B, logistic regression")
imp_A <- perm_imp(fits[["A Logistic regression"]], TIER$A) %>% mutate(model = "Tier A, logistic regression")
imp <- bind_rows(imp_A, imp_B) %>% mutate(var = VLAB[var])
fig_imp <- imp %>% group_by(model) %>% slice_max(imp, n = 12) %>% ungroup() %>%
  mutate(var = reorder(paste(var, model, sep = "___"), imp)) %>%
  ggplot(aes(imp, var)) + geom_vline(xintercept = 0, colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, colour = "grey50", orientation = "y") + geom_point(colour = PAL["blue"], size = 2) +
  facet_wrap(~model, scales = "free") + scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  labs(x = "Drop in test AUC when the variable is permuted (mean, 5–95% range of 30 permutations)", y = NULL,
       title = "Permutation importance (top 12 variables)")
save_fig(fig_imp, "fig_ml_importance", width = 7, height = 3.9)

# ---- Secondary outcome: anxiety (BAI >= 16), tier A & B, all learners ---------
anx_res <- map_dfr(c("A", "B"), function(t) map_dfr(names(LEARNERS), function(l) {
  vars <- setdiff(TIER[[t]], "anxiety")
  f <- LEARNERS[[l]](vars, y = "anx")
  r <- roc(test$anx, f$predict(test), quiet = TRUE, direction = "<"); ci <- ci.auc(r)
  tibble(tier = t, learner = l, AUC = as.numeric(auc(r)), lo = ci[1], hi = ci[3])
}))
save_tab(anx_res %>% transmute(Tier = TIER_LAB[tier], Learner = learner, `Test AUC [95% CI]` = sprintf("%s [%s, %s]", f3(AUC), f3(lo), f3(hi))),
         "tab_ml_anxiety", font = "footnotesize",
         caption = "Secondary outcome: hold-out discrimination for probable anxiety (BAI $\\geq$ 16).", escape = FALSE)

# ---- Export tier-B design matrices for the Python cross-check ------------------
export <- function(df, file) {
  X <- as.data.frame(mm(df, TIER$B, train)); names(X) <- make.names(names(X)); X$dep <- df$dep
  write_csv(X, P("data/processed", file))
}
export(train, "ml_train_tierB.csv"); export(test, "ml_test_tierB.csv")

# ---- Macros ------------------------------------------------------------------
macro("MlTrainN", fint(nrow(train))); macro("MlTestN", fint(nrow(test)))
g <- function(t, l, col, fn = f3) fn(res[[col]][res$tier == t & res$learner == l])
for (t in names(TIER)) for (l in names(LEARNERS)) {
  key <- paste0("Auc", t, c(`Logistic regression` = "Logit", `Elastic net` = "Enet", `Random forest` = "Rf", XGBoost = "Xgb")[[l]])
  macro(key, g(t, l, "AUC"))
}
for (t in names(TIER)) {
  b <- best_by_tier %>% filter(tier == t)
  macro(paste0("Best", t, "Learner"), b$learner); macro(paste0("Best", t, "Auc"), f3(b$AUC))
  macro(paste0("Best", t, "Lo"), f3(b$AUC_lo)); macro(paste0("Best", t, "Hi"), f3(b$AUC_hi))
  macro(paste0("Best", t, "Brier"), f3(b$Brier)); macro(paste0("Best", t, "Sens"), pct(b$Sens, 0)); macro(paste0("Best", t, "Spec"), pct(b$Spec, 0))
  macro(paste0("Best", t, "Slope"), f2(b$Cal_slope))
  pv <- best_vs_logit$p_vs_logit[best_vs_logit$tier == t]; macro(paste0("Best", t, "VsLogitP"), ifelse(is.na(pv), "n/a", fp(pv)))
}
macro("TierBvsALogitP", fp(tierB_vs_A["Logistic regression"])); macro("TierBvsAXgbP", fp(tierB_vs_A["XGBoost"]))
macro("TierBvsARfP", fp(tierB_vs_A["Random forest"])); macro("TierBvsAEnetP", fp(tierB_vs_A["Elastic net"]))
macro("AucRangeLo", f3(min(res$AUC))); macro("AucRangeHi", f3(max(res$AUC)))
macro("AucSpreadB", f3(diff(range(res$AUC[res$tier == "B"]))))
macro("ImpTopB", imp_B %>% slice_max(imp, n = 1) %>% pull(var) %>% VLAB[.])
macro("ImpTopBval", f3(max(imp_B$imp)))
macro("ImpSecondB", imp_B %>% arrange(desc(imp)) %>% slice(2) %>% pull(var) %>% VLAB[.])
macro("ImpTopA", imp_A %>% slice_max(imp, n = 1) %>% pull(var) %>% VLAB[.])
macro("ImpSecondA", imp_A %>% arrange(desc(imp)) %>% slice(2) %>% pull(var) %>% VLAB[.])
macro("ImpThirdA", imp_A %>% arrange(desc(imp)) %>% slice(3) %>% pull(var) %>% VLAB[.])
ab <- anx_res %>% group_by(tier) %>% slice_max(AUC, n = 1)
macro("AnxBestA", f3(ab$AUC[ab$tier == "A"])); macro("AnxBestB", f3(ab$AUC[ab$tier == "B"]))
enB <- fits[["B Elastic net"]]; macro("EnetAlphaB", enB$alpha)
nz <- coef(enB$model, s = "lambda.min"); macro("EnetNzB", sum(nz[-1] != 0)); macro("EnetPB", length(nz) - 1)
gapB <- res %>% filter(tier == "B") %>% mutate(g = cv_auc - AUC)
macro("CvGapLo", f2(min(gapB$g))); macro("CvGapHi", f2(max(gapB$g)))
cal_best <- res %>% group_by(tier) %>% slice_min(abs(Cal_slope - 1), n = 1) %>% ungroup()
for (t in c("A", "B", "C")) { l <- cal_best$learner[cal_best$tier == t]; macro(paste0("CalBest", t), ifelse(l == "XGBoost", l, tolower(l))) }
macro("SensMlLo", f3(min(clean_auc))); macro("SensMlHi", f3(max(clean_auc))); macro("SensMlNtest", sum(!test$suspect))
LK <- c(`Logistic regression` = "Logit", `Elastic net` = "Enet", `Random forest` = "Rf", XGBoost = "Xgb")
for (i in seq_len(nrow(rep_sum))) {
  k <- paste0(rep_sum$tier[i], LK[[rep_sum$learner[i]]])
  macro(paste0("RepAuc", k), f3(rep_sum$mean[i])); macro(paste0("RepSd", k), f3(rep_sum$sd[i]))
  macro(paste0("RepWin", k), pct(rep_sum$win[i], 0)); macro(paste0("RepDiff", k), sprintf("%+.3f", rep_sum$d[i]))
}
for (i in seq_len(nrow(res))) macro(paste0("Slope", res$tier[i], LK[[res$learner[i]]]), f2(res$Cal_slope[i]))
macro("TestCiWidthLo", f2(min(res$AUC_hi - res$AUC_lo))); macro("TestCiWidthHi", f2(max(res$AUC_hi - res$AUC_lo)))
lb <- res %>% filter(tier == "B", learner == "Logistic regression"); macro("LogitSensB", f2(lb$Sens)); macro("LogitSpecB", f2(lb$Spec))
for (t in c("A", "B")) { a <- anx_res$AUC[anx_res$tier == t]; macro(paste0("Anx", t, "Lo"), f3(min(a))); macro(paste0("Anx", t, "Hi"), f3(max(a))) }
ce <- rep_sum %>% filter(tier == "C", learner == "Elastic net")
macro("RepDiffCEnetLo", sprintf("%+.3f", ce$d_lo)); macro("RepDiffCEnetHi", sprintf("%+.3f", ce$d_hi))
macro("RepN", N_REP); macro("RepSdMin", f3(min(rep_sum$sd))); macro("RepSdMax", f3(max(rep_sum$sd)))
macro("RepMaxAbsDiff", f3(max(abs(rep_sum$d))))
for (t in names(TIER)) { b <- rep_sum %>% filter(tier == t) %>% slice_max(mean, n = 1)
  macro(paste0("RepBest", t, "Learner"), ifelse(b$learner == "XGBoost", b$learner, tolower(b$learner))); macro(paste0("RepBest", t, "Auc"), f3(b$mean)) }
ag <- function(nm) f3(abl_sum$mean[abl_sum$learner == nm])
macro("AblA", ag("Tier A (15 variables)")); macro("AblANoDiv", ag("Tier A without division"))
macro("AblAOpt", ag("Tier A + optimism")); macro("AblAPess", ag("Tier A + pessimism")); macro("AblAMind", ag("Tier A + mindfulness"))
macro("AblB", ag("Tier B (A + all three scales)")); macro("AblMind", ag("Mindfulness alone")); macro("AblOpt", ag("Optimism alone"))
macro("AblPhqTwo", ag("Reference: PHQ-2 score alone"))
write_macros("mh_ml")
saveRDS(list(rep = rep_sum, abl = abl_sum, res = res, best = best_vs_logit, tierB_vs_A = tierB_vs_A, imp = imp, anx = anx_res), P("data/processed/ml_results.rds"))
print(res %>% select(tier, learner, cv_auc, AUC, AUC_lo, AUC_hi, Brier, Cal_slope, tuning), width = 200)
print(best_vs_logit %>% select(tier, learner, AUC, p_vs_logit)); print(tierB_vs_A); print(rep_sum, width = 200); print(abl_sum)
