# 06_ml_crosscheck.R — compare R (tier B) with the independent scikit-learn re-implementation.
source(here::here("scripts", "00_setup.R"))
r <- read_csv(P("output/ml_results_r.csv"), show_col_types = FALSE) %>% filter(tier == "B")
py <- read_csv(P("output/ml_results_python.csv"), show_col_types = FALSE)
pyi <- read_csv(P("output/ml_importance_python.csv"), show_col_types = FALSE)
pair <- tribble(~family, ~r_learner, ~py_learner,
  "Linear (logistic)", "Elastic net", "Logistic regression (L2)",
  "Random forest", "Random forest", "Random forest",
  "Gradient boosting", "XGBoost", "Hist. gradient boosting")
cmp <- pair %>% left_join(r %>% select(r_learner = learner, r_cv = cv_auc, r_auc = AUC, r_lo = AUC_lo, r_hi = AUC_hi), by = "r_learner") %>%
  left_join(py %>% select(py_learner = learner, py_cv = cv_auc, py_auc = test_auc, py_lo = lo, py_hi = hi), by = "py_learner")
save_tab(cmp %>% transmute(`Model family` = family, `R implementation` = r_learner,
                           `R test AUC [95% CI]` = sprintf("%s [%s, %s]", f3(r_auc), f3(r_lo), f3(r_hi)),
                           `Python implementation` = py_learner,
                           `Python test AUC [95% CI]` = sprintf("%s [%s, %s]", f3(py_auc), f3(py_lo), f3(py_hi)),
                           `$\\Delta$AUC` = f3(py_auc - r_auc)),
         "tab_ml_crosscheck", font = "footnotesize", escape = FALSE,
         caption = "Cross-implementation check for tier B: R (glmnet, ranger, xgboost) versus Python (scikit-learn) on the identical train/test split.",
         notes = "R intervals: DeLong; Python intervals: 2,000 bootstrap resamples of the test set.")
top_py <- pyi %>% group_by(learner) %>% slice_max(imp, n = 1)
macro("PyMaxDiff", f3(max(abs(cmp$py_auc - cmp$r_auc))))
macro("PyLogitAuc", f3(cmp$py_auc[1])); macro("PyRfAuc", f3(cmp$py_auc[2])); macro("PyGbAuc", f3(cmp$py_auc[3]))
macro("PyTopFeature", ifelse(all(top_py$feature == "mindfulness"), "mindfulness", paste(unique(top_py$feature), collapse = ", ")))
write_macros("ml_crosscheck")
print(cmp); print(top_py)
