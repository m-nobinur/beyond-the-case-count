# 03_mh_explore.R — Tier 3 (individual): cleaning, psychometrics, prevalence and explanatory models
# for the May–June 2021 Bangladesh mental-health survey (Islam et al., 2023 data).
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(psych) })

raw <- read_sav(P("data/raw/mental_health_2021/mh2021.sav"))
stopifnot(nrow(raw) == 1899)
n_third <- sum(raw$Q2 == 3)

# ---- Verify the file's scale totals against the items (guards against mis-scoring)
it <- function(prefix, d = raw) d %>% select(starts_with(prefix)) %>% zap_labels() %>% as.matrix()
stopifnot(all(rowSums(it("Q24_")) == raw$Anxiety),            # Beck Anxiety Inventory (21 items, 0-63)
          all(rowSums(it("Q27_")) == raw$Depression),         # PHQ-9 (0-27)
          all(rowSums(it("Q26_")) == raw$Stress),             # PSS-4 as scored by the authors
          all(rowSums(it("Q29_")) == raw$Mindfulness),        # MAAS-5 (higher = more mindful)
          all(rowSums(zap_labels(raw[, c("Q25_1", "Q25_4", "Q25_10")])) == raw$Optimism),
          all(rowSums(zap_labels(raw[, c("Q25_3", "Q25_7", "Q25_9")])) == raw$Pessimism))

# ---- Clean analysis frame -----------------------------------------------------
# Third-gender respondents (n = 2) are excluded, as in the original study: two cases cannot support estimation.
mh <- raw %>% zap_labels() %>% filter(Q2 %in% c(0, 1)) %>%
  transmute(
    id = NewID,
    sex = factor(Q2, 0:1, c("Male", "Female")),
    age = Q3,
    education = factor(case_when(Q4 <= 2 ~ "SSC or below", Q4 == 3 ~ "HSC", Q4 == 4 ~ "Bachelor", Q4 == 5 ~ "Postgraduate"),
                       c("SSC or below", "HSC", "Bachelor", "Postgraduate")),
    mh_history = factor(Q5, c(0, 1, 3), c("No", "Yes", "Prefer not to say")),
    division = factor(Q6, 1:8, c("Dhaka", "Barishal", "Rajshahi", "Rangpur", "Sylhet", "Khulna", "Chattogram", "Mymensingh")),
    residence = factor(Q7, 0:1, c("Rural", "Urban")),
    living = factor(case_when(Q8 == 1 ~ "Family", Q8 == 2 ~ "Hall/hostel", Q8 == 3 ~ "Mess/flat", TRUE ~ "Alone/other"),
                    c("Family", "Hall/hostel", "Mess/flat", "Alone/other")),
    profession = factor(Q9, 1:4, c("Paid service", "Business/farming", "Non-paid service", "Student")),
    frontliner = factor(Q10, 0:1, c("No", "Yes")),
    ses = Q11,
    married = factor(Q13, 0:1, c("No", "Yes")),
    children = factor(coalesce(Q14, 0), 0:1, c("No", "Yes")),             # structurally missing for unmarried
    covid_self = factor(case_when(Q17 == 0 ~ "Not infected", Q19 == 1 ~ "Infected, hospital", TRUE ~ "Infected, home"),
                        c("Not infected", "Infected, home", "Infected, hospital")),
    covid_family = factor(case_when(Q20 == 0 ~ "None", coalesce(Q21, 0) == 1 ~ "Infected, a death", TRUE ~ "Infected, no death"),
                          c("None", "Infected, no death", "Infected, a death")),
    vaccine = factor(Q22, c(3, 1, 2), c("None", "One dose", "Two doses")),
    anxiety = Anxiety, depression = Depression, stress = Stress,
    optimism = Optimism, pessimism = Pessimism, mindfulness = Mindfulness,
    phq2 = Q27_1 + Q27_2,                  # PHQ-2 (anhedonia + low mood): a reference ceiling only, never a feature
    dep = as.integer(Depression >= 10),    # PHQ-9 >= 10 (Manea et al., 2012)
    anx = as.integer(Anxiety >= 16)        # BAI >= 16, moderate or worse (Beck et al., 1988)
  )
items <- raw %>% filter(Q2 %in% c(0, 1)) %>% zap_labels()
# Data-quality flag: an identical response block shared by >= 5 respondents is implausible for genuine
# answers (templated entries). The block uses only non-outcome items (all 10 LOT-R items incl. fillers,
# PSS-4, MAAS-5) so that the flag does not condition on depression or anxiety.
flag_block <- function(cols, k = 5) { b <- apply(items[, cols], 1, paste, collapse = ""); as.vector(table(b)[b] >= k) }
NONOUT <- grep("^Q2[569]_", names(items), value = TRUE)
stopifnot(length(NONOUT) == 19)
mh$suspect <- flag_block(NONOUT)
flag_alt <- tibble(
  Definition = c("Non-outcome items: LOT-R, PSS-4, MAAS (primary)", "Primary + PHQ-9 items and impairment item",
                 "Primary + BAI items"),
  Items = c(19, 29, 40),
  flag = list(mh$suspect, flag_block(c(NONOUT, paste0("Q27_", 1:9), "Q28")), flag_block(c(NONOUT, paste0("Q24_", 1:21))))) %>%
  mutate(Flagged = map_int(flag, sum), `From Chattogram` = map_int(flag, ~ sum(mh$division[.x] == "Chattogram")),
         `PHQ-9 ≥ 10` = map_int(flag, ~ sum(mh$dep[.x])), `Overlap with primary` = map_int(flag, ~ sum(.x & mh$suspect)))
n_patterns <- length(unique(apply(items[mh$suspect, NONOUT], 1, paste, collapse = "")))
stopifnot(nrow(mh) == 1897, !anyNA(mh %>% select(-id)))
write_csv(mh, P("data/processed/mh_clean.csv"))
saveRDS(mh, P("data/processed/mh_clean.rds"))

save_tab(flag_alt %>% select(-flag), "tab_suspect_defs", font = "footnotesize",
         caption = "Records flagged as templated under three definitions of the response block (identical block shared by at least five respondents).",
         notes = "The primary definition excludes all outcome items, so the flag cannot select on depression or anxiety by construction.")

# ---- Reliability ---------------------------------------------------------------
scales_items <- list(
  "Beck Anxiety Inventory (BAI)" = paste0("Q24_", 1:21),
  "Patient Health Questionnaire (PHQ-9)" = paste0("Q27_", 1:9),
  "Perceived Stress Scale (PSS-4)" = paste0("Q26_", 1:4),
  "LOT-R optimism" = c("Q25_1", "Q25_4", "Q25_10"),
  "LOT-R pessimism" = c("Q25_3", "Q25_7", "Q25_9"),
  "Mindful Attention Awareness (MAAS-5)" = paste0("Q29_", 1:5))
rel <- imap_dfr(scales_items, function(v, nm) {
  x <- items[, v]
  a <- suppressWarnings(psych::alpha(x, check.keys = FALSE, warnings = FALSE))
  om <- if (length(v) >= 3) suppressWarnings(suppressMessages(psych::omega(x, nfactors = 1, plot = FALSE)))$omega.tot else NA
  tibble(Scale = nm, Items = length(v), alpha = a$total$raw_alpha, omega = om,
         `Mean (SD)` = sprintf("%s (%s)", f2(mean(rowSums(x))), f2(sd(rowSums(x)))),
         Range = sprintf("%d–%d", min(rowSums(x)), max(rowSums(x))))
})
save_tab(rel %>% rename(`Cronbach's α` = alpha, `McDonald's ω` = omega), "tab_reliability", digits = 2,
         caption = "Internal consistency and descriptive statistics of the six psychometric scales (N = 1,897).",
         notes = "α ≥ 0.70 is conventionally acceptable. ω total from a one-factor model.")

# ---- Table 1: sample description by depression status ---------------------------
cat_vars <- c(sex = "Sex", education = "Education", mh_history = "History of mental disturbance",
              division = "Division", residence = "Residence", living = "Living arrangement",
              profession = "Profession", frontliner = "Frontline worker", married = "Married",
              children = "Has children", covid_self = "Own COVID-19 infection",
              covid_family = "Family COVID-19 experience", vaccine = "Vaccination status")
t1 <- map_dfr(names(cat_vars), function(v) {
  tab <- table(mh[[v]], mh$dep)
  p <- suppressWarnings(chisq.test(tab))$p.value
  tibble(Variable = c(cat_vars[[v]], paste0("\\quad ", levels(mh[[v]]))),
         `Total n (\\%)` = c("", sprintf("%d (%s)", as.vector(table(mh[[v]])), pct(prop.table(table(mh[[v]]))))),
         `Depressed \\% (row)` = c("", pct(tab[, 2] / rowSums(tab))),
         p = c(fp(p), rep("", nlevels(mh[[v]]))))
})
t1 <- bind_rows(
  tibble(Variable = "Age, mean (SD)", `Total n (\\%)` = sprintf("%s (%s)", f1(mean(mh$age)), f1(sd(mh$age))),
         `Depressed \\% (row)` = "", p = fp(t.test(age ~ dep, mh)$p.value)),
  tibble(Variable = "SES ladder (0–10), mean (SD)", `Total n (\\%)` = sprintf("%s (%s)", f1(mean(mh$ses)), f1(sd(mh$ses))),
         `Depressed \\% (row)` = "", p = fp(t.test(ses ~ dep, mh)$p.value)),
  t1)
writeLines(as.character(
  kbl(t1, format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE, linesep = "",
      caption = "Characteristics of the analytic sample and row percentage meeting the PHQ-9 $\\geq$ 10 depression threshold (N = 1,897).",
      label = "tab_table1", align = "lrrr") %>%
    kable_styling(latex_options = c("repeat_header"), font_size = 9) %>%
    footnote(general = "p from Pearson $\\\\chi^2$ tests (categorical) or Welch t-tests (continuous).", general_title = "Note.", escape = FALSE)),
  P("output", "tables", "tab_table1.tex"))

# ---- Prevalence with Wilson CIs -------------------------------------------------
prev_overall <- tibble(outcome = c("Depression (PHQ-9 ≥ 10)", "Anxiety (BAI ≥ 16)", "Both", "Either"),
                       x = c(sum(mh$dep), sum(mh$anx), sum(mh$dep & mh$anx), sum(mh$dep | mh$anx))) %>%
  bind_cols(map_dfr(.$x, ~ wilson(.x, nrow(mh))))
grp_vars <- c("sex", "profession", "residence", "mh_history", "covid_self", "covid_family", "vaccine", "division")
prev_grp <- map_dfr(grp_vars, function(v) {
  mh %>% group_by(level = .data[[v]]) %>%
    summarise(n = n(), Depression = sum(dep), Anxiety = sum(anx)) %>%
    pivot_longer(c(Depression, Anxiety), names_to = "outcome", values_to = "x") %>%
    mutate(var = cat_vars[[v]]) %>% bind_cols(map2_dfr(.$x, .$n, wilson))
})
fig_prev <- prev_grp %>%
  mutate(var = factor(var, unique(var)), level = fct_rev(fct_inorder(as.character(level)))) %>%
  ggplot(aes(p, level, colour = outcome)) +
  geom_vline(data = prev_overall %>% filter(outcome %in% c("Depression (PHQ-9 ≥ 10)", "Anxiety (BAI ≥ 16)")) %>%
               mutate(outcome = word(outcome)), aes(xintercept = p, colour = outcome), linetype = 2, linewidth = .3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(.5), size = .2) +
  facet_grid(var ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(labels = percent) +
  scale_colour_manual(values = c(Depression = PAL[["blue"]], Anxiety = PAL[["orange"]]), name = NULL) +
  labs(x = "Prevalence (95% Wilson CI)", y = NULL, title = "Probable depression and anxiety by subgroup") +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7.5), strip.placement = "outside",
        axis.text.y = element_text(size = 7.5), panel.spacing = unit(2, "pt"))
save_fig(fig_prev, "fig_mh_prevalence", width = 6.5, height = 8.2)

# ---- Symptom-level view: PHQ-9 and top BAI items -------------------------------
phq_lab <- c("Anhedonia", "Low mood", "Sleep problems", "Fatigue", "Appetite change", "Worthlessness",
             "Poor concentration", "Psychomotor change", "Suicidal ideation")
phq_items <- items %>% select(paste0("Q27_", 1:9)) %>% set_names(phq_lab) %>%
  pivot_longer(everything()) %>% count(name, value) %>% group_by(name) %>% mutate(prop = n / sum(n)) %>%
  mutate(value = factor(value, 3:0, c("Nearly every day", "More than half the days", "Several days", "Not at all")))
ord <- phq_items %>% filter(value != "Not at all") %>% group_by(name) %>% summarise(s = sum(prop)) %>% arrange(s) %>% pull(name)
fig_phq <- ggplot(phq_items, aes(prop, factor(name, ord), fill = value)) + geom_col(width = .75) +
  scale_fill_manual(values = c("#08306B", "#2171B5", "#9ECAE1", "#EEEEEE"), name = NULL) +
  scale_x_continuous(labels = percent) + guides(fill = guide_legend(reverse = TRUE, nrow = 2)) +
  labs(x = "Share of respondents (past two weeks)", y = NULL, title = "PHQ-9 symptom endorsement")
save_fig(fig_phq, "fig_mh_phq_items", width = 6.3, height = 3.6)
suicidal <- mean(items$Q27_9 > 0)

# ---- Correlations among psychological constructs -------------------------------
sc <- mh %>% select(Anxiety = anxiety, Depression = depression, Stress = stress,
                    Optimism = optimism, Pessimism = pessimism, Mindfulness = mindfulness, Age = age, SES = ses)
cm <- cor(sc, method = "spearman")
fig_cor <- ggcorrplot::ggcorrplot(cm, type = "lower", lab = TRUE, lab_size = 2.8, outline.color = "white",
                                  colors = c(PAL[["red"]], "white", PAL[["blue"]]), legend.title = "Spearman ρ") +
  labs(title = "Correlations among mental-health and resource measures") + theme(plot.title = element_text(face = "bold", size = 11))
save_fig(fig_cor, "fig_mh_correlation", width = 5.6, height = 5)

# ---- Explanatory multivariable logistic regression -----------------------------
fA <- ~ sex + age + education + mh_history + division + residence + living + profession + frontliner +
  ses + married + children + covid_self + covid_family + vaccine
fB <- update(fA, ~ . + optimism + pessimism + mindfulness)
fit_or <- function(y, f) {
  m <- glm(update(f, paste(y, "~ .")), data = mh, family = binomial)
  ci <- suppressMessages(confint.default(m))
  tibble(term = names(coef(m)), OR = exp(coef(m)), lo = exp(ci[, 1]), hi = exp(ci[, 2]),
         p = summary(m)$coefficients[, 4]) %>% filter(term != "(Intercept)") %>%
    mutate(q = p.adjust(p, "BH")) %>% list(tab = ., model = m)
}
dA <- fit_or("dep", fA); dB <- fit_or("dep", fB); aB <- fit_or("anx", fB)
pretty_term <- function(t) {
  t %>% str_replace("^sex", "Sex: ") %>% str_replace("^education", "Education: ") %>%
    str_replace("^mh_history", "MH history: ") %>% str_replace("^division", "Division: ") %>%
    str_replace("^residence", "Residence: ") %>% str_replace("^living", "Living: ") %>%
    str_replace("^profession", "Profession: ") %>% str_replace("^frontliner", "Frontliner: ") %>%
    str_replace("^married", "Married: ") %>% str_replace("^children", "Children: ") %>%
    str_replace("^covid_self", "Own COVID: ") %>% str_replace("^covid_family", "Family COVID: ") %>%
    str_replace("^vaccine", "Vaccine: ") %>%
    recode(age = "Age (per year)", ses = "SES ladder (per step)", optimism = "Optimism (per point)",
           pessimism = "Pessimism (per point)", mindfulness = "Mindfulness (per point)")
}
or_tab <- dB$tab %>% select(term, OR, lo, hi, q) %>%
  left_join(aB$tab %>% select(term, OR_a = OR, lo_a = lo, hi_a = hi, q_a = q), by = "term") %>%
  transmute(Predictor = pretty_term(term),
            `Depression aOR [95\\% CI]` = sprintf("%s [%s, %s]", f2(OR), f2(lo), f2(hi)), `q (dep.)` = fp(q),
            `Anxiety aOR [95\\% CI]` = sprintf("%s [%s, %s]", f2(OR_a), f2(lo_a), f2(hi_a)), `q (anx.)` = fp(q_a))
writeLines(as.character(
  kbl(or_tab, format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE, linesep = "", align = "lrrrr",
      caption = "Adjusted odds ratios from multivariable logistic regressions of probable depression and anxiety on sociodemographic, COVID-19 exposure and psychological-scale variables (N = 1,897).",
      label = "tab_or") %>%
    kable_styling(latex_options = "repeat_header", font_size = 8.5) %>%
    footnote(general = "q = Benjamini--Hochberg adjusted p-value. Reference categories: male, SSC or below, no MH history, Dhaka, rural, living with family, paid service, not a frontliner, unmarried, no children, not infected, no family infection, unvaccinated.",
             general_title = "Note.", threeparttable = TRUE, escape = FALSE)),
  P("output", "tables", "tab_or.tex"))

fig_or <- dB$tab %>% mutate(model = "Depression") %>% bind_rows(aB$tab %>% mutate(model = "Anxiety")) %>%
  mutate(term = pretty_term(term), sig = q < .05) %>%
  filter(term %in% (bind_rows(dB$tab, aB$tab) %>% filter(q < .1) %>% pull(term) %>% pretty_term())) %>%
  ggplot(aes(OR, fct_rev(factor(term)), colour = model, alpha = sig)) +
  geom_vline(xintercept = 1, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(.55), size = .25) +
  scale_x_log10() + scale_alpha_manual(values = c(`FALSE` = .35, `TRUE` = 1), guide = "none") +
  scale_colour_manual(values = c(Depression = PAL[["blue"]], Anxiety = PAL[["orange"]]), name = NULL) +
  labs(x = "Adjusted odds ratio (log scale, 95% CI)", y = NULL,
       title = "Adjusted associations (predictors with q < 0.10 in either model)",
       caption = "Faded estimates: q ≥ 0.05 in that model.")
save_fig(fig_or, "fig_mh_or", width = 6.5, height = 5.2)

or_get <- function(tab, t, what = "OR") f2(tab$tab[[what]][tab$tab$term == t])
# ---- Prevalence ratios (modified Poisson, robust SE): odds ratios overstate risk at 45% prevalence ----
fit_pr <- function(y, f) {
  m <- glm(update(f, paste(y, "~ .")), data = mh, family = poisson)
  V <- sandwich::vcovHC(m, type = "HC0"); se <- sqrt(diag(V))
  tibble(term = names(coef(m)), PR = exp(coef(m)), lo = exp(coef(m) - 1.96 * se), hi = exp(coef(m) + 1.96 * se)) %>%
    filter(term != "(Intercept)")
}
prA <- fit_pr("dep", fA); prB <- fit_pr("dep", fB)
pr_terms <- c("mh_historyYes", "professionStudent", "ses", "covid_familyInfected, a death", "covid_selfInfected, hospital",
              "vaccineTwo doses", "sexFemale", "optimism", "mindfulness")
ci_fmt <- function(tab, t, est) { r <- tab[tab$term == t, ]; if (!nrow(r)) return("--"); sprintf("%s [%s, %s]", f2(r[[est]]), f2(r$lo), f2(r$hi)) }
tab_pr <- tibble(Predictor = pretty_term(pr_terms),
                 `aOR, model A` = map_chr(pr_terms, ~ ci_fmt(dA$tab, .x, "OR")),
                 `aOR, model B` = map_chr(pr_terms, ~ ci_fmt(dB$tab, .x, "OR")),
                 `aPR, model A` = map_chr(pr_terms, ~ ci_fmt(prA, .x, "PR")),
                 `aPR, model B` = map_chr(pr_terms, ~ ci_fmt(prB, .x, "PR")))
save_tab(tab_pr, "tab_mh_pr", font = "footnotesize", escape = FALSE,
         caption = "Probable depression: adjusted odds ratios (aOR) and adjusted prevalence ratios (aPR) for key predictors, without (model A) and with (model B) the psychological scales.",
         notes = "aPR from modified Poisson regression with robust (HC0) standard errors. Model A adjusts for all sociodemographic and COVID-19 exposure variables; model B adds optimism, pessimism and mindfulness, which may lie on the causal path. Values are estimate [95\\\\% CI].")
prg <- function(tab, t) f2(tab$PR[tab$term == t])
macro("PrFamDeathA", prg(prA, "covid_familyInfected, a death")); macro("PrFamDeathB", prg(prB, "covid_familyInfected, a death"))
macro("PrHistB", prg(prB, "mh_historyYes")); macro("PrTwoDoseB", prg(prB, "vaccineTwo doses")); macro("PrTwoDoseA", prg(prA, "vaccineTwo doses"))
macro("PrStudentB", prg(prB, "professionStudent")); macro("PrMindB", prg(prB, "mindfulness"))
macro("OrTwoDoseA", or_get(dA, "vaccineTwo doses")); macro("OrHistA", or_get(dA, "mh_historyYes")); macro("OrStudentA", or_get(dA, "professionStudent"))
macro("OrHospB", or_get(dB, "covid_selfInfected, hospital")); macro("OrHospLo", or_get(dB, "covid_selfInfected, hospital", "lo"))
macro("OrHospHi", or_get(dB, "covid_selfInfected, hospital", "hi")); macro("MhHospN", sum(mh$covid_self == "Infected, hospital"))

# Incremental value of psychological resources (likelihood-ratio test and Nagelkerke R2)
lr <- anova(dA$model, dB$model, test = "LRT")
nag <- function(m) { n <- nobs(m); l0 <- -m$null.deviance / 2; l1 <- as.numeric(logLik(m))
  (1 - exp(2 * (l0 - l1) / n)) / (1 - exp(2 * l0 / n)) }

# ---- Sensitivity: refit model B without the flagged records -------------------
dB_s <- { m <- glm(update(fB, dep ~ .), data = filter(mh, !suspect), family = binomial)
  ci <- confint.default(m); tibble(term = names(coef(m)), OR = exp(coef(m)), lo = exp(ci[, 1]), hi = exp(ci[, 2])) }
aB_s <- { m <- glm(update(fB, anx ~ .), data = filter(mh, !suspect), family = binomial)
  ci <- confint.default(m); tibble(term = names(coef(m)), OR = exp(coef(m)), lo = exp(ci[, 1]), hi = exp(ci[, 2])) }
key_terms <- c("mh_historyYes", "professionStudent", "ses", "divisionChattogram", "covid_familyInfected, a death",
               "vaccineTwo doses", "optimism", "mindfulness")
fmt_or <- function(t, tt) { r <- tt[tt$term == t, ]; sprintf("%s [%s, %s]", f2(r$OR), f2(r$lo), f2(r$hi)) }
sens <- tibble(Predictor = pretty_term(key_terms),
               `Depression: all` = map_chr(key_terms, fmt_or, dB$tab), `Depression: excl.` = map_chr(key_terms, fmt_or, dB_s),
               `Anxiety: all` = map_chr(key_terms, fmt_or, aB$tab), `Anxiety: excl.` = map_chr(key_terms, fmt_or, aB_s))
save_tab(sens, "tab_mh_sensitivity", font = "footnotesize",
         caption = sprintf("Sensitivity of key adjusted odds ratios (model B) to excluding the %d respondents with a templated response block (N = %s vs.\\ %s).", sum(mh$suspect), fint(nrow(mh)), fint(sum(!mh$suspect))),
         escape = FALSE, notes = "Values are aOR [95\\\\% CI]. ``Excl.'' omits respondents whose 19 LOT-R, PSS-4 and MAAS item responses were identical to those of at least four others.")
macro("SuspectN", sum(mh$suspect)); macro("SuspectChatto", sum(mh$suspect & mh$division == "Chattogram"))
macro("SuspectTwoDose", pct(mean(mh$vaccine[mh$suspect] == "Two doses")))
macro("CleanN", fint(sum(!mh$suspect)))
macro("SuspectPatterns", n_patterns); macro("SuspectMaxPhq", max(mh$depression[mh$suspect]))
macro("SuspectPhqZero", sum(mh$depression[mh$suspect] == 0))
macro("SuspectAltPhq", flag_alt$Flagged[2]); macro("SuspectAltBai", flag_alt$Flagged[3])
macro("SuspectDep", sum(mh$dep[mh$suspect])); macro("SuspectAnx", sum(mh$anx[mh$suspect]))
macro("SensChattoAll", f2(dB$tab$OR[dB$tab$term == "divisionChattogram"])); macro("SensChattoExcl", f2(dB_s$OR[dB_s$term == "divisionChattogram"]))
macro("SensVaxExcl", f2(dB_s$OR[dB_s$term == "vaccineTwo doses"]))
oth <- setdiff(key_terms, c("divisionChattogram", "vaccineTwo doses"))
macro("SensMaxOther", f2(max(abs(round(dB$tab$OR[match(oth, dB$tab$term)], 2) - round(dB_s$OR[match(oth, dB_s$term)], 2)))))  # as shown in the table

# ---- Macros ------------------------------------------------------------------
macro("MhRaw", fint(nrow(raw))); macro("MhN", fint(nrow(mh))); macro("MhThird", n_third)
macro("MhFemale", pct(mean(mh$sex == "Female"))); macro("MhAgeMean", f1(mean(mh$age))); macro("MhAgeSD", f1(mean(sd(mh$age))))
macro("MhStudents", pct(mean(mh$profession == "Student"))); macro("MhUrban", pct(mean(mh$residence == "Urban")))
macro("MhDhaka", pct(mean(mh$division == "Dhaka"))); macro("MhInfected", pct(mean(mh$covid_self != "Not infected")))
macro("MhFamInf", pct(mean(mh$covid_family != "None"))); macro("MhFamDeath", pct(mean(mh$covid_family == "Infected, a death")))
macro("MhTwoDose", pct(mean(mh$vaccine == "Two doses"))); macro("MhNoVax", pct(mean(mh$vaccine == "None")))
macro("MhHistYes", pct(mean(mh$mh_history == "Yes")))
po <- prev_overall
macro("PrevDep", pct(po$p[1])); macro("PrevDepLo", pct(po$lo[1])); macro("PrevDepHi", pct(po$hi[1]))
macro("PrevAnx", pct(po$p[2])); macro("PrevAnxLo", pct(po$lo[2])); macro("PrevAnxHi", pct(po$hi[2]))
macro("PrevBoth", pct(po$p[3])); macro("PrevEither", pct(po$p[4]))
macro("Suicidal", pct(suicidal))
macro("AlphaBAI", f2(rel$alpha[1])); macro("AlphaPHQ", f2(rel$alpha[2])); macro("AlphaPSS", f2(rel$alpha[3]))
macro("AlphaOpt", f2(rel$alpha[4])); macro("AlphaPess", f2(rel$alpha[5])); macro("AlphaMAAS", f2(rel$alpha[6]))
macro("RhoAnxDep", f2(cm["Anxiety", "Depression"])); macro("RhoMindDep", f2(cm["Mindfulness", "Depression"]))
macro("RhoOptDep", f2(cm["Optimism", "Depression"]))
gp <- function(var, lev, out) prev_grp %>% filter(var == cat_vars[[!!var]], level == lev, outcome == out) %>% pull(p) %>% pct()
macro("DepFemale", gp("sex", "Female", "Depression")); macro("DepMale", gp("sex", "Male", "Depression"))
macro("DepStudent", gp("profession", "Student", "Depression")); macro("DepPaid", gp("profession", "Paid service", "Depression"))
macro("DepHistYes", gp("mh_history", "Yes", "Depression")); macro("DepHistNo", gp("mh_history", "No", "Depression"))
macro("DepFamDeath", gp("covid_family", "Infected, a death", "Depression")); macro("DepFamNone", gp("covid_family", "None", "Depression"))
macro("DepHosp", gp("covid_self", "Infected, hospital", "Depression"))
macro("DepVaxTwo", gp("vaccine", "Two doses", "Depression")); macro("DepVaxNone", gp("vaccine", "None", "Depression"))
divp <- prev_grp %>% filter(var == "Division", outcome == "Depression") %>% arrange(p)
for (d in levels(mh$division)) macro(paste0("DepDiv", d), pct(divp$p[divp$level == d]))
mid <- divp %>% filter(!level %in% c("Chattogram", "Rangpur", "Rajshahi"))
macro("DepDivMidLo", pct(min(mid$p), 0)); macro("DepDivMidHi", pct(max(mid$p), 0))
macro("DivNLo", min(table(mh$division)[setdiff(levels(mh$division), c("Dhaka", "Chattogram"))]))
macro("DivNHi", max(table(mh$division)[setdiff(levels(mh$division), c("Dhaka", "Chattogram"))]))
macro("OrMind", or_get(dB, "mindfulness")); macro("OrMindLo", or_get(dB, "mindfulness", "lo")); macro("OrMindHi", or_get(dB, "mindfulness", "hi"))
macro("OrOpt", or_get(dB, "optimism")); macro("OrOptLo", or_get(dB, "optimism", "lo")); macro("OrOptHi", or_get(dB, "optimism", "hi"))
macro("OrHist", or_get(dB, "mh_historyYes")); macro("OrHistLo", or_get(dB, "mh_historyYes", "lo")); macro("OrHistHi", or_get(dB, "mh_historyYes", "hi"))
macro("OrFemaleA", or_get(dA, "sexFemale")); macro("OrFemaleB", or_get(dB, "sexFemale"))
macro("OrFamDeathA", or_get(dA, "covid_familyInfected, a death")); macro("OrFamDeathB", or_get(dB, "covid_familyInfected, a death"))
macro("OrFamDeathBq", fp(dB$tab$q[dB$tab$term == "covid_familyInfected, a death"]))
macro("OrTwoDose", or_get(dB, "vaccineTwo doses")); macro("OrTwoDoseq", fp(dB$tab$q[dB$tab$term == "vaccineTwo doses"]))
macro("OrStudent", or_get(dB, "professionStudent")); macro("OrStudentq", fp(dB$tab$q[dB$tab$term == "professionStudent"]))
macro("OrAnxMind", or_get(aB, "mindfulness")); macro("OrAnxHist", or_get(aB, "mh_historyYes"))
or3 <- function(nm, tab, t) { macro(nm, or_get(tab, t)); macro(paste0(nm, "Lo"), or_get(tab, t, "lo")); macro(paste0(nm, "Hi"), or_get(tab, t, "hi")) }
or3("OrFamInf", dB, "covid_familyInfected, no death"); or3("OrHome", dB, "covid_selfInfected, home")
or3("OrSes", dB, "ses"); or3("OrPess", dB, "pessimism")
or3("OrAnxTwoDose", aB, "vaccineTwo doses"); or3("OrAnxFemale", aB, "sexFemale")
macro("LrChi", f1(lr$Deviance[2])); macro("LrDf", lr$Df[2]); macro("LrP", fp(lr$`Pr(>Chi)`[2]))
macro("NagA", f3(nag(dA$model))); macro("NagB", f3(nag(dB$model)))
macro("NsigDep", sum(dB$tab$q < .05)); macro("NsigAnx", sum(aB$tab$q < .05))
pc <- cor(items %>% select(starts_with("Q26_")))
pss_phq <- cor(items %>% select(starts_with("Q26_")), items$Depression)
pss_flip <- items %>% select(starts_with("Q26_")) %>% mutate(Q26_2 = 4 - Q26_2, Q26_3 = 4 - Q26_3)
macro("PssPhqMin", f2(min(pss_phq))); macro("PssPhqMax", f2(max(pss_phq)))
macro("PssAlphaFlip", f2(suppressWarnings(psych::alpha(pss_flip, check.keys = FALSE, warnings = FALSE))$total$raw_alpha))
macro("PssRlo", f2(min(pc[c(1, 4), c(2, 3)]))); macro("PssRhi", f2(max(pc[c(1, 4), c(2, 3)]))); macro("PssRin", f2(pc[2, 3]))
macro("FamDeathMissing", sum(items$Q20 == 1 & is.na(items$Q21))); macro("TreatMissing", sum(items$Q17 == 1 & is.na(items$Q19)))
macro("OrMindPct", pct(1 - dB$tab$OR[dB$tab$term == "mindfulness"], 0)); macro("OrOptPct", pct(1 - dB$tab$OR[dB$tab$term == "optimism"], 0))
write_macros("mh_explore")
saveRDS(list(rel = rel, dA = dA$tab, dB = dB$tab, aB = aB$tab, prev = prev_overall), P("data/processed/mh_explore.rds"))
cat("03_mh_explore done. Dep:", pct(po$p[1]), " Anx:", pct(po$p[2]), "\n"); print(rel)
print(dB$tab %>% filter(q < .05)); print(aB$tab %>% filter(q < .05))
