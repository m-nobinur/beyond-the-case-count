# 02_panel.R — Tier 2 (household): women's diet diversity before (Nov 2019) and during (Nov 2020) COVID-19.
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(sandwich); library(lmtest) })

raw <- read_csv(P("data/raw/women_panel_2019_2020/DATA_BNA_BASELINE_COVID_FOLLOWUP_11March21_SH.csv"),
                show_col_types = FALSE, guess_max = 5000)
FG <- c(grain = "Grains, roots, tubers", pulses = "Pulses", nuts = "Nuts and seeds", dairy = "Dairy",
        meatfish = "Meat, poultry, fish", egg = "Eggs", darkveg = "Dark green leafy veg.",
        vita = "Vitamin-A rich fruit/veg.", othveg = "Other vegetables", othfruit = "Other fruits")

base <- raw %>% filter(!is.na(E_Attrited)) %>%
  mutate(retained = E_Attrited == 0 & !is.na(E_nbrfood),
         quintile = case_when(QUINT1 == 1 ~ 1, QUINT2 == 1 ~ 2, QUINT3 == 1 ~ 3, QUINT4 == 1 ~ 4, QUINT5 == 1 ~ 5),
         district = factor(district), works_out_b = as.integer(w6_1 == "Yes"))
pan <- base %>% filter(retained)
stopifnot(all(pan$E_nbrfood - pan$nbrfood == pan$EB_nbrfood, na.rm = TRUE))  # derived change is consistent

# ---- Attrition check (baseline covariates, retained vs. lost) -----------------
att_vars <- c(age_w = "Age (years)", EDUHIGH = "Education above primary", CHILD_W = "Number of children",
              NOHB = "No husband in household", wealthsc = "Wealth score (Equity Tool)",
              works_out_b = "Works outside home", DECFOOD = "Food-purchase decision index",
              nbrfood = "Food groups consumed (0-10)", mdd_w = "Met MDD-W")
att <- map_dfr(names(att_vars), function(v) {
  x1 <- base[[v]][base$retained]; x0 <- base[[v]][!base$retained]
  tt <- t.test(x1, x0)
  smd <- (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
  tibble(Variable = att_vars[[v]], Retained = mean(x1, na.rm = TRUE), Lost = mean(x0, na.rm = TRUE),
         SMD = smd, p = tt$p.value)
})
save_tab(att %>% mutate(p = fp(p)), "tab_panel_attrition", digits = 3, escape = FALSE,
         caption = sprintf("Baseline (2019) characteristics of women retained in the 2020 phone follow-up (n = %d) versus those lost to follow-up (n = %d).", sum(base$retained), sum(!base$retained)),
         notes = "SMD = standardised mean difference; |SMD| < 0.1 is conventionally treated as negligible imbalance. p from Welch t-tests.")

# ---- Paired change in diet diversity ------------------------------------------
dd <- pan %>% filter(!is.na(nbrfood))
tt <- t.test(dd$E_nbrfood, dd$nbrfood, paired = TRUE)
wx <- wilcox.test(dd$E_nbrfood, dd$nbrfood, paired = TRUE, exact = FALSE)
dz <- mean(dd$EB_nbrfood) / sd(dd$EB_nbrfood)
mc_tab <- table(Before = dd$mdd_w, After = dd$E_MDD_W)
mc <- mcnemar.test(mc_tab)

fg <- map_dfr(names(FG), function(g) {
  b <- dd[[g]]; a <- dd[[paste0("E_", g)]]; ok <- !is.na(a) & !is.na(b)
  t2 <- table(factor(b[ok], 0:1), factor(a[ok], 0:1))
  p <- if (t2[1, 2] + t2[2, 1] > 0) mcnemar.test(t2)$p.value else 1
  tibble(group = FG[[g]], key = g, before = mean(b[ok]), after = mean(a[ok]), gained = t2[1, 2], lost = t2[2, 1], p = p)
}) %>% mutate(p_bh = p.adjust(p, "BH"), diff = after - before)

fig_fg <- fg %>% mutate(group = fct_reorder(group, diff)) %>%
  ggplot(aes(y = group)) +
  geom_segment(aes(x = before, xend = after, yend = group), colour = "grey60",
               arrow = arrow(length = unit(.12, "cm"), type = "closed")) +
  geom_point(aes(x = before, colour = "Nov 2019 (in person)"), size = 2) +
  geom_point(aes(x = after, colour = "Nov 2020 (phone)"), size = 2) +
  geom_text(aes(x = pmax(before, after) + .04, label = ifelse(p_bh < .05, "*", "")), size = 4) +
  scale_x_continuous(labels = percent, limits = c(0, 1.05)) +
  scale_colour_manual(values = unname(PAL[c("blue", "orange")]), name = NULL) +
  labs(x = "Share of women consuming the food group (previous 24 h)", y = NULL,
       title = "Food-group consumption before and during COVID-19",
       caption = "* McNemar test, Benjamini–Hochberg adjusted p < 0.05. n = paired women.")
save_fig(fig_fg, "fig_panel_foodgroups", width = 6.3, height = 3.8)

save_tab(fg %>% transmute(`Food group` = group, `2019` = pct(before), `2020` = pct(after),
                          `Gained` = gained, `Lost` = lost, `p (BH)` = fp(p_bh)),
         "tab_panel_foodgroups", escape = FALSE,
         caption = "Percentage of women consuming each of the ten MDD-W food groups in 2019 and 2020, with the number of women who started (gained) or stopped (lost) consuming each group.",
         notes = "Paired McNemar tests with Benjamini–Hochberg false-discovery-rate adjustment across the ten groups.")

# ---- Heterogeneity: who lost diet diversity? (ANCOVA, HC3) --------------------
reg <- dd %>%
  mutate(quintile = factor(quintile), youngw = factor(youngw, 0:1, c("19+", "18-19")),
         EB_WWORK = factor(EB_WWORK, c(0, -1, 1), c("No change", "Stopped", "Started")),
         hb_lowinc = factor(EB_HBLOWINCO, 0:1, c("No", "Yes")),
         hb_lowinc = fct_na_value_to_level(hb_lowinc, "No husband / NA")) %>%
  select(E_nbrfood, nbrfood, quintile, district, age_w, EDUHIGH, CHILD_HH, E_MKTCHALL, E_PPPCHALL,
         hb_lowinc, EB_WWORK, EB_DECFOOD) %>% drop_na()
m_anc <- lm(E_nbrfood ~ nbrfood + quintile + district + age_w + EDUHIGH + CHILD_HH +
              E_MKTCHALL + E_PPPCHALL + hb_lowinc + EB_WWORK + EB_DECFOOD, data = reg)
ct <- coeftest(m_anc, vcov = vcovHC(m_anc, type = "HC3"))
lab_terms <- c("(Intercept)" = "Intercept", nbrfood = "Baseline food groups (2019)",
               quintile2 = "Wealth quintile 2", quintile3 = "Wealth quintile 3", quintile4 = "Wealth quintile 4",
               quintile5 = "Wealth quintile 5 (richest)", districtPatuakhali = "District: Patuakhali",
               age_w = "Age (years)", EDUHIGH = "Education above primary", CHILD_HH = "Children in household",
               E_MKTCHALL = "Market-access constraint (2020)", E_PPPCHALL = "Purchasing-power constraint (2020)",
               `hb_lowincYes` = "Husband's income fell", `hb_lowincNo husband / NA` = "No husband / not reported",
               EB_WWORKStopped = "Stopped outside work", EB_WWORKStarted = "Started outside work",
               EB_DECFOOD = "Change in food-decision index")
tab_anc <- tidy(ct, conf.int = TRUE) %>%
  mutate(term = recode(term, !!!lab_terms),
         `95% CI` = sprintf("[%s, %s]", f2(conf.low), f2(conf.high)), p = fp(p.value)) %>%
  select(Term = term, B = estimate, SE = std.error, `95% CI`, p)
save_tab(tab_anc, "tab_panel_ancova", digits = 2, escape = FALSE,
         caption = "ANCOVA of the number of food groups consumed in 2020, adjusting for 2019 diet diversity (OLS with HC3 robust standard errors).",
         notes = sprintf("n = %d; $R^2$ = %s; reference categories: wealth quintile 1 (poorest), Faridpur, husband's income did not fall, no change in outside work.", nobs(m_anc), f3(summary(m_anc)$r.squared)))

# ---- Change-score model: did the 2019-2020 gain differ by wealth? -------------
# With a baseline coefficient far below 1 the ANCOVA mostly re-expresses the 2019 level gap
# (Lord's paradox), so the gain itself is modelled as well.
m_chg <- lm(I(E_nbrfood - nbrfood) ~ quintile + district, data = mutate(dd, quintile = factor(quintile), district = factor(district)))  # all paired women, as in the quintile table
v_chg <- vcovHC(m_chg, type = "HC3")
qi <- grep("^quintile", names(coef(m_chg)))
w_chg <- as.numeric(t(coef(m_chg)[qi]) %*% solve(v_chg[qi, qi]) %*% coef(m_chg)[qi])
p_chg <- pchisq(w_chg, df = length(qi), lower.tail = FALSE)   # robust joint Wald test of the quintile terms

# by wealth quintile
byq <- dd %>% group_by(quintile) %>%
  summarise(n = n(), b = mean(nbrfood), a = mean(E_nbrfood), gain = mean(E_nbrfood - nbrfood),
            gain_se = sd(E_nbrfood - nbrfood) / sqrt(n()), mdd_b = mean(mdd_w), mdd_a = mean(E_MDD_W),
            ppp = mean(E_PPPCHALL, na.rm = TRUE), mkt = mean(E_MKTCHALL, na.rm = TRUE))
fig_q <- byq %>% select(quintile, `MDD-W 2019` = mdd_b, `MDD-W 2020` = mdd_a,
                        `Purchasing-power constraint 2020` = ppp) %>%
  pivot_longer(-quintile) %>%
  ggplot(aes(factor(quintile), value, fill = name)) + geom_col(position = position_dodge(.8), width = .75) +
  scale_fill_manual(values = unname(PAL[c("sky", "blue", "red")]), name = NULL) +
  scale_y_continuous(labels = percent) +
  labs(x = "Baseline wealth quintile (1 = poorest)", y = "Share of women", title = "Diet adequacy and economic stress by wealth")
save_fig(fig_q, "fig_panel_quintile", width = 6.3, height = 3.4)

save_tab(byq %>% transmute(`Wealth quintile` = quintile, n,
                          `Food groups 2019` = f2(b), `Food groups 2020` = f2(a),
                          `Change [95% CI]` = sprintf("%s [%s, %s]", f2(gain), f2(gain - 1.96 * gain_se), f2(gain + 1.96 * gain_se)),
                          `MDD-W 2019` = pct(mdd_b), `MDD-W 2020` = pct(mdd_a)),
         "tab_panel_byquintile", escape = FALSE,
         caption = "Diet diversity and MDD-W attainment by baseline wealth quintile, 2019 and 2020.",
         notes = sprintf("Robust (HC3) joint Wald test that the 2019--2020 change is equal across quintiles, adjusting for district: $\\\\chi^2$(4) = %s, p = %s.", f2(w_chg), f3(p_chg)))

# ---- Macros ------------------------------------------------------------------
macro("PanBase", fint(nrow(base))); macro("PanN", fint(nrow(dd))); macro("PanLost", fint(sum(!base$retained)))
macro("PanRetention", pct(mean(base$retained)))
macro("PanMaxSMD", f2(max(abs(att$SMD))))
macro("PanFgBefore", f2(mean(dd$nbrfood))); macro("PanFgAfter", f2(mean(dd$E_nbrfood)))
macro("PanFgDiff", f2(mean(dd$EB_nbrfood))); macro("PanFgT", f2(tt$statistic)); macro("PanFgP", fp(tt$p.value))
macro("PanFgCIlo", f2(tt$conf.int[1])); macro("PanFgCIhi", f2(tt$conf.int[2])); macro("PanDz", f2(dz))
macro("PanWilcoxP", fp(wx$p.value))
macro("PanMddBefore", pct(mean(dd$mdd_w))); macro("PanMddAfter", pct(mean(dd$E_MDD_W)))
macro("PanMddGain", mc_tab["0", "1"]); macro("PanMddLoss", mc_tab["1", "0"])
macro("PanMcChi", f2(mc$statistic)); macro("PanMcP", fp(mc$p.value))
macro("PanDeclined", pct(mean(dd$EB_nbrfood < 0))); macro("PanImproved", pct(mean(dd$EB_nbrfood > 0)))
macro("PanPPP", pct(mean(dd$E_PPPCHALL, na.rm = TRUE))); macro("PanMKT", pct(mean(dd$E_MKTCHALL, na.rm = TRUE)))
macro("PanHbLow", pct(mean(dd$EB_HBLOWINCO == 1, na.rm = TRUE)))
macro("PanAncN", nobs(m_anc)); macro("PanAncR", f3(summary(m_anc)$r.squared))
macro("PanPPPb", f2(coef(m_anc)["E_PPPCHALL"])); macro("PanPPPp", fp(ct["E_PPPCHALL", 4]))
macro("PanMKTb", f2(coef(m_anc)["E_MKTCHALL"])); macro("PanMKTp", fp(ct["E_MKTCHALL", 4]))
macro("PanQFiveb", f2(coef(m_anc)["quintile5"])); macro("PanQFivep", fp(ct["quintile5", 4]))
macro("PanBaseb", f2(coef(m_anc)["nbrfood"]))
macro("PanRaw", fint(nrow(raw)))
macro("PanChgP", fp(p_chg)); macro("PanChgW", f2(w_chg))
for (q in 1:5) macro(paste0("PanGainQ", c("One", "Two", "Three", "Four", "Five")[q]), f2(byq$gain[q]))
macro("PanQOneMddB", pct(byq$mdd_b[1])); macro("PanQFiveMddB", pct(byq$mdd_b[5]))
macro("PanQOneFgB", f2(byq$b[1])); macro("PanQFiveFgB", f2(byq$b[5]))
top_fg <- fg %>% slice_max(diff, n = 1); low_fg <- fg %>% slice_min(diff, n = 1)
macro("PanTopFg", top_fg$group); macro("PanTopFgB", pct(top_fg$before)); macro("PanTopFgA", pct(top_fg$after))
macro("PanLowFg", low_fg$group); macro("PanLowFgB", pct(low_fg$before)); macro("PanLowFgA", pct(low_fg$after))
macro("PanNsigFg", sum(fg$p_bh < .05))
macro("PanQOneMddA", pct(byq$mdd_a[1])); macro("PanQFiveMddA", pct(byq$mdd_a[5]))
macro("PanQOnePPP", pct(byq$ppp[1])); macro("PanQFivePPP", pct(byq$ppp[5]))
write_macros("panel")
cat("02_panel done: n =", nrow(dd), " diff =", mean(dd$EB_nbrfood), " change-by-quintile p =", p_chg, "\n")
print(fg); print(byq)
