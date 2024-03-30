# 01_national.R — Tier 1 (macro): the national epidemic, policy and mobility context, 2020–2022.
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({ library(changepoint); library(sandwich); library(lmtest) })

# ---- Load & merge ------------------------------------------------------------
owid <- read_csv(P("data/raw/national/owid_bgd.csv"), show_col_types = FALSE) %>%
  transmute(date, cases = new_cases_smoothed, deaths = new_deaths_smoothed,
            total_cases, total_deaths, positivity = positive_rate,
            fully_vax = people_fully_vaccinated_per_hundred, population)
ox <- read_csv(P("data/raw/national/oxcgrt_bgd.csv"), show_col_types = FALSE) %>%
  filter(Jurisdiction == "NAT_TOTAL") %>%
  transmute(date = ymd(Date), stringency = StringencyIndex_Average,
            econ_support = EconomicSupportIndex)
mob <- read_csv(P("data/raw/national/mobility_BD_national.csv"), show_col_types = FALSE) %>%
  transmute(date,
            Retail = retail_and_recreation_percent_change_from_baseline,
            Grocery = grocery_and_pharmacy_percent_change_from_baseline,
            Parks = parks_percent_change_from_baseline,
            Transit = transit_stations_percent_change_from_baseline,
            Workplaces = workplaces_percent_change_from_baseline,
            Residential = residential_percent_change_from_baseline)
nat <- owid %>% full_join(ox, by = "date") %>% full_join(mob, by = "date") %>%
  filter(date >= ymd("2020-03-01"), date <= ymd("2022-10-15")) %>% arrange(date)
stopifnot(!anyDuplicated(nat$date))
write_csv(nat, P("data/processed/national_daily.csv"))

SURVEY <- c(ymd("2021-05-01"), ymd("2021-06-30"))   # mental-health survey window
PANEL_FU <- c(ymd("2020-11-01"), ymd("2020-11-30"))  # women's panel phone follow-up
LOCKDOWNS <- tribble(~start, ~end, ~label,
  ymd("2020-03-26"), ymd("2020-05-30"), "General holiday",
  ymd("2021-04-05"), ymd("2021-08-10"), "2021 restrictions")

# ---- Wave detection: PELT change points on weekly log cases --------------------
# Weekly means remove day-of-week reporting artefacts; daily data over-segment because of autocorrelation.
cp_dat <- nat %>% filter(!is.na(cases), cases > 0)
cp_wk <- cp_dat %>% mutate(week = floor_date(date, "week", week_start = 1)) %>%
  group_by(week) %>% summarise(cases = mean(cases)) %>% mutate(lc = log(cases))
cp <- cpt.meanvar(cp_wk$lc, method = "PELT", penalty = "MBIC", minseglen = 6)
cps <- cpts(cp)
# Sensitivity of the number of change points to the penalty and minimum segment length
# (PELT assumes independent errors, so the count is descriptive rather than a test result).
cp_sweep <- expand_grid(penalty = c("MBIC", "BIC", "SIC"), minseglen = c(4, 6, 8)) %>%
  mutate(n = map2_int(penalty, minseglen, ~ length(cpts(cpt.meanvar(cp_wk$lc, method = "PELT", penalty = .x, minseglen = .y)))))
cp_dates <- cp_wk$week[cps]
segs <- tibble(seg = seq_len(length(cps) + 1),
               start = cp_wk$week[c(1, cps + 1)], end = cp_wk$week[c(cps, nrow(cp_wk))] + 6) %>%
  rowwise() %>%
  mutate(mean_cases = mean(cp_dat$cases[cp_dat$date >= start & cp_dat$date <= end]),
         peak_cases = max(cp_dat$cases[cp_dat$date >= start & cp_dat$date <= end]),
         days = as.integer(end - start) + 1) %>% ungroup()

# Waves are labelled descriptively: the four largest peaks of smoothed cases at least 90 days apart.
# The number four is a choice, not an estimate.
pk <- cp_dat %>% arrange(desc(cases))
waves <- tibble()
for (i in seq_len(nrow(pk))) {
  if (nrow(waves) == 0 || all(abs(as.numeric(pk$date[i] - waves$date)) > 90)) waves <- bind_rows(waves, pk[i, ])
  if (nrow(waves) == 4) break
}
waves <- waves %>% arrange(date) %>% mutate(wave = paste("Wave", row_number()))
death_peak <- nat %>% filter(deaths == max(deaths, na.rm = TRUE)) %>% slice(1)

# ---- Figure: national timeline -------------------------------------------------
shade <- function() list(
  annotate("rect", xmin = SURVEY[1], xmax = SURVEY[2], ymin = -Inf, ymax = Inf, fill = PAL["orange"], alpha = .25),
  annotate("rect", xmin = PANEL_FU[1], xmax = PANEL_FU[2], ymin = -Inf, ymax = Inf, fill = PAL["green"], alpha = .25))
xlim_all <- range(nat$date)
p_cases <- ggplot(nat, aes(date, cases)) + shade() +
  geom_area(fill = PAL["blue"], alpha = .7) +
  geom_vline(xintercept = cp_dates, linetype = 3, colour = "grey30", linewidth = .3) +
  geom_point(data = waves, aes(date, cases), colour = PAL["red"]) +
  geom_text(data = waves, aes(date, cases, label = wave), vjust = -.6, size = 2.8) +
  scale_y_continuous(labels = comma, expand = expansion(c(0, .15))) +
  labs(x = NULL, y = "New cases / day", title = "A. Confirmed cases (7-day mean) with PELT change points")
p_pos <- ggplot(nat %>% filter(!is.na(positivity)), aes(date, positivity)) + shade() +
  geom_line(colour = PAL["red"]) + scale_y_continuous(labels = percent) +
  labs(x = NULL, y = "Test positivity", title = "B. Test positivity rate")
p_str <- ggplot(nat %>% filter(!is.na(stringency)), aes(date)) + shade() +
  geom_line(aes(y = stringency, colour = "Stringency index")) +
  geom_line(aes(y = econ_support, colour = "Economic support index")) +
  geom_line(data = nat %>% filter(!is.na(fully_vax)), aes(y = fully_vax, colour = "Fully vaccinated (%)")) +
  scale_colour_manual(values = unname(PAL[c("green", "purple", "blue")]), name = NULL) +
  labs(x = NULL, y = "Index / %", title = "C. Government response (OxCGRT) and vaccination")
p_mob <- nat %>% select(date, Workplaces, Residential, Transit) %>%
  pivot_longer(-date) %>% filter(!is.na(value)) %>%
  group_by(name) %>% mutate(value = zoo::rollmean(value, 7, fill = NA)) %>%
  ggplot(aes(date, value, colour = name)) + shade() + geom_hline(yintercept = 0, colour = "grey50") +
  geom_line() + scale_colour_manual(values = unname(PAL[c("grey", "orange", "blue")]), name = NULL) +
  labs(x = NULL, y = "% vs. baseline", title = "D. Google community mobility (7-day mean)")
fig_timeline <- (p_cases / p_pos / p_str / p_mob) &
  scale_x_date(limits = xlim_all, date_breaks = "3 months", date_labels = "%b\n%Y") &
  theme(legend.position = "right", plot.title = element_text(size = 9.5))
fig_timeline <- fig_timeline + plot_annotation(
  caption = "Orange band: mental-health survey (May–Jun 2021). Green band: women's panel follow-up (Nov 2020). Dotted lines: PELT change points.")
save_fig(fig_timeline, "fig_national_timeline", width = 7, height = 9)

# ---- Mobility PCA --------------------------------------------------------------
mob_cc <- nat %>% select(date, Retail:Residential) %>% drop_na()
pca <- prcomp(mob_cc %>% select(-date), scale. = TRUE)
pve <- summary(pca)$importance[2, ]
load <- as_tibble(pca$rotation[, 1:2], rownames = "Category")
if (load$PC1[load$Category == "Residential"] < 0) { pca$rotation[, 1] <- -pca$rotation[, 1]; pca$x[, 1] <- -pca$x[, 1] }
# orient PC1 so that higher = more time at home ("retreat index")
load <- as_tibble(pca$rotation[, 1:2], rownames = "Category")
mob_cc$retreat <- pca$x[, 1]
save_tab(load %>% mutate(across(where(is.numeric), ~round(., 3))) %>%
           rename(`PC1 (retreat)` = PC1, PC2 = PC2),
         "tab_mobility_pca",
         caption = sprintf("Loadings of the first two principal components of the six Google mobility series (PC1 explains %s\\%% and PC2 %s\\%% of variance).", pct(pve[1]), pct(pve[2])),
         digits = 3, escape = FALSE)

# ---- Policy -> behaviour: does stringency move people less over time? ---------
wk <- nat %>% left_join(mob_cc %>% select(date, retreat), by = "date") %>%
  mutate(week = floor_date(date, "week", week_start = 1)) %>%
  group_by(week) %>%
  summarise(cases = mean(cases, na.rm = TRUE), stringency = mean(stringency, na.rm = TRUE),
            retreat = mean(retreat, na.rm = TRUE), positivity = mean(positivity, na.rm = TRUE)) %>%
  filter(is.finite(retreat), is.finite(stringency), cases > 0) %>%
  mutate(period = factor(case_when(week < ymd("2021-01-01") ~ "2020",
                                   week < ymd("2022-01-01") ~ "2021", TRUE ~ "2022")),
         growth = log(lead(cases)) - log(cases))           # next-week log growth of cases
# Levels model (descriptive only): residuals are strongly autocorrelated, so inference uses first differences.
m_lev <- lm(retreat ~ stringency * period, data = wk)
lev_rho <- cor(head(resid(m_lev), -1), tail(resid(m_lev), -1))
lev_dw <- sum(diff(resid(m_lev))^2) / sum(resid(m_lev)^2)
# Primary model: week-on-week change in retreat on week-on-week change in stringency, slope by year.
wkd <- wk %>% mutate(d_retreat = retreat - lag(retreat), d_str = stringency - lag(stringency)) %>% filter(!is.na(d_retreat))
m_pol <- lm(d_retreat ~ d_str * period, data = wkd)
ct_pol <- coeftest(m_pol, vcov = NeweyWest(m_pol, lag = 4, prewhite = FALSE))
dif_rho <- cor(head(resid(m_pol), -1), tail(resid(m_pol), -1))
slopes <- tibble(period = levels(wkd$period)) %>%
  mutate(slope = map_dbl(period, ~ { b <- coef(m_pol); b["d_str"] + ifelse(.x == "2020", 0, b[paste0("d_str:period", .x)]) }))
# Robustness: drop March 2020, when the first lockdown moved both series at once.
m_pol_x <- lm(d_retreat ~ d_str * period, data = filter(wkd, week >= ymd("2020-04-01")))
ct_pol_x <- coeftest(m_pol_x, vcov = NeweyWest(m_pol_x, lag = 4, prewhite = FALSE))
fig_pol <- ggplot(wk, aes(stringency, retreat, colour = period)) +
  geom_point(alpha = .6, size = 1.4) + geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  scale_colour_manual(values = unname(PAL[c("blue", "orange", "green")]), name = "Year") +
  labs(x = "OxCGRT stringency index (weekly mean)", y = "Mobility retreat index (PC1)",
       title = "Policy stringency vs. behavioural retreat, by year")
save_fig(fig_pol, "fig_policy_mobility", width = 6.3, height = 3.6)

tab_pol <- tidy(ct_pol) %>%
  mutate(term = recode(term, "(Intercept)" = "Intercept", d_str = "$\\Delta$ Stringency",
                       period2021 = "Year 2021", period2022 = "Year 2022",
                       `d_str:period2021` = "$\\Delta$ Stringency × 2021", `d_str:period2022` = "$\\Delta$ Stringency × 2022"),
         p.value = fp(p.value)) %>%
  rename(Term = term, Estimate = estimate, `SE (HAC)` = std.error, t = statistic, p = p.value)
save_tab(tab_pol, "tab_policy_mobility",
         caption = "Week-on-week change in the mobility retreat index regressed on week-on-week change in policy stringency, with the slope allowed to differ by year (Newey–West HAC standard errors, 4 lags).",
         digits = 3, escape = FALSE,
         notes = sprintf("n = %d weekly differences; $R^2$ = %s. Residual lag-1 autocorrelation is %s, against %s in the same model fitted in levels (Durbin–Watson %s), which is why the model is fitted in differences.",
                         nrow(wkd), f3(summary(m_pol)$r.squared), f2(dif_rho), f2(lev_rho), f2(lev_dw)))

# ---- Behaviour -> transmission: lagged mobility vs. case growth ---------------
lag_fit <- map_dfr(0:4, function(L) {
  d <- wk %>% mutate(x = lag(retreat, L), s = lag(stringency, L)) %>% filter(!is.na(growth), !is.na(x))
  m <- lm(growth ~ x, data = d)
  ct <- coeftest(m, vcov = NeweyWest(m, lag = 4, prewhite = FALSE))
  tibble(lag = L, beta = ct["x", 1], se = ct["x", 2], p = ct["x", 4], r2 = summary(m)$r.squared, n = nrow(d))
})
best_lag <- lag_fit %>% slice_min(p, n = 1)

# ---- Survey-window context numbers --------------------------------------------
sw <- nat %>% filter(date >= SURVEY[1], date <= SURVEY[2])
pf <- nat %>% filter(date >= PANEL_FU[1], date <= PANEL_FU[2])
last <- nat %>% filter(!is.na(total_cases)) %>% slice_tail(n = 1)
macro("NatTotalCases", fint(last$total_cases)); macro("NatTotalDeaths", fint(last$total_deaths))
macro("NatEndDate", format(last$date, "%d %B %Y"))
macro("NatCFR", f2(100 * last$total_deaths / last$total_cases))
for (i in seq_len(nrow(waves))) {
  nm <- c("One", "Two", "Three", "Four")[i]
  macro(paste0("Wave", nm, "Date"), format(waves$date[i], "%d %B %Y"))
  macro(paste0("Wave", nm, "Cases"), fint(waves$cases[i]))
}
macro("NWaves", nrow(waves))
macro("DeathPeakDate", format(death_peak$date, "%d %B %Y")); macro("DeathPeak", fint(death_peak$deaths))
macro("NChangepoints", length(cps))
macro("NChangepointsLo", min(cp_sweep$n)); macro("NChangepointsHi", max(cp_sweep$n))
macro("SurveyStringency", f1(mean(sw$stringency, na.rm = TRUE)))
macro("SurveyPositivity", f1(100 * mean(sw$positivity, na.rm = TRUE)))
macro("SurveyCases", fint(mean(sw$cases, na.rm = TRUE)))
macro("SurveyCasesStart", fint(sw$cases[1])); macro("SurveyCasesMin", fint(min(sw$cases, na.rm = TRUE)))
macro("SurveyCasesEnd", fint(tail(sw$cases, 1)))
macro("SurveyPosStart", f1(100 * sw$positivity[1])); macro("SurveyPosEnd", f1(100 * tail(na.omit(sw$positivity), 1)))
macro("SurveyVaxEnd", f1(max(sw$fully_vax, na.rm = TRUE)))
macro("SurveyResidential", f1(mean(sw$Residential, na.rm = TRUE)))
macro("SurveyWorkplaces", f1(mean(sw$Workplaces, na.rm = TRUE)))
macro("PanelStringency", f1(mean(pf$stringency, na.rm = TRUE)))
macro("PanelPositivity", f1(100 * mean(pf$positivity, na.rm = TRUE)))
macro("MaxStringency", f1(max(nat$stringency, na.rm = TRUE)))
macro("MaxStringencyDate", format(nat$date[which.max(nat$stringency)], "%d %B %Y"))
macro("MinWorkplaces", f1(min(zoo::rollmean(mob_cc$Workplaces, 7))))
macro("PCAOne", pct(pve[1])); macro("PCATwo", pct(pve[2]))
macro("PolSlopeA", f3(slopes$slope[1])); macro("PolSlopeB", f3(slopes$slope[2])); macro("PolSlopeC", f3(slopes$slope[3]))
macro("PolRsquared", f3(summary(m_pol)$r.squared)); macro("PolNweeks", nrow(wkd))
macro("PolRatio", f1(slopes$slope[1] / slopes$slope[3]))
macro("PolLevRho", f2(lev_rho)); macro("PolLevDw", f2(lev_dw)); macro("PolDifRho", f2(dif_rho))
macro("PolIntCpNoMarch", fp(ct_pol_x["d_str:period2022", 4]))
macro("PolLevSlopeA", f3(coef(m_lev)["stringency"])); macro("PolLevSlopeC", f3(coef(m_lev)["stringency"] + coef(m_lev)["stringency:period2022"]))
macro("PolIntBp", fp(ct_pol["d_str:period2021", 4])); macro("PolIntCp", fp(ct_pol["d_str:period2022", 4]))
macro("LagBest", best_lag$lag); macro("LagBeta", f3(best_lag$beta)); macro("LagP", fp(best_lag$p))
macro("LagZeroBeta", f3(lag_fit$beta[1])); macro("LagZeroP", fp(lag_fit$p[1]))
macro("LagTwoBeta", f3(lag_fit$beta[3])); macro("LagTwoP", fp(lag_fit$p[3]))
long_seg <- segs %>% slice_max(days, n = 1)
macro("RegimeLongStart", format(long_seg$start, "%d %B %Y")); macro("RegimeLongEnd", format(long_seg$end, "%d %B %Y"))
between <- nat %>% filter(date > waves$date[2], date < waves$date[3], !is.na(positivity))
macro("PosTrough", f1(100 * min(between$positivity)))
macro("PosPeakThree", f1(100 * max(nat$positivity[nat$date > waves$date[2] & nat$date < waves$date[4]], na.rm = TRUE)))
macro("PosLowEarly", f1(100 * min(nat$positivity[nat$date >= ymd("2021-01-01") & nat$date < waves$date[2]], na.rm = TRUE)))
macro("PosLastDate", format(max(nat$date[!is.na(nat$positivity)]), "%B %Y"))
write_macros("national")

save_tab(segs %>% mutate(start = format(start, "%d %b %Y"), end = format(end, "%d %b %Y"),
                         mean_cases = fint(mean_cases), peak_cases = fint(peak_cases)) %>%
           rename(Segment = seg, Start = start, End = end, `Mean cases/day` = mean_cases,
                  `Peak cases/day` = peak_cases, Days = days),
         "tab_changepoints", digits = 0,
         caption = "Descriptive segmentation of the epidemic curve by PELT (changes in mean and variance of log weekly-mean cases, MBIC penalty, minimum segment six weeks).")
saveRDS(list(lag_fit = lag_fit, slopes = slopes, waves = waves, segs = segs), P("data/processed/national_results.rds"))
cat("01_national done: ", length(cps), "change points;", nrow(wk), "weeks\n")
print(lag_fit); print(slopes); print(cp_sweep); print(ct_pol); print(ct_pol_x); print(waves %>% select(date, cases, wave))
