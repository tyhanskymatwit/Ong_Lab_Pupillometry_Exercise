#Package setup
suppressPackageStartupMessages({ library(survival); library(dplyr) })

#Data output directory setup
args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "out"
pat <- readRDS(file.path(out_dir, "patients_clean.rds"))
pup <- readRDS(file.path(out_dir, "pupil_clean.rds"))

# Splits each patient's follow-up at every exam, carrying forward the latest valid pupil values (AI USED TO CORRECT SYNTAX)
build_cp <- function(lag_h = 0) {
  pup %>% mutate(avail_t = hours_from_admit + lag_h) %>%
    filter(avail_t < time) %>%
    group_by(patient_id) %>% arrange(avail_t, .by_group = TRUE) %>%
    mutate(lv = cummax(ifelse(is.na(npi_min), 0L, row_number())),
           npi_min = ifelse(lv > 0, npi_min[pmax(lv, 1)], NA), npi_asym = ifelse(lv > 0, npi_asym[pmax(lv, 1)], NA)) %>%
    filter(lv > 0) %>%
    mutate(tstart = avail_t, tstop = c(avail_t[-1], time[1]),
           status = as.integer(event == 1 & row_number() == n())) %>%
    ungroup() %>% filter(tstop > tstart) %>%
    select(patient_id, tstart, tstop, status, npi_min, npi_asym, last_failed = failed) %>%
    left_join(pat %>% select(patient_id, age, sex, diagnosis, admit_gcs, admit_midline_shift_mm,
                             sedation_infusion_day1, hx_cataract_surgery), by = "patient_id") %>%
    mutate(npi_low = -npi_min, age10 = age / 10)
}

# Builds the interval data with and without a 2 h lag and reports events left out
cp <- build_cp(0); cp2 <- build_cp(2)
cat(sprintf("Counting-process data: %d intervals, %d patients, %d events (of %d total)\n",
            nrow(cp), n_distinct(cp$patient_id), sum(cp$status), sum(pat$event)))
lost <- pat %>% filter(event == 1, !patient_id %in% cp$patient_id[cp$status == 1])
cat("Events not in pupil risk set (no exam or event before first valid exam):", paste(lost$patient_id, collapse = ", "), "\n")

# Fits admission-only and admission + pupillometry Cox models, with and without the 2 h lag
adm_terms <- "age10 + sex + diagnosis + admit_gcs + admit_midline_shift_mm + sedation_infusion_day1 + hx_cataract_surgery"
fit <- function(extra = "", d = cp) coxph(as.formula(paste("Surv(tstart, tstop, status) ~", adm_terms, extra)), data = d)
m_adm <- fit(); m_pup <- fit("+ npi_low + npi_asym + last_failed")
m_adm2 <- fit(d = cp2); m_pup2 <- fit("+ npi_low + npi_asym + last_failed", cp2)

# Extracts hazard ratios with 95% CIs and labels each variable
hr_tab <- function(m, model) {
  s <- summary(m)
  data.frame(model = model, term = rownames(s$conf.int), HR = s$conf.int[, 1], lo = s$conf.int[, 3], hi = s$conf.int[, 4],
             p = s$coefficients[, "Pr(>|z|)"], row.names = NULL)
}
lab <- c(age10 = "Age (per 10 y)", sexM = "Male sex", diagnosisICH = "ICH vs AIS",
         diagnosisSAH = "SAH vs AIS", diagnosisTBI = "TBI vs AIS",
         admit_gcs = "Admission GCS (per point)", admit_midline_shift_mm = "Midline shift (per mm)",
         sedation_infusion_day1 = "Sedative infusion day 1", hx_cataract_surgery = "Prior cataract surgery",
         npi_low = "Worst-eye NPi (per 1-point lower)", npi_asym = "NPi asymmetry L-R (per 1 point)",
         last_failed = "Last pupil exam could not be done")

# Saves and prints the hazard ratio table for each model
hrs <- bind_rows(hr_tab(m_adm,  "A: admission only"),
                 hr_tab(m_pup,  "B: admission + pupillometry (primary)"),
                 hr_tab(m_pup2, "S: as B, pupil data lagged 2 h")) %>%
  mutate(label = lab[term], `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi), p = signif(p, 2))
write.csv(hrs, file.path(out_dir, "cox_hazard_ratios.csv"), row.names = FALSE)
for (mm in unique(hrs$model)) { cat("\n==", mm, "==\n"); print(hrs %>% filter(model == mm) %>% select(label, `HR (95% CI)`, p), row.names = FALSE) }

# Tests whether pupillometry adds to admission data (likelihood ratio, AIC, C-index)
cmp <- function(a, b, tag) {
  lr <- anova(a, b)
  cd <- concordance(a, b)
  dv <- cd$concordance[2] - cd$concordance[1]; se <- sqrt(c(1, -1) %*% cd$var %*% c(1, -1))
  data.frame(comparison = tag, LR_chisq = lr$Chisq[2], df = lr$Df[2], LR_p = lr$`Pr(>|Chi|)`[2],
             AIC_adm = AIC(a), AIC_pup = AIC(b), C_adm = cd$concordance[1], C_pup = cd$concordance[2],
             C_diff = dv, C_diff_lo = dv - 1.96 * se, C_diff_hi = dv + 1.96 * se, row.names = NULL)
}
added <- bind_rows(cmp(m_adm, m_pup, "Primary (no lag)"), cmp(m_adm2, m_pup2, "Pupil data lagged 2 h"))
cat("\n== Added value of pupillometry ==\n"); print(added, digits = 3, row.names = FALSE)
write.csv(added, file.path(out_dir, "cox_added_value.csv"), row.names = FALSE)

# Checks the proportional hazards assumption for the pupillometry model
cat("\n== Proportional hazards check (model B) ==\n"); print(cox.zph(m_pup))