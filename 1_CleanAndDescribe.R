#Package installation/setup
install.packages(c("dplyr", "survival", "ranger", "pROC", "ggplot2", "sandwich", "lmtest"))
suppressPackageStartupMessages({ library(dplyr); library(sandwich); library(lmtest) })
args     <- commandArgs(trailingOnly = TRUE)

#Data output directory setup
data_dir <- if (length(args) >= 1) args[1] else "."
out_dir  <- if (length(args) >= 2) args[2] else "out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.csv(file.path(data_dir, f), stringsAsFactors = FALSE, na.strings = c("", "NA"))
pat <- rd("patients.csv"); pup <- rd("pupillometry.csv")

#Stop if duplicate entries detected
stopifnot(!anyDuplicated(pat$patient_id))

# Fixes inconsistent sex format
pat$sex <- dplyr::recode(pat$sex, female = "F", Male = "M")

# Defines the survival outcome (AI USED TO CORRECT SYNTAX)
pat <- pat %>% mutate(
  sex       = factor(sex, levels = c("F", "M")),
  diagnosis = factor(diagnosis, levels = c("AIS", "ICH", "SAH", "TBI")),
  event     = neuro_deterioration,
  time      = ifelse(event == 1, deterioration_hrs, icu_discharge_hrs))

#Invalid exams turned into empty rows
pup$exam_note[is.na(pup$exam_note)] <- ""

# Removes duplicate exam rows
pup <- pup[!duplicated(pup), ]

# Replaces the -99 missing code in left NPi with NA
is99 <- pup$npi_left %in% -99
pup$npi_left[is99] <- NA

# Fixes right NPi values missing their decimal point
hi <- which(pup$npi_right > 5)
pup$npi_right[hi] <- pup$npi_right[hi] / 10

# Adds each patient's follow-up to exams
pup <- pup %>% left_join(pat %>% select(patient_id, icu_discharge_hrs, event, time), by = "patient_id")

# Builds per-exam pupil features
pup <- pup %>% mutate(
  exam_status = ifelse(exam_note == "", "completed", exam_note),
  failed      = as.integer(exam_note != ""),
  both_npi    = !is.na(npi_left) & !is.na(npi_right),
  npi_min     = ifelse(both_npi, pmin(npi_left, npi_right), NA),
  npi_asym    = ifelse(both_npi, abs(npi_left - npi_right), NA),
  size_asym   = abs(size_left_mm - size_right_mm),
  cv_min      = pmin(cv_left_mms, cv_right_mms)) %>% arrange(patient_id, hours_from_admit)

#Output of cohorts
cat(sprintf("Patients: %d; deteriorated: %d (%.1f%%); median time to deterioration %.1f h (IQR %.1f-%.1f)\n",
            nrow(pat), sum(pat$event), 100 * mean(pat$event), median(pat$deterioration_hrs, na.rm = TRUE),
            quantile(pat$deterioration_hrs, .25, na.rm = TRUE), quantile(pat$deterioration_hrs, .75, na.rm = TRUE)))
print(pat %>% group_by(deteriorated = event) %>% summarise(n = n(), age = median(age), gcs = median(admit_gcs),
      midline_shift = median(admit_midline_shift_mm), sedation_pct = round(100 * mean(sedation_infusion_day1))))

#Output of missingness
pup_risk <- pup %>% filter(hours_from_admit <= icu_discharge_hrs) %>%
  mutate(pre_event_12h = as.integer(event == 1 & (time - hours_from_admit) <= 12),
         day1 = as.integer(hours_from_admit <= 24))
cat(sprintf("Exams (after de-dup): %d in %d patients; median %.0f exams/patient; median gap %.2f h\n",
            nrow(pup), n_distinct(pup$patient_id), median(table(pup$patient_id)),
            median(ave(pup$hours_from_admit, pup$patient_id, FUN = function(x) c(NA, diff(x))), na.rm = TRUE)))
miss_reason <- pup %>% count(exam_status) %>% mutate(pct = round(100 * n / sum(n), 1))
print(miss_reason)
cat(sprintf("Any pupil measurement missing: %.1f%% of exams; whole exam failed: %.1f%%; left NPi only missing: %d exams\n",
            100 * mean(is.na(pup$npi_min)), 100 * mean(pup$failed), sum(!pup$both_npi & pup$failed == 0)))
pfail <- tapply(pup$failed, pup$patient_id, mean)
cat(sprintf("Per patient failed-exam share: median %.1f%%, IQR %.1f-%.1f%%, max %.1f%%; %d%% of patients had >=1 failed exam\n",
            100 * median(pfail), 100 * quantile(pfail, .25), 100 * quantile(pfail, .75), 100 * max(pfail), round(100 * mean(pfail > 0))))

# Is missingness random? Exam-level logistic regression of exam failure on
# admission variables, time since admission and proximity to deterioration,
# with patient-clustered (sandwich) standard errors.
md <- pup_risk %>% left_join(pat %>% select(patient_id, age, sex, diagnosis, admit_gcs, admit_midline_shift_mm,
                                            sedation_infusion_day1, hx_cataract_surgery), by = "patient_id")
fit_m <- glm(failed ~ age + sex + diagnosis + admit_gcs + admit_midline_shift_mm +
               sedation_infusion_day1 + hx_cataract_surgery + day1 + pre_event_12h, data = md, family = binomial)
ct <- coeftest(fit_m, vcov. = vcovCL(fit_m, cluster = ~patient_id))
or <- data.frame(term = rownames(ct), OR = exp(ct[, 1]), lo = exp(ct[, 1] - 1.96 * ct[, 2]),
                 hi = exp(ct[, 1] + 1.96 * ct[, 2]), p = ct[, 4])[-1, ]
cat("\nOdds of a failed exam (patient-clustered 95% CI):\n"); print(or, digits = 3, row.names = FALSE)
write.csv(or, file.path(out_dir, "missingness_model.csv"), row.names = FALSE)

print(md %>% group_by(within_12h_before_event = pre_event_12h) %>% summarise(n = n(), fail_pct = round(100 * mean(failed), 1)))

saveRDS(pat, file.path(out_dir, "patients_clean.rds"))
saveRDS(pup, file.path(out_dir, "pupil_clean.rds"))
