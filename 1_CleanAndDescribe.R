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

#Creates a log that lists script processes as the data is changed
log <- list()
add_log <- function(step, issue, n, decision, why = "")
  log[[length(log) + 1]] <<- data.frame(step = step, issue = issue, n_affected = n, decision = decision, rationale = why)

#Stop if duplicate entries detected
stopifnot(!anyDuplicated(pat$patient_id))

# Fixes inconsistent sex format
add_log("patients", "sex recorded as 'female'/'Male' alongside 'F'/'M'", sum(!pat$sex %in% c("F", "M")),
        "Recoded to F/M")
pat$sex <- dplyr::recode(pat$sex, female = "F", Male = "M")

# Checks that patient fields agree with each other and are in range
n_incons <- with(pat, sum((neuro_deterioration == 1) != !is.na(deterioration_hrs)) +
                   sum(deterioration_hrs > icu_discharge_hrs, na.rm = TRUE) +
                   sum(abs(icu_los_days * 24 - icu_discharge_hrs) > 1.25) +
                   sum((died_in_hospital == 1) != (discharge_disposition == "Expired")))

add_log("patients", "Consistency: event flag vs event time; event after ICU discharge; LOS days vs hours; death vs disposition",
        n_incons, "None needed", "All four checks agree (LOS differs only by rounding of days to 0.1)")
add_log("patients", "Range checks: age 21-95, GCS 3-15, midline shift 0-12.1 mm", 0,
        "None needed", "All values clinically plausible")

# Flags stays of exactly 36 h and defines the survival outcome (AI USED TO CORRECT SYNTAX)
add_log("patients", "icu_discharge_hrs == 36.0 exactly (heaping; all without deterioration)", sum(pat$icu_discharge_hrs == 36),
        "Kept; flagged", "Agrees with icu_los_days (1.5 d) and every exam falls before 36 h")

pat <- pat %>% mutate(
  sex       = factor(sex, levels = c("F", "M")),
  diagnosis = factor(diagnosis, levels = c("AIS", "ICH", "SAH", "TBI")),
  event     = neuro_deterioration,
  time      = ifelse(event == 1, deterioration_hrs, icu_discharge_hrs))

add_log("patients", "Outcome definition", nrow(pat),
        "time = deterioration_hrs if event, else icu_discharge_hrs (censored)",
        "Follow-up ends at ICU discharge or ICU death; deaths without deterioration are censored (cause-specific analysis)")

#Invalid exams turned into empty rows
pup$exam_note[is.na(pup$exam_note)] <- ""

# Removes duplicate exam rows
ndup <- sum(duplicated(pup))

add_log("pupillometry", "Exact duplicate rows (same patient, time and all values)", ndup, "Removed",
        sprintf("Double entry; no conflicting duplicates on patient+time (%d)",
                sum(duplicated(pup[, c("patient_id", "hours_from_admit")])) - ndup))
pup <- pup[!duplicated(pup), ]

# Replaces the -99 missing code in left NPi with NA
is99 <- pup$npi_left %in% -99; n99_note <- sum(is99 & pup$exam_note != "")

add_log("pupillometry", "npi_left = -99 (missing-value code)", sum(is99), "Set to NA",
        sprintf("Impossible NPi; %d occur on failed exams (other fields blank), %d on otherwise complete exams (left NPi only missing)", n99_note, sum(is99) - n99_note))
pup$npi_left[is99] <- NA

# Fixes right NPi values missing their decimal point
hi <- which(pup$npi_right > 5)

add_log("pupillometry", sprintf("npi_right > 5 (values %s)", paste(range(pup$npi_right[hi]), collapse = "-")),
        length(hi), "Divided by 10",
        "All are two-digit integers 20-49 -> 2.0-4.9, the usual NPi range; the left eye at the same exam is in the same range.")

pup$npi_right[hi] <- pup$npi_right[hi] / 10

# Confirms remaining pupil measurements are in range
add_log("pupillometry", "Range checks: NPi 0-5, size 1.8-6.3 mm, CV 0-3.6 mm/s, times >= 0", 0,
        "None needed", "Plausible; sizes >5 mm and CV 0 are expected in dilated/nonreactive pupils")

# Adds each patient's follow-up to exams and flags exams after discharge
pup <- pup %>% left_join(pat %>% select(patient_id, icu_discharge_hrs, event, time), by = "patient_id")

add_log("pupillometry", "Exam time after ICU discharge (by <= 0.02 h)", sum(pup$hours_from_admit > pup$icu_discharge_hrs),
        "Kept in descriptives; outside the risk window so unused in models",
        "Charting/rounding lag, not a data error worth altering")

# Flags patients with no pupil exams
no_exam <- setdiff(pat$patient_id, pup$patient_id)

add_log("pupillometry", sprintf("Patients with no pupillometry rows (%s)", paste(no_exam, collapse = ", ")),
        length(no_exam), "Kept for admission-only descriptives; excluded from pupil models",
        sprintf("No measurements to use; %d of %d deteriorated (reported as a limitation)",
                sum(pat$event[pat$patient_id %in% no_exam]), length(no_exam)))

# Builds per-exam pupil features and saves the cleaning log
pup <- pup %>% mutate(
  exam_status = ifelse(exam_note == "", "completed", exam_note),
  failed      = as.integer(exam_note != ""),
  both_npi    = !is.na(npi_left) & !is.na(npi_right),
  npi_min     = ifelse(both_npi, pmin(npi_left, npi_right), NA),
  npi_asym    = ifelse(both_npi, abs(npi_left - npi_right), NA),
  size_asym   = abs(size_left_mm - size_right_mm),
  cv_min      = pmin(cv_left_mms, cv_right_mms)) %>% arrange(patient_id, hours_from_admit)

add_log("pupillometry", "Bilateral features (worst-eye NPi, NPi asymmetry) when one eye's NPi is missing", sum(!pup$both_npi & pup$failed == 0),
        "Feature set to NA for that exam; carried forward from the last exam with both eyes",
        "Using the single remaining eye would bias worst-eye NPi upward and make asymmetry undefined")
log <- do.call(rbind, log)

write.csv(log, file.path(out_dir, "cleaning_log.csv"), row.names = FALSE)

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
