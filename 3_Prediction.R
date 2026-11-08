#Package setup
suppressPackageStartupMessages({ library(dplyr); library(ranger); library(pROC) })

#Data output directory setup
args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "out"
set.seed(20261107)
pat <- readRDS(file.path(out_dir, "patients_clean.rds"))
pup <- readRDS(file.path(out_dir, "pupil_clean.rds"))

# Sets the 24 h prediction window and a helper for saving results
H <- 24
wcsv <- function(x, f) write.csv(x, file.path(out_dir, f), row.names = FALSE)

# Helpers that carry values forward, look back in time and summarise recent exams
locf <- function(x) { i <- cummax(ifelse(is.na(x), 0L, seq_along(x))); ifelse(i > 0, x[pmax(i, 1)], NA) }
lag_val <- function(t, x, w) { idx <- findInterval(t - w, t); ifelse(idx > 0, x[pmax(idx, 1)], NA) }
roll <- function(t, x, w, f) sapply(seq_along(t), function(i) {
  v <- x[t >= t[i] - w & t <= t[i]]; v <- v[!is.na(v)]; if (length(v)) f(v) else NA })

# Builds one prediction row per exam from admission data, current pupil values and recent trends
feat <- pup %>%
  filter(hours_from_admit < time) %>%
  group_by(patient_id) %>% arrange(hours_from_admit, .by_group = TRUE) %>%
  mutate(t = hours_from_admit, npi_min_c = locf(npi_min), npi_asym_c = locf(npi_asym),
         size_asym_c = locf(size_asym), cv_min_c = locf(cv_min),
         d_npi_6h = coalesce(npi_min_c - lag_val(t, npi_min_c, 6), 0), d_npi_12h = coalesce(npi_min_c - lag_val(t, npi_min_c, 12), 0),
         min_npi_12h = roll(t, npi_min, 12, min), max_asym_12h = coalesce(roll(t, npi_asym, 12, max), npi_asym_c),
         fail_frac_12h = roll(t, failed, 12, mean), hrs_since_valid = t - locf(ifelse(is.na(npi_min), NA, t))) %>%
  ungroup() %>%
  filter(!is.na(npi_min_c)) %>%
  left_join(pat %>% select(patient_id, age, sex, diagnosis, admit_gcs, admit_midline_shift_mm,
                           sedation_infusion_day1, hx_cataract_surgery), by = "patient_id") %>%
  mutate(y = factor(as.integer(event == 1 & time - t <= H), levels = 0:1, labels = c("no", "yes")),
         size_asym_c = coalesce(size_asym_c, 0), cv_min_c = coalesce(cv_min_c, median(cv_min_c, na.rm = TRUE)))

# Lists the model features and checks none are missing
x_all <- c("age", "sex", "diagnosis", "admit_gcs", "admit_midline_shift_mm", "sedation_infusion_day1",
           "hx_cataract_surgery", "t",
           "npi_min_c", "npi_asym_c", "size_asym_c", "cv_min_c", "d_npi_6h", "d_npi_12h",
           "min_npi_12h", "max_asym_12h", "failed", "fail_frac_12h", "hrs_since_valid")
stopifnot(!anyNA(feat[, x_all]))
cat(sprintf("Prediction rows: %d exams from %d patients; %.1f%% labelled 'deterioration within 24 h'\n",
            nrow(feat), n_distinct(feat$patient_id), 100 * mean(feat$y == "yes")))

# Splits patients 70/30 into training and test sets
pts <- feat %>% distinct(patient_id, event)
test_ids <- unlist(lapply(split(pts$patient_id, pts$event), function(v) sample(v, round(0.3 * length(v)))))
train <- feat %>% filter(!patient_id %in% test_ids); test <- feat %>% filter(patient_id %in% test_ids)
ev <- function(test) sum(pts$event[(pts$patient_id %in% test_ids) == test])
cat(sprintf("Train: %d patients (%d events), %d rows | Test: %d patients (%d events), %d rows\n",
            n_distinct(train$patient_id), ev(FALSE), nrow(train), n_distinct(test$patient_id), ev(TRUE), nrow(test)))

# Tunes the random forest with 5-fold cross-validation grouped by patient
tr_ids <- unique(train$patient_id); fold_of <- setNames(sample(rep(1:5, length.out = length(tr_ids))), tr_ids)
train$fold <- fold_of[train$patient_id]
grid <- expand.grid(mtry = c(3, 5, 8), min.node.size = c(25, 100, 300))
fit_rf <- function(d, mtry, mns, trees = 500)
  ranger(y ~ ., data = d[, c("y", x_all)], probability = TRUE, num.trees = trees, mtry = mtry,
         min.node.size = mns, num.threads = 4)
p_yes <- function(m, d) predict(m, d[, x_all])$predictions[, "yes"]
cv_oof <- lapply(seq_len(nrow(grid)), function(g) {
  oof <- numeric(nrow(train))
  for (k in 1:5) oof[train$fold == k] <- p_yes(fit_rf(train[train$fold != k, ], grid$mtry[g], grid$min.node.size[g]), train[train$fold == k, ])
  oof
})
grid$cv_auc <- sapply(cv_oof, function(o) as.numeric(auc(train$y, o, quiet = TRUE))); print(grid)
best <- which.max(grid$cv_auc); oof <- cv_oof[[best]]
cat(sprintf("Chosen: mtry=%d, min.node.size=%d (CV AUROC %.3f)\n", grid$mtry[best], grid$min.node.size[best], grid$cv_auc[best]))

# Recalibrates the forest's predictions using the out-of-fold results
lg <- function(p) qlogis(pmin(pmax(p, 1e-4), 1 - 1e-4))
recal <- glm(I(y == "yes") ~ lg(oof), data = data.frame(y = train$y, oof = oof), family = binomial)

# Fits the final forest and an admission-only logistic model for comparison
rf <- fit_rf(train, grid$mtry[best], grid$min.node.size[best], trees = 1000)
lr_adm <- glm(I(y == "yes") ~ age + sex + diagnosis + admit_gcs + admit_midline_shift_mm +
                sedation_infusion_day1 + hx_cataract_surgery + splines::ns(t, 3), data = train, family = binomial)

# Predicts 24 h risk for test-set exams with both models
test <- test %>% mutate(p_rf = predict(recal, newdata = data.frame(oof = p_yes(rf, test)), type = "response"),
                        p_adm = predict(lr_adm, newdata = test, type = "response"), yy = as.integer(y == "yes"))

# Measures discrimination and calibration on the test set with bootstrap 95% CIs (AI USED TO CORRECT SYNTAX)
auprc <- function(y, p) { y <- y[order(-p)]; sum((cumsum(y) / seq_along(y))[y == 1]) / sum(y) }
metrics <- function(y, p) {
  lp <- lg(p)
  c(auroc = as.numeric(auc(y, p, quiet = TRUE)), auprc = auprc(y, p),
    cal_intercept = unname(coef(glm(y ~ offset(lp), family = binomial))[1]),
    cal_slope = unname(coef(glm(y ~ lp, family = binomial))[2]), brier = mean((p - y)^2),
    mean_pred = mean(p), obs_rate = mean(y))
}
all_metrics <- function(d) {
  m1 <- metrics(d$yy, d$p_rf); m0 <- metrics(d$yy, d$p_adm)
  c(rf = m1, adm = m0, d_auroc = m1["auroc"] - m0["auroc"])
}
point <- all_metrics(test)
rows <- split(seq_len(nrow(test)), test$patient_id)
ci <- apply(replicate(1000, all_metrics(test[unlist(rows[sample(names(rows), replace = TRUE)]), ])),
            1, quantile, c(.025, .975), na.rm = TRUE)
perf <- data.frame(metric = names(point), estimate = point, lo = ci[1, ], hi = ci[2, ], row.names = NULL)
cat("\n== Test-set performance (patient-cluster bootstrap 95% CI) ==\n"); print(perf, digits = 3)
wcsv(perf, "prediction_performance.csv")

# Compares predicted and observed risk by decile
cal_tab <- test %>% mutate(bin = ntile(p_rf, 10)) %>% group_by(bin) %>%
  summarise(n = n(), mean_pred = mean(p_rf), obs = mean(yy),
            lo = binom.test(sum(yy), n())$conf.int[1], hi = binom.test(sum(yy), n())$conf.int[2])
cat("\n== Calibration by decile (recalibrated RF) ==\n"); print(cal_tab, digits = 3)
wcsv(cal_tab, "prediction_calibration_deciles.csv")

# Saves the prediction rows for the figure
saveRDS(feat, file.path(out_dir, "prediction_rows.rds"))
