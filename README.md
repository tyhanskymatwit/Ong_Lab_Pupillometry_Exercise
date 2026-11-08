# Pupillometry and neurologic deterioration in the Neuro ICU

This is the take home assessment provided by Dr. Ong and her lab to assess skills in survival analysis, data cleaning, and clinical reporting. It took the better part of 2 days to compile the code and report in a coherent and effective format after a lot of research.

## Required R Packages

```r
install.packages(c("dplyr", "survival", "ranger", "pROC", "ggplot2", "sandwich", "lmtest"))
```
dplyr(1.1.4): Data manipulation
survival(3.5.8): Cox models with time-varying covariates, PH check, C-index 
ranger(0.16.0): Random forest prediction model 
PROC(1.18.5): AUROC 
ggplot2(3.4.4): Clinician figure 
sandwich(3.1.0), lmtest(0.9.40): Patient-clustered standard errors in the missingness model 


##CSVs Exported

cleaning_log.csv: Every cleaning decision, the number of rows affected and the reason
missingness_model.csv: Odds ratios for a failed exam (is the missingness random?)
cox_hazard_ratios.csv: Hazard ratios with 95% CIs: admission-only, admission + pupillometry, and the 2-h lag sensitivity model
cox_added_value.csv: Likelihood-ratio test, AIC and C-index: admission-only vs + pupillometry
prediction_performance.csv: Test-set AUROC, AUPRC, calibration intercept and slope, and Brier score, each with bootstrap 95% CIs
prediction_calibration_deciles.csv: Predicted vs observed risk by decile
figure_npi_24h_risk.png: The clinician figure
