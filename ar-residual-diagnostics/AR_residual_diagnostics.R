
#####################################
### AR RESIDUAL DIAGNOSTICS #########
#####################################

library(dplyr)
library(FKF)
library(tseries) 
library(tidyverse)

# set wd 
setwd("...")

# load data
data <- read.csv("data_FEEL_Study_1.csv")

# get helper function
source("get_ARWN_residuals.R")


# get all uuids
uuids <- unique(data$UUID) # 179 participants

# draw 1/3 of sample for testings
set.seed(2310)
uuids.sub <- sample(uuids, round(length(uuids) / 3)) # 60

# check if seeding worked
uuids.sub[3] # "5b4774af-4900-4d02-80ad-dc3a334e868f"
uuids.sub[23] # "4ea6afa1-e202-45eb-80c5-ad90f1ce6ce5"

# retrieve sub data
data.sub <- data[data$UUID %in% uuids.sub, ]

# get affect items
items.neg <- c("ANG_ES", "SAD_ES", "STR_ES")
items.pos <- c("CONF_ES", "HAP_ES", "RLX_ES") 
emotions <- c(items.pos, items.neg) 

## --------------------------------
# PROCEDURE -----------------------
## --------------------------------

## --------------------------------
# Setup parallel ------------------

n_cores <- availableCores() - 1
plan(multisession, workers = n_cores)

# Prepare data for parallel processing

task_data <- expand_grid(
  uuid = uuids.sub,
  emotion = emotions
) %>%
  mutate(
    ts = map2(uuid, emotion, function(u, e) {
      temp <- data.sub %>% filter(UUID == u)
      na.omit(temp[[e]])
    })
  )

## -------------------------------
# Helper function ----------------

fit_ar1wn_residuals_parallel <- function(ts, uuid, emotion) {
  tryCatch({
    # Frequentist
    freq_result <- fit_ar1wn_residuals(ts, verbose = FALSE)
    
    # Bayesian
    bayes_result <- fit_ar1wn_residuals_bayesian(ts, verbose = FALSE)
    
    # Return both
    list(
      uuid = uuid,
      emotion = emotion,
      frequentist = freq_result,
      bayesian = bayes_result,
      success = TRUE
    )
  }, error = function(e) {
    list(
      uuid = uuid,
      emotion = emotion,
      error = e$message,
      success = FALSE
    )
  })
}

## ------------------------------
# Run parallel loop -------------

start_time <- Sys.time()

# Single
results_parallel <- task_data %>%
  mutate(
    result = future_pmap(
      list(ts, uuid, emotion),
      fit_ar1wn_residuals_parallel,
      .progress = TRUE
    )
  )

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

# Save
saveRDS(results_parallel, file = "ARWN_results.rds")

# # Setup progress bar and output data frame
# total_iterations <- length(uuids.sub) * length(emotions)
# pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
# iteration <- 0
# 
# # Initialize empty result list
# results.ARWN <- list()
# results.ARWN.bayesian <- list()
# 
# # loop
# for (uuid in uuids.sub) {
#   temp <- data.sub %>% filter(UUID == uuid)
#   
#   for (e in emotions) {
#     iteration <- iteration + 1
#     setTxtProgressBar(pb, iteration)
#     
#     ts <- na.omit(temp[[e]])
#     
#     # --- Fit AR(1)+WN model from helper functions
#     
#     # Frequentist:
#     arwn.res <- fit_ar1wn_residuals(ts)
#     
#     # Bayesian:
#     arwn.res.bayes <- fit_ar1wn_residuals_bayesian(ts)
#     
#     results.ARWN[[uuid]][[e]] <- arwn.res
#     results.ARWN.bayesian[[uuid]][[e]] <- arwn.res.bayes
#   }
# }
# 
# close(pb)
# 
# ## --------------------------------
# # Check Heywood Cases -------------
# 
# # Check for Heywood cases across all time series
# heywood_summary <- data.frame()
# 
# for (uuid in names(results.ARWN)) {
#   for (emotion in names(results.ARWN[[uuid]])) {
#     res <- results.ARWN[[uuid]][[emotion]]
#     
#     heywood_summary <- rbind(heywood_summary, data.frame(
#       UUID = uuid,
#       Emotion = emotion,
#       convergence = res$convergence,
#       heywood_ivar = res$heywood_ivar,
#       heywood_evar = res$heywood_evar,
#       ivar = res$ivar,
#       evar = res$evar
#     ))
#   }
# }
# 
# 
# ## --------------------------------
# # BDS Test on Residuals -----------
# 
# # Calculate total iterations
# total_iterations <- sum(sapply(results.ARWN, function(x) length(x)))
# 
# # Setup progress bar
# pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
# iteration <- 0
# 
# # Initialize BDS results data frame with bootstrap
# bds_results <- data.frame()
# 
# # Loop through all time series and perform bootstrap BDS test
# for (uuid in names(results.ARWN)) {
#   for (emotion in names(results.ARWN[[uuid]])) {
#     
#     # Progress bar update
#     iteration <- iteration + 1
#     setTxtProgressBar(pb, iteration)
#     
#     res <- results.ARWN[[uuid]][[emotion]]
#     
#     # Extract measurement-error-free residuals
#     residuals <- res$residuals
#     
#     # Perform bootstrap BDS test
#     tryCatch({
#       bds_boot_result <- bds_test_bootstrap(residuals, m = 2, eps = 1, n_bootstrap = 499)
#       
#       # Get Heywood status
#       heywood_row <- heywood_summary[heywood_summary$UUID == uuid & heywood_summary$Emotion == emotion, ]
#       has_heywood <- heywood_row$heywood_ivar | heywood_row$heywood_evar
#       
#       # Add to results
#       bds_results <- rbind(bds_results, data.frame(
#         UUID = uuid,
#         Emotion = emotion,
#         n_obs = length(residuals),
#         bds_statistic = bds_boot_result$bds_statistic,
#         bds_pvalue_parametric = bds_boot_result$bds_pvalue_parametric,
#         bds_pvalue_empirical = bds_boot_result$bds_pvalue_empirical,
#         heywood_ivar = heywood_row$heywood_ivar,
#         heywood_evar = heywood_row$heywood_evar,
#         has_heywood = has_heywood,
#         residuals_mean = mean(residuals),
#         residuals_sd = sd(residuals),
#         stringsAsFactors = FALSE
#       ))
#     }, error = function(e) {
#       cat("Error in bootstrap BDS test for", uuid, "-", emotion, ":", e$message, "\n")
#     })
#   }
# }
# 
# close(pb)
# 
# # Check results
# head(bds_results)
# 
# 
# # --- Compare Heywood cases and valid cases ---
# heywood_cases <- bds_results[bds_results$has_heywood == TRUE,] # 219
# valid_cases <- bds_results[bds_results$has_heywood == FALSE,] # 141
# 
# # Test 1: Difference of BDS statistic
# t.test(valid_cases$bds_statistic, heywood_cases$bds_statistic)
# 
# # Welch Two Sample t-test
# # 
# # data:  valid_cases$bds_statistic and heywood_cases$bds_statistic
# # t = -1.5242, df = 356.39, p-value = 0.1284
# # alternative hypothesis: true difference in means is not equal to 0
# # 95 percent confidence interval:
# #   -1.2065712  0.1529391
# # sample estimates:
# #   mean of x mean of y 
# # 3.547910  4.074726 
# 
# ### -> No difference in BDS statistic.
# 
# # Test 2: Difference in probability of significant Bootstrap BDS tests
# t <- matrix(c(
#   sum(valid_cases$bds_pvalue_empirical < 0.05),
#   sum(valid_cases$bds_pvalue_empirical >= 0.05),
#   sum(heywood_cases$bds_pvalue_empirical < 0.05),
#   sum(heywood_cases$bds_pvalue_empirical >= 0.05)
# ), nrow = 2, byrow = TRUE)
# 
# rownames(t) <- c("Valid Cases", "Heywood Cases")
# colnames(t) <- c("Reject H0", "Fail to Reject H0")
# chisq.test(t)
# # Pearson's Chi-squared test with Yates' continuity correction
# # 
# # data:  t
# # X-squared = 1.4954e-30, df = 1, p-value = 1
# 
# ### -> No difference in probability of significant bootstrap test.