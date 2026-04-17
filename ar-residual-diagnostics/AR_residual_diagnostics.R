
#####################################
### AR RESIDUAL DIAGNOSTICS #########
#####################################

library(dplyr)
library(FKF)
library(tseries) 
library(tidyverse)
library(future)
library(furrr)

# set wd 
# setwd("...")
setwd("C:/Users/Sebastian Küppers/Desktop/Formal Theory of Co-Occuring Emotions (DFG project)/_PhD/_PhD_Study_1/Research Exchange/PROJECT/Preregistration/github/ar-residual-diagnostics")

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
# 
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

log_message <- sprintf(
  "[%s] Dauer: %.2f Minuten\n",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  as.numeric(duration)
)
cat(log_message, file = "log.txt", append = TRUE)

# Save
saveRDS(results_parallel, file = "ARWN_results.rds")

plan(sequential)


# Read data
results <- readRDS("ARWN_results.rds")

## --------------------------------
# CONVERGENCE AND HEYWOOD CASES ---
## --------------------------------

nonvalid_rows_bayesian <- c()

for (i in 1:length(results$result)) {
  bayesian <- results$result[[i]]$bayesian
  if (is.null(bayesian)) nonvalid_rows_bayesian <- c(nonvalid_rows_bayesian, i)
}

length(nonvalid_rows_bayesian) # 20

# ---> 20 rows with problems in Bayesian estimation. Can possibly be reduced by
# ---> tweaking init parameters.


## --------------------------------
# HEYWOOD CASES -------------------

# Frequentist Heywood (any)
freq_heywood_idx <- which(map_lgl(results$result, 
                                  ~.$frequentist$heywood_ivar | .$frequentist$heywood_evar))
length(freq_heywood_idx) # 219
length(freq_heywood_idx) / nrow(results) # 0.608

# Bayesian Heywood (any, only successful)
bayes_heywood_idx <- which(map_lgl(results$result, 
                                   ~!is.null(.$bayesian) && (.$bayesian$heywood_ivar | .$bayesian$heywood_evar)))

length(bayes_heywood_idx) # 0
