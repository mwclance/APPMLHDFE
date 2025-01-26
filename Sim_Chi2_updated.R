# --------------------------------------------------------------
# Revised Simulation Script Converted from Stata to R
# --------------------------------------------------------------
# Date: 2024-12-21
# Description: This script replicates a Stata simulation in R,
#              performing expectile regressions at specified levels,
#              utilizing a custom `appml` function for model estimation,
#              and enabling parallel processing for enhanced performance.
# --------------------------------------------------------------

# --------------------------------------------------------------
# Generating Data Outside the Loop & Labeling Expectiles
# --------------------------------------------------------------

# -----------------------------
# 1. Setup and Initialization
# -----------------------------

rm(list = ls())

required_packages <- c(
  "data.table",
  "fixest",
  "parallel",
  "doParallel",
  "foreach",
  "log4r",
  "stringr"
)

install_if_missing <- function(packages) {
  installed <- rownames(installed.packages())
  for (pkg in packages) {
    if (!(pkg %in% installed)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

install_if_missing(required_packages)

library(data.table)
library(fixest)
library(parallel)
library(doParallel)
library(foreach)
library(log4r)
library(stringr)

# Set seed for reproducibility
set.seed(10122024)

# Logging
logger <- create.logger()
logfile(logger) <- "Sims_Chi2.log"
level(logger) <- "INFO"
info(logger, "Simulation started")

# Chi-square threshold
invchi2_threshold <- qchisq(0.15, df = 4)  # ~15% prob with df=4

# -----------------------------
# 2. Generate Data Structure Outside the Loop
# -----------------------------

generate_data_structure <- function(n, t) {
  obs <- n^2
  IDi <- ceiling(1:obs / n)
  IDe <- ((1:obs - 1) %% n) + 1
  pair <- (IDi - 1) * n + IDe
  
  total_obs <- obs * t
  IDi_expanded <- rep(IDi, each = t)
  IDe_expanded <- rep(IDe, each = t)
  pair_expanded <- rep(pair, each = t)
  T_time <- rep(1:t, times = obs)
  
  data.table(IDi = IDi_expanded,
             IDe = IDe_expanded,
             pair = pair_expanded,
             T = T_time)
}

# -----------------------------
# 3. Define `appml` for Expectile Poisson
# -----------------------------

appml <- function(formula, data, expectile = 0.5, tolerance = 1e-12, iterate = 50, 
                  start = NULL, nolog = FALSE, residual_name = NULL, 
                  absorb = NULL, cluster = NULL, maxiter = 200, strict = FALSE) {
  
  # Convert data to data.table
  if (!is.data.table(data)) data <- as.data.table(data)
  
  y_var <- as.character(formula[[2]])
  
  # Initial values
  if (is.null(start)) {
    initial_fit <- fixest::feglm(
      fml = formula,
      data = data,
      family = poisson(),
      fixef = if (!is.null(absorb)) terms(as.formula(absorb)) else NULL,
      cluster = if (!is.null(cluster)) as.character(as.formula(cluster))[[2]] else NULL,
      maxiter = maxiter
    )
    fitted_values <- fitted(initial_fit, na.rm = FALSE)
    data[, residuals := get(y_var) - fitted_values]
  } else {
    data[, residuals := start]
  }
  
  # Weights based on expectile
  data[, weights := abs(expectile - (residuals < 0))]
  
  cv <- 100
  count <- 0
  bold <- NULL
  
  while (cv > tolerance && count < iterate) {
    count <- count + 1
    
    fit <- tryCatch({
      fixest::feglm(
        fml = formula,
        data = data,
        family = poisson(),
        weights = data$weights,
        fixef = if (!is.null(absorb)) terms(as.formula(absorb)) else NULL,
        cluster = if (!is.null(cluster)) as.character(as.formula(cluster))[[2]] else NULL,
        maxiter = maxiter
      )
    }, error = function(e) {
      cat("Error at iteration", count, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(fit)) {
      cat("Null fit at iteration", count, "\n")
      break
    }
    
    fitted_values <- fitted(fit, na.rm = FALSE)
    new_residuals <- data[[y_var]] - fitted_values
    b <- coef(fit)
    
    if (count == 1) {
      bold <- rep(0, length(b))
    }
    
    current_vcov <- vcov(fit)
    invV <- tryCatch({
      solve(current_vcov)
    }, error = function(e) {
      cat("Error inverting vcov:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(invV)) break
    
    diff_b <- (b - bold)
    cv <- as.numeric(t(diff_b) %*% invV %*% diff_b)
    
    data[, residuals := new_residuals]
    data[, weights := abs(expectile - (residuals < 0))]
    
    if (!nolog) {
      cat("Iteration", count, ": objective function =", cv, "\n")
    }
    
    bold <- b
  }
  
  converged <- cv <= tolerance
  neg_share <- mean(data$residuals < 0)
  
  if (!is.null(residual_name)) data[[residual_name]] <- data$residuals
  
  list(
    fit = fit,
    converged = converged,
    iterations = count,
    tolerance = tolerance,
    objective_function = cv,
    negative_residuals = neg_share,
    expectile = expectile,
    coefficients = coef(fit),
    vcov = vcov(fit),
    residuals = data$residuals,
    num_obs = nrow(data)
  )
}

# -----------------------------
# 4. Precompute Data Structures for Each (n, t)
# -----------------------------

# Parameter grid
betas <- c(0.2, 0)
ns <- c(45, 90)
ts <- c(15, 30, 45)
param_grid <- CJ(beta = betas, n = ns, t = ts)

# Precompute data structures
data_structures <- list()
for (i in 1:nrow(param_grid)) {
  n_i <- param_grid$n[i]
  t_i <- param_grid$t[i]
  key <- paste0(n_i, "_", t_i)
  data_structures[[key]] <- generate_data_structure(n_i, t_i)
}

# -----------------------------
# 5. Replication Function Using Precomputed Data
# -----------------------------

run_replication <- function(beta, data_structure, rep_id) {
  options(warn = -1)
  
  # Copy the precomputed structure
  data <- copy(data_structure)
  
  # Generate random components each replication
  data[, eta1 := rnorm(.N)]
  data[, me1 := mean(eta1), by = IDi]
  data[, fe1 := me1]
  
  data[, eta2 := rnorm(.N)]
  data[, me2 := mean(eta2), by = IDe]
  data[, fe2 := me2]
  
  data[, eta3 := rnorm(.N)]
  data[, me3 := mean(eta3), by = pair]
  data[, fe3 := me3]
  
  # Generate dummy variable d
  data[, d := as.integer(rchisq(.N, df = 4) < invchi2_threshold)]
  if (length(unique(data$d)) < 2) {
    return(c(b_0.5=NA, t_0.5=NA, b_0.1=NA, t_0.1=NA, b_0.3=NA, t_0.3=NA, 
             b_0.7=NA, t_0.7=NA, b_0.9=NA, t_0.9=NA))
  }
  
  # Generate y
  data[, mu := exp(beta * d + 0.4 * (fe1^2 + fe2^2 + fe3^2))]
  data[, y := mu * (rchisq(.N, df = 2)) / 1e3]
  
  # Custom name vector to reflect actual expectile levels
  labels_map <- c("b_0.5","t_0.5","b_0.1","t_0.1","b_0.3","t_0.3","b_0.7","t_0.7","b_0.9","t_0.9")
  result <- numeric(length(labels_map))
  names(result) <- labels_map
  
  # Expectile order
  expectile_levels <- c(0.5, 0.1, 0.3, 0.7, 0.9)
  
  # Store residuals from 0.5 if needed
  residuals_0.5 <- NULL
  
  for (e in expectile_levels) {
    formula_fixed <- as.formula("y ~ d | IDe:T + IDi:T + pair")
    
    # Decide on 'start' values
    if (e == 0.5) {
      start_vals <- NULL
    } else if (!is.null(residuals_0.5)) {
      start_vals <- residuals_0.5
    } else {
      start_vals <- NULL
    }
    
    fit <- appml(
      formula = formula_fixed,
      data = data,
      expectile = e,
      tolerance = 1e-4,
      iterate = 50,
      start = start_vals,
      nolog = TRUE,
      residual_name = NULL,
      absorb = NULL,
      cluster = NULL,
      maxiter = 200,
      strict = FALSE
    )
    
    if (!is.null(fit$converged) && fit$converged && !any(is.na(fit$coefficients))) {
      coef_d <- fit$coefficients["d"]
      se_d <- sqrt(fit$vcov["d","d"])
      t_stat <- (coef_d - beta) / se_d
      
      # Store in named vector
      result[paste0("b_", e)] <- coef_d
      result[paste0("t_", e)] <- t_stat
      
      # If e == 0.5, keep residuals
      if (e == 0.5) {
        residuals_0.5 <- fit$residuals
      }
    } else {
      # No convergence or invalid
      result[paste0("b_", e)] <- NA
      result[paste0("t_", e)] <- NA
    }
  }
  
  return(result)
}

# --------------------------------------------------------------
# 6. Main Simulation Loop
# --------------------------------------------------------------

replications <- 1000  # For testing; set higher for actual runs

num_cores <- detectCores() - 1
#num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)

all_results <- list()

for (i in 1:nrow(param_grid)) {
  beta <- param_grid$beta[i]
  n_val <- param_grid$n[i]
  t_val <- param_grid$t[i]
  
  # Key for the precomputed data structure
  key <- paste0(n_val, "_", t_val)
  data_structure <- data_structures[[key]]
  
  param_start_time <- Sys.time()
  
  info(logger, paste("Starting simulations for beta =", beta, 
                     ", n =", n_val, ", t =", t_val))
  
  replication_results <- foreach(r = 1:replications, .combine = rbind,
                                 .packages = c("data.table", "fixest")) %dopar% {
                                   run_replication(beta, data_structure, r)
                                 }
  
  replication_dt <- as.data.table(replication_results)
  
  # Example coverage logic
  # coverage for each expectile: abs(t_stat) < 1.96
  # Columns match the order in our 'labels_map'
  # [b_0.5, t_0.5, b_0.1, t_0.1, b_0.3, t_0.3, b_0.7, t_0.7, b_0.9, t_0.9]
  replication_dt[, c_0.5 := abs(replication_dt$t_0.5) < 1.96]
  replication_dt[, c_0.1 := abs(replication_dt$t_0.1) < 1.96]
  replication_dt[, c_0.3 := abs(replication_dt$t_0.3) < 1.96]
  replication_dt[, c_0.7 := abs(replication_dt$t_0.7) < 1.96]
  replication_dt[, c_0.9 := abs(replication_dt$t_0.9) < 1.96]
  
  summary_stats <- list(
    beta = beta,
    n = n_val,
    t = t_val,
    # Round medians, coverage, etc. to 5 decimal places
    median_b_0.5 = round(median(replication_dt$b_0.5, na.rm = TRUE), 5),
    median_b_0.1 = round(median(replication_dt$b_0.1, na.rm = TRUE), 5),
    median_b_0.3 = round(median(replication_dt$b_0.3, na.rm = TRUE), 5),
    median_b_0.7 = round(median(replication_dt$b_0.7, na.rm = TRUE), 5),
    median_b_0.9 = round(median(replication_dt$b_0.9, na.rm = TRUE), 5),
    coverage_0.5 = round(mean(replication_dt$c_0.5, na.rm = TRUE), 5),
    coverage_0.1 = round(mean(replication_dt$c_0.1, na.rm = TRUE), 5),
    coverage_0.3 = round(mean(replication_dt$c_0.3, na.rm = TRUE), 5),
    coverage_0.7 = round(mean(replication_dt$c_0.7, na.rm = TRUE), 5),
    coverage_0.9 = round(mean(replication_dt$c_0.9, na.rm = TRUE), 5),
    duration_secs = NA,
    failure_rate = NA
  )
  
  param_end_time <- Sys.time()
  duration_secs <- as.numeric(difftime(param_end_time, param_start_time, units = "secs"))
  summary_stats$duration_secs <- duration_secs
  
  failed_replications <- replication_dt[
    is.na(b_0.5) | is.na(b_0.1) | is.na(b_0.3) | is.na(b_0.7) | is.na(b_0.9)
  ]
  failure_rate <- nrow(failed_replications) / replications
  summary_stats$failure_rate <- failure_rate
  
  info(logger, paste("Completed simulations for beta =", beta, 
                     ", n =", n_val, ", t =", t_val))
  print(summary_stats)
  
  all_results[[i]] <- list(replications = replication_dt, summary = summary_stats)
}

stopCluster(cl)
info(logger, "Parallel backend stopped.")

summary_table <- rbindlist(lapply(all_results, function(x) x$summary))
fwrite(summary_table, "Simulation_Summary.csv")
info(logger, "Summary statistics saved to Simulation_Summary.csv")

info(logger, "Simulation completed successfully.")
info(logger, "Log file saved as Sims_Chi2.log")

# Load the CSV file into a data.table
summary_table <- fread("Simulation_Summary.csv")

# Reorder the columns
reordered_columns <- c(
  "beta", "n", "t",
  "median_b_0.1", "coverage_0.1", 
  "median_b_0.3", "coverage_0.3", 
  "median_b_0.5", "coverage_0.5", 
  "median_b_0.7", "coverage_0.7", 
  "median_b_0.9", "coverage_0.9", 
  "duration_secs", "failure_rate"
)
# Reorder columns in the data.table
summary_table <- summary_table[, ..reordered_columns]

# Sum the duration_secs column
total_duration_secs <- sum(summary_table$duration_secs, na.rm = TRUE)

# Convert the total duration to hours and minutes
total_duration_hours <- floor(total_duration_secs / 3600)
total_duration_minutes <- round((total_duration_secs %% 3600) / 60)

# Print the total duration
cat("Total duration:", total_duration_hours, "hours and", total_duration_minutes, "minutes\n")

# Save the reordered table to a new CSV file
fwrite(summary_table, "Reordered_Simulation_Summary.csv")

# Print a message indicating completion
cat("Reordered simulation summary saved as 'Reordered_Simulation_Summary.csv'\n")
