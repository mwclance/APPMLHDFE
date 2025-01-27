rm(list = ls())
gc()


# Load necessary packages
library(fixest)
library(data.table)
library(readstata13)
library(ggplot2)

# Set working directory based on your environment variable
computer_name <- Sys.getenv("COMPUTERNAME")

if (computer_name == "MWCLANCE") {
  setwd("D:/Work/International Economics/Expectiles/Empirical/Complete/R")
} else if (computer_name == "LEV1") {
  setwd("D:/Work/International Economics/Expectiles/Empirical/Complete/R")
} else {
  stop("Unknown computer name: ", computer_name)
}

# Source the appml function
source("appml_r2.R")

# Load the data
data <- read.dta13("new_data.dta")

# Convert to data.table for consistency
setDT(data)

# Subset the data to include only the required variables
required_vars <- c("iso_o", "iso_d", "year", "pairid", "trade_x", "gdp_o", "EIA", "EIAn", 
                   "expyear2", "impyear2")
data <- data[, ..required_vars]

# Data preparation
data[, INTER := as.integer(iso_o != iso_d)]
setorder(data, year, pairid)
data[, trade := trade_x]
data[is.na(trade), trade := 0]
data[, agg_exports := sum(trade), by = .(iso_o, year)]
data[iso_o == iso_d, trade := gdp_o - agg_exports]
data[is.na(trade) | trade < 0, trade := 0]
data[, trade_all := trade / 1e9]
data <- data[year >= 1962]
data[iso_o == iso_d, rta := 0]

# Create necessary variables
data[, `:=`(
  CUCMECU = ifelse(EIAn >= 4 & !is.na(EIAn), 1, 0),
  EIAp = ifelse(EIAn > 0 & !is.na(EIAn), 1, 0),
  EIA123 = ifelse(EIAn > 0 & EIAn < 4 & !is.na(EIAn), 1, 0),
  EIA12 = ifelse(EIAn > 0 & EIAn < 3 & !is.na(EIAn), 1, 0),
  EIA45 = ifelse(EIAn > 3 & EIAn < 6 & !is.na(EIAn), 1, 0),
  EIA1 = ifelse(EIAn == 1 & !is.na(EIAn), 1, 0),
  EIA2 = ifelse(EIAn == 2 & !is.na(EIAn), 1, 0),
  EIA3 = ifelse(EIAn == 3 & !is.na(EIAn), 1, 0),
  EIA6 = ifelse(EIAn == 6 & !is.na(EIAn), 1, 0)
)]

# Remove rows with missing EIA values
data <- data[!is.na(EIAp)]

# Create interaction term for INTER and year
data[, INTER_YEAR := interaction(INTER, year, sep = "_")]

# Set the formula without fixed effects in the formula itself
# Only the dependent variable and regressors go here:
y <- "trade_all"
regs <- "EIAp"
full_formula <- as.formula(paste(y, "~", regs))

# Define fixed effects as a character vector
fixef_vars <- c("INTER_YEAR", "expyear2", "impyear2", "pairid")

# Initial model to determine the sample
# Use feglm directly to get initial fitted values:
initial_fit <- feglm(fml = full_formula,
                     data = data,
                     family = poisson(),
                     fixef = fixef_vars)

# Predict fitted values and subset data
data$fitted <- fitted(initial_fit, na.rm = FALSE)
data_est <- data[!is.na(fitted)]

# Calculate initial residuals
data_est[, initial_residuals := trade_all - fitted]

# Define the expectiles to be tested
expectiles_down <- seq(0.50, 0.02, by = -0.01)
expectiles_up <- seq(0.50, 0.99, by = 0.01)

# Function to run appml for a single expectile
run_expectile <- function(expectile, formula, data, start_values = NULL) {
  tryCatch({
    result <- appml(
      formula = formula,
      data = data,
      expectile = expectile,
      start = start_values,
      absorb = fixef_vars,  # Absorb the fixed effects via this argument
      cluster = ~pairid,    # Cluster by pairid if needed
      maxiter = 200
    )
    return(list(expectile = expectile, result = result))
  }, error = function(e) {
    return(list(expectile = expectile, error = e$message))
  })
}

# Process expectiles in sequence
process_expectiles_sequential <- function(expectiles, data, formula) {
  results <- list()
  start_values <- NULL
  for (expectile in expectiles) {
    res <- run_expectile(expectile, formula, data, start_values)
    if (!is.null(res$result)) {
      start_values <- res$result$residuals
    }
    results[[as.character(expectile)]] <- res
  }
  return(results)
}

# Process expectiles
start.time <- Sys.time()
results_down <- process_expectiles_sequential(expectiles_down, data_est, full_formula)
results_up <- process_expectiles_sequential(expectiles_up, data_est, full_formula)
end.time <- Sys.time()

# Combine and save results
result_df <- data.table(expectile = numeric(), coefficient = numeric(), se = numeric(), ci_low = numeric(), ci_high = numeric())

all_results <- c(results_down, results_up)
for (res in all_results) {
  if (!is.null(res$result)) {
    result <- res$result
    coef_val <- result$coefficients[regs]
    cat("Expectile:", res$expectile, "Coefficient:", coef_val, "\n") # Debug
    se_val <- sqrt(result$vcov[1, 1])
    cat("SE for", regs, ":", se_val, "\n")
    result_df <- rbind(result_df, data.table(
      expectile = res$expectile,
      coefficient = coef_val,
      se = se_val,
      ci_low = coef_val - 1.96 * se_val,
      ci_high = coef_val + 1.96 * se_val
    ))
  }
}

fwrite(result_df, "results__mono.csv", sep = ";")

# Plot results
ggplot(result_df, aes(x = expectile, y = coefficient)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), fill = "grey80") +
  geom_line(color = "blue") +
  labs(title = "Expectile Regression Results", x = "Expectile", y = "Coefficient") +
  theme_minimal()

ggsave("expectile_regression_results.png")

time.taken <- end.time - start.time
print(paste("Time taken: ", round(time.taken, 2)))

# Extract specific expectiles and compare
compare_expectiles <- function(df, target, tol = .Machine$double.eps^0.5) {
  df[abs(expectile - target) < tol]
}

expectile_10th <- compare_expectiles(result_df, 0.1)
expectile_50th <- compare_expectiles(result_df, 0.5)
expectile_90th <- compare_expectiles(result_df, 0.9)

selected_expectiles <- rbind(expectile_10th, expectile_50th, expectile_90th)
selected_expectiles <- unique(selected_expectiles, by = "expectile")

print(selected_expectiles)

print_table <- selected_expectiles[, .(
  Expectile = c("10th", "50th", "90th"),
  Coefficient = round(coefficient, 6),
  SE = round(se, 9)
)]

cat("Results for Selected Expectiles\n")
cat("-------------------------------------------------\n")
cat(sprintf("%-10s %-15s %-15s\n", "Expectile", "Coefficient", "SE"))
cat("-------------------------------------------------\n")

for (i in 1:nrow(print_table)) {
  cat(sprintf("%-10s %-15.6f %-15.9f\n", 
              print_table$Expectile[i], 
              print_table$Coefficient[i], 
              print_table$SE[i]))
}

cat("-------------------------------------------------\n")
cat("Note: Clustered standard errors by country-pair.\n")
