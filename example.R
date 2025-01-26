# Load necessary packages
library(data.table)
library(fixest)

# Load your custom appml function
source("appml_r2.R")

set.seed(42)
n <- 1000
data <- data.table(
  id = rep(1:(n/10), each = 10),        # Fixed effect group
  time = rep(1:10, n/10),               # Time fixed effect
  x = rnorm(n, 5, 2),                   # Independent variable
  z = rbinom(n, 1, 0.5)                 # Binary treatment variable
)

# Add a strong relationship between y, x, and z
beta_x <- 0.3  # Coefficient for x
beta_z <- 1.5  # Coefficient for z
noise_sd <- 2  # Standard deviation of noise

# Generate y with a meaningful relationship
data[, y := rpois(n, lambda = exp(beta_x * x + beta_z * z + rnorm(n, 0, noise_sd)))]

# Define model formula
y_var <- "y"
x_vars <- c("z", "x")
formula <- as.formula(paste(y_var, "~", paste(x_vars, collapse = " + ")))

# Define fixed effects as character vector (correct type for fixest)
fixef_vars <- c("id", "time")

# Initialize expectile range: starts at 0.5, moves outward
expectiles_down <- seq(0.5, 0.02, by = -0.01)
expectiles_up <- seq(0.5, 0.98, by = 0.01)
expectiles <- c(expectiles_down, setdiff(expectiles_up, 0.5))

# Run appml function for each expectile
results <- list()
start_values <- NULL # Initial residuals for iterative estimation

for (tau in expectiles) {
  cat("\nRunning expectile regression for τ =", tau, "...\n")
  
  result <- tryCatch({
    appml(
      formula = formula,
      data = data,
      expectile = tau,
      start = start_values,   # Use residuals from the previous step
      absorb = fixef_vars,    # Pass fixed effects as character vector
      cluster = ~id           # Cluster by id for robust standard errors
    )
  }, error = function(e) {
    cat("Error for τ =", tau, ":", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    results[[as.character(tau)]] <- result
    # Update start_values for the next iteration
    start_values <- result$residuals
  }
}

# Extract and summarize results
summary_df <- data.table(
  Expectile = numeric(),
  Coefficient_z = numeric(),
  Coefficient_x = numeric(),
  Iterations = numeric(),
  Objective = numeric()
)

for (tau in names(results)) {
  res <- results[[tau]]
  if (!is.null(res)) {
    summary_df <- rbind(summary_df, data.table(
      Expectile = as.numeric(tau),
      Coefficient_z = res$coefficients["z"],
      Coefficient_x = res$coefficients["x"],
      Iterations = res$iterations,
      Objective = res$objective_function
    ))
  }
}

# Print the summary table
print(summary_df)

# Plot coefficients across expectiles
library(ggplot2)
ggplot(summary_df, aes(x = Expectile)) +
  geom_line(aes(y = Coefficient_z), color = "blue", size = 1) +
  geom_point(aes(y = Coefficient_z), color = "blue") +
  geom_line(aes(y = Coefficient_x), color = "red", size = 1) +
  geom_point(aes(y = Coefficient_x), color = "red") +
  labs(title = "Expectile Regression Coefficients",
       x = "Expectile (τ)",
       y = "Coefficient") +
  theme_minimal()
