# The software is provided "as is," without warranty of any kind, express or implied, 
# including but not limited to the warranties of merchantability, fitness for a particular 
# purpose, and noninfringement. 
# 
# In no event shall the authors be liable for any claim, damages, or other liability, 
# whether in an action of contract, tort, or otherwise, arising from, out of, or in connection 
# with the software or its use, or other dealings with the software.

library(fixest)
library(data.table)

appml <- function(formula, data, expectile = 0.5, tolerance = 1e-12, iterate = 50, 
                  start = NULL, nolog = FALSE, residual_name = NULL, 
                  absorb = NULL, # specify as ~fe_var if needed
                  cluster = NULL, # specify as ~cluster_var if needed
                  maxiter = 200, strict = FALSE) {
  
  # Convert data to data.table if not already
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }
  
  # Extract the dependent variable name
  y_var <- as.character(formula[[2]])
  
  # Initial values
  if (is.null(start)) {
    # Initial fit without weights
    # Adjust the formula with absorb and cluster if needed
    # For example, if absorb = ~fe and cluster = ~pairid:
    initial_fit <- fixest::feglm(fml = formula,
                         data = data,
                         family = poisson(),
                         fixef = if (!is.null(absorb)) absorb else NULL,
                         cluster = cluster,
                         maxiter = maxiter)
    fitted_values <- fitted(initial_fit, na.rm = FALSE)
    data[, residuals := get(y_var) - fitted_values]
  } else {
    data[, residuals := start]
  }
  
  # Initial weights
  data[, weights := abs(expectile - (residuals < 0))]
  
  # For checking convergence
  cv <- 100
  count <- 0
  bold <- NULL
  prev_vcov <- NULL
  
  # Iteration
  while (cv > tolerance && count < iterate) {
    count <- count + 1
    
    # Fit with current weights
    fit <- tryCatch({
      fixest::feglm(fml = formula,
            data = data,
            family = poisson(),
            weights = data$weights,
            fixef = fixef_vars,
            cluster = cluster,
            maxiter = maxiter)
    }, error = function(e) {
      cat("Error during feglm fitting at iteration", count, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(fit)) {
      cat("Null fit at iteration", count, "\n")
      break
    }
    
    # Compute new residuals
    fitted_values <- fitted(fit, na.rm = FALSE)
    new_residuals <- data[[y_var]] - fitted_values
    
    # Extract coefficients
    b <- coef(fit)
    # The ado code sets the last coefficient to zero each time:
    # mat `b'[1,`k']=0
    # If you need to enforce this condition (depends on your normalization),
    # you can do something like:
    # b[length(b)] <- 0
    # But this may not always make sense in R. Consider leaving as-is.
    
    # Compute the objective function 'cv'
    # In the ado code:
    # cv = (b - bold)*invsym(e(V))*(b - bold)'
    # On the first iteration, if we don't have bold, set bold to zero vector.
    if (count == 1) {
      # If start was provided, bold is zero vector initially
      # Otherwise, also start bold as zero vector
      bold <- rep(0, length(b))
    }
    
    # Get vcov matrix from the model
    current_vcov <- vcov(fit)
    
    # Inverse of vcov
    invV <- tryCatch({
      solve(current_vcov)
    }, error = function(e) {
      cat("Error inverting vcov at iteration", count, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(invV)) break
    
    diff_b <- (b - bold)
    cv <- as.numeric(t(diff_b) %*% invV %*% diff_b)
    
    # Update residuals and weights
    data[, residuals := new_residuals]
    data[, weights := abs(expectile - (residuals < 0))]
    
    # Print iteration info if not suppressed
    if (!nolog) {
      cat("Iteration", count, ": objective function =", cv, "\n")
    }
    
    # Update bold and vcov for next iteration
    bold <- b
    prev_vcov <- current_vcov
  }
  
  # Check convergence
  converged <- cv <= tolerance
  
  # Final output
  negative_share <- mean(data$residuals < 0)
  
  # If user wants residuals saved under a particular name
  if (!is.null(residual_name)) {
    data[[residual_name]] <- data$residuals
  }
  
  result <- list(
    fit = fit,
    converged = converged,
    iterations = count,
    tolerance = tolerance,
    objective_function = cv,
    negative_residuals = negative_share,
    expectile = expectile,
    coefficients = coef(fit),
    vcov = vcov(fit),
    residuals = data$residuals,
    num_obs = nrow(data)
  )
  
  # Print final summary
  cat("\nNumber of obs =", result$num_obs, "\n")
  cat("Iterations =", count, "\n")
  cat("Tolerance =", tolerance, "\n")
  cat("Objective function =", cv, "\n")
  cat("% negative residuals =", round(100 * negative_share, 3), "%\n")
  cat("Expectile =", expectile, "expectile regression\n")
  
  return(result)
}
