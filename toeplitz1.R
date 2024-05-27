library(glmnet)
library(ggplot2)
library(MASS)
library(Matrix)

# Function to generate Toeplitz matrix (with PD fix using nearPD)
generate_toeplitz <- function(t1, t2, t3) {
  T <- matrix(c(t1, t2, t3, t2, t1, t2, t3, t2, t1), nrow = 3, byrow = TRUE)
  # Ensure positive definiteness using nearPD
  if (Matrix::is.positive.definite(T)) {
    return(T)
  } else {
    return(Matrix::nearPD(T)$mat)  # Use nearPD if not positive definite
  }
}

# Function to compute SURE estimate of GBEC score (with singularity check)
compute_sure_gbec <- function(X, Y, lambda_seq, Sigma) {
  n <- nrow(X)
  model <- glmnet(X, Y, alpha = 0, lambda = lambda_seq, intercept = FALSE)
  df <- model$df
  sigma2_hat <- colMeans((Y - predict(model, X))^2) / (1 - df/n)
  sure <- colSums((Y - predict(model, X))^2) + 2 * sigma2_hat * df - n * sigma2_hat
  best_lambda <- model$lambda[which.min(sure)]
  
  J <- matrix(1, ncol(X), 1)
  
  # Singular Value Decomposition of Sigma to assess rank
  svd_result <- svd(Sigma)
  
  # If matrix is close to singular return NA
  if (min(svd_result$d) / max(svd_result$d) < 1e-7) {
    return(NA) # No warning here
  } else {
    gamma_sure <- solve(Sigma + best_lambda * diag(ncol(X))) %*% J / as.numeric(t(J) %*% solve(Sigma + best_lambda * diag(ncol(X))) %*% J)
    gbec_score <- as.numeric(t(gamma_sure) %*% t(X))
    return(gbec_score)
  }
}

# Function to compute ridge estimate of GBEC score (with singularity check)
compute_ridge_gbec <- function(X, Y, lambda, Sigma) {
  # Singular Value Decomposition of Sigma to assess rank
  svd_result <- svd(Sigma)
  # If matrix is close to singular return NA
  if (min(svd_result$d) / max(svd_result$d) < 1e-7) {
    return(NA) # No warning here
  } 
  
  model <- glmnet(X, Y, alpha = 0, lambda = lambda, intercept = FALSE)
  gamma_ridge <- coef(model)[-1]
  gbec_score <- as.numeric(t(gamma_ridge) %*% t(X))
  return(gbec_score)
}

# Function to compute BEC using pseudo-inverse (with error handling and improved inversion)
compute_bec <- function(X, Y, tol = 1e-5, Sigma) {
  
  # Ensure Sigma is a matrix (not a dataframe)
  Sigma <- as.matrix(Sigma)
  
  # Ensure J is a column vector and has the correct dimensions
  J <- matrix(1, ncol(Sigma), 1)
  
  # Calculate pseudo-inverse using SVD with tolerance
  Sigma_pseudo_inv <- tryCatch(
    {
      ginv(Sigma, tol = tol)
    },
    error = function(e) {
      warning(paste("Error inverting Sigma in BEC:", e$message))
      return(NA)
    }
  )
  
  # Check if Sigma_pseudo_inv was calculated successfully
  if (any(is.na(Sigma_pseudo_inv))) {
    return(rep(NA, nrow(Sigma_pseudo_inv))) # Return NAs if there's an error
  }
  
  # Explicitly ensure dimensions match
  if (round(ncol(Sigma_pseudo_inv)) != round(nrow(J))) {
    stop("Dimension mismatch between Sigma_pseudo_inv and J (even after rounding).")
  }
  
  c_opt <- tryCatch(
    {
      as.numeric(1 / (t(J) %*% Sigma_pseudo_inv %*% J))  # Force numeric result
    },
    error = function(e) {
      warning(paste("Error calculating c_opt in BEC:", e$message))
      NA  # Return NA if there's an error
    }
  )
  
  # If c_opt is not NA, calculate gamma
  if (!is.na(c_opt)) {
    gamma <- c_opt * Sigma_pseudo_inv %*% J
  } else {
    gamma <- rep(NA, nrow(Sigma_pseudo_inv))  # Return NAs for gamma if c_opt is NA
  }
  
  return(drop(gamma))
}


# Simulation settings
set.seed(123)
n <- 100
num_simulations <- 100 
t1_values <- seq(0.8, 1.2, by = 0.05) 
t2 <- 1.0
t3 <- 0.9
lambda_seq <- seq(0.01, 1, length.out = 100)


# Function to calculate the MSE (Mean Squared Error)
calculate_mse <- function(estimator_func, ..., true_gbec) {
  mse_values <- replicate(num_simulations, {
    # Obtain the estimated coefficients (gamma) or Gbec score(y)
    estimated_values <- tryCatch(
      estimator_func(...),
      error = function(e) {
        rep(NA, num_simulations) # Return NAs if there's an error
      }
    )
    
    # Filter out NA values before calculating MSE
    valid_estimates <- estimated_values[!is.na(estimated_values)]
    if(length(valid_estimates)>0) {
      mean((valid_estimates - true_gbec)^2, na.rm = TRUE) # Calculate MSE only on valid estimates
    } else {
      NA
    }
  })
  # Check if mse_values has all NAs
  if(all(is.na(mse_values))){
    return(NA)
  } else {
    return(mean(mse_values, na.rm = TRUE)) # Return the average MSE (ignoring NAs)
  }
}

# Create an empty dataframe to store the results
results_mse <- data.frame(t1 = numeric(), Method = character(), MSE = numeric())
# Perform the simulation and calculate MSEs for SURE
for(t1 in t1_values) {
  Sigma <- generate_toeplitz(t1, t2, t3)
  X <- mvrnorm(n, rep(0, 3), Sigma)
  Y <- rowSums(X) + rnorm(n, sd = 0.5)
  
  # Obtain the true GBEC score (adjusted for actual mean of X)
  true_gbec <- 3/2 * mean(X[, 1]) 
  
  # Estimate MSEs
  mse_sure <- calculate_mse(compute_sure_gbec, X, Y, lambda_seq, Sigma, true_gbec = true_gbec)
  mse_ridge <- calculate_mse(compute_ridge_gbec, X, Y, lambda = 0.1, Sigma, true_gbec = true_gbec)
  mse_bec <- tryCatch(
    {
      gamma_bec <- compute_bec(X, Y, Sigma = Sigma)
      mean((gamma_bec - rep(1, 3))^2) # Calculate MSE directly for gamma_bec
    },
    error = function(e) {
      warning(paste("Error computing BEC:", e$message))
      NA
    }
  )
  
  # Add results to the data frame
  results_mse <- rbind(results_mse, data.frame(t1 = t1, Method = "SURE", MSE = mse_sure))
  results_mse <- rbind(results_mse, data.frame(t1 = t1, Method = "Ridge", MSE = mse_ridge))
  results_mse <- rbind(results_mse, data.frame(t1 = t1, Method = "BEC", MSE = mse_bec)) # Added this line
}

# Create plot (highlighting the true value)
ggplot(results_mse, aes(x = t1, y = MSE, color = Method, linetype = Method)) +
  geom_line(size = 1, na.rm = TRUE) +
  labs(title = "Comparison of SURE, Ridge, and BEC (MSE)",
       x = expression(t[1]~Value), 
       y = "Mean Squared Error") +
  scale_color_manual(values = c("BEC" = "black", "SURE" = "blue", "Ridge" = "red")) +  
  scale_linetype_manual(values = c("BEC" = "solid", "SURE" = "dashed", "Ridge" = "dashed")) +  
  theme_minimal() 

