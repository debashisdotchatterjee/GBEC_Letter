# Load necessary libraries
library(glmnet)
library(ggplot2)
library(MASS)
library(Matrix)

# Function to generate Toeplitz matrix (with PD fix using nearPD)
generate_toeplitz <- function(t1, t2, t3) {
  T <- matrix(c(t1, t2, t3, t2, t1, t2, t3, t2, t1), nrow = 3, byrow = TRUE)
  # Ensure positive definiteness using nearPD
  if (isSymmetric(T) && all(eigen(T)$values > 0)) {
    return(T)
  } else {
    return(nearPD(T)$mat)  # Use nearPD if not positive definite
  }
}

# Function to compute SURE estimate of GBEC score (with singularity check)
# Function to compute SURE estimate of GBEC score (with singularity check)
compute_sure_gbec <- function(X, Y, lambda_seq) {
  n <- nrow(X)
  Sigma <- cov(X)
  model <- glmnet(X, Y, alpha = 0, lambda = lambda_seq, intercept = FALSE)
  df <- model$df
  sigma2_hat <- colMeans((Y - predict(model, X))^2) / (1 - df/n)
  sure <- colSums((Y - predict(model, X))^2) + 2 * sigma2_hat * df - n * sigma2_hat
  best_lambda <- model$lambda[which.min(sure)]
  
  # Check singularity (with higher tolerance for stability)
  if (rcond(cov(X) + best_lambda * diag(ncol(X))) < 1e-6) {
    return(NA)
  } 
  
  gamma_sure <- as.vector(coef(model, s = best_lambda)[-1])
  gbec_score <- as.numeric(t(gamma_sure) %*% t(X))
  return(gbec_score)
}


# Function to compute ridge estimate of GBEC score (with singularity check)
compute_ridge_gbec <- function(X, Y, lambda) {
  Sigma <- cov(X)
  # Check singularity (with higher tolerance for stability)
  if (rcond(cov(X) + lambda * diag(ncol(X))) < 1e-6) {
    return(NA)
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
    return(rep(NA, ncol(Sigma))) # Return NAs if there's an error
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
t2_values <- seq(0.1, 1.1, by = 0.05)  # Vary t2
t1 <- 1.0  # Fixed t1
t3 <- 0.9  # Fixed t3
lambda_seq <- seq(0.01, 1, length.out = 100)

# Store results
results_mse <- data.frame(t2 = numeric(), Method = character(), MSE = numeric())
true_gbec_values <- numeric(length(t2_values))

# Main Simulation Loop
for (i in 1:length(t2_values)) {
  t2 <- t2_values[i]
  Sigma <- generate_toeplitz(t1, t2, t3)
  true_gbec <- mean(diag(Sigma))
  true_gbec_values[i] <- true_gbec
  
  for (sim in 1:num_simulations) {
    X <- mvrnorm(n, rep(0, 3), Sigma)
    Y <- rowSums(X) + rnorm(n, sd = 0.5)
    
    # Calculate GBEC scores for SURE and Ridge
    gbec_sure <- compute_sure_gbec(X, Y, lambda_seq)
    gbec_ridge <- compute_ridge_gbec(X, Y, lambda = 0.1)  
    
    # Calculate MSEs for SURE and Ridge
    mse_sure <- ifelse(is.na(gbec_sure), NA, (gbec_sure - true_gbec)^2) 
    mse_ridge <- ifelse(is.na(gbec_ridge), NA, (gbec_ridge - true_gbec)^2) 
    
    results_mse <- rbind(results_mse, data.frame(t2 = t2, Method = "SURE", MSE = mse_sure))
    results_mse <- rbind(results_mse, data.frame(t2 = t2, Method = "Ridge", MSE = mse_ridge))
  }
  
  # Calculate MSE for BEC (outside of the num_simulations loop)
  gamma_bec <- compute_bec(X, Y, Sigma = Sigma)
  mse_bec <- ifelse(all(is.na(gamma_bec)), NA, mean((gamma_bec - rep(1, 3))^2)) # Calculate MSE directly for gamma_bec
  results_mse <- rbind(results_mse, data.frame(t2 = t2, Method = "BEC", MSE = mse_bec))
}


# Load the dplyr package
library(dplyr) 

# Compute average MSE for each method
results_mse <- results_mse %>%
  group_by(t2, Method) %>%
  summarize(MSE = mean(MSE, na.rm = TRUE), .groups = "drop")

# Create the plot
ggplot(results_mse, aes(x = t2, y = MSE, color = Method, linetype = Method)) +
  geom_line(size = 1, na.rm = TRUE) +
  labs(title = "Comparison of SURE, Ridge, and BEC (MSE) for Exchangeable Correlation Matrix",
       x = expression(t[2]~Value),  # Update x-axis label
       y = "Mean Squared Error") +
  scale_color_manual(values = c("BEC" = "black", "SURE" = "blue", "Ridge" = "red")) +
  scale_linetype_manual(values = c("BEC" = "solid", "SURE" = "dashed", "Ridge" = "dashed")) +
  theme_minimal() + 
  ylim(c(0, 0.6)) # For better visualization of the plot
