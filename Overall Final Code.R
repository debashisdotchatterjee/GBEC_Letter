# Load necessary libraries
library(MASS)
library(glmnet)
library(ggplot2)

# Function to generate simulated data
generate_data <- function(n, p, sigma_noise) {
  X <- mvrnorm(n, rep(0, p), diag(p))
  beta <- runif(p, -1, 1)
  f <- rnorm(n)
  epsilon <- rnorm(n, 0, sigma_noise)
  Y <- X %*% beta + f + epsilon
  return(list(X = X, Y = Y, beta = beta))
}

# Function to compute BEC using pseudo-inverse (with error handling and improved inversion)
compute_bec <- function(X, Y, tol = 1e-7) {
  Sigma <- cov(X)
  
  # Ensure Sigma is a matrix (not a dataframe)
  Sigma <- as.matrix(Sigma)
  
  # Ensure J is a column vector and has the correct dimensions
  J <- matrix(1, ncol(Sigma), 1)
  
  # Check if Sigma is invertible with tolerance
  if (rcond(Sigma) < tol) {
    warning("Covariance matrix is ill-conditioned, BEC results may be unreliable.")
  }
  
  # Calculate pseudo-inverse using SVD with tolerance
  Sigma_pseudo_inv <- ginv(Sigma, tol = tol) 
  
  # Explicitly ensure dimensions match (with rounding)
  if (round(ncol(Sigma_pseudo_inv)) != round(nrow(J))) {
    stop("Dimension mismatch between Sigma_pseudo_inv and J (even after rounding).")
  }
  
  # Use tryCatch to handle potential errors in c_opt calculation
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


# Function to compute GBEC using ridge regression
compute_gbec_ridge <- function(X, Y, lambda) {
  model <- glmnet(X, Y, alpha = 0, lambda = lambda, intercept = FALSE)
  gamma <- as.vector(coef(model)[-1])  # Remove intercept coefficient
  return(gamma)
}

# Function to compute GBEC using SURE (Corrected)
compute_gbec_sure <- function(X, Y) {
  n <- nrow(X)
  p <- ncol(X)
  
  model <- glmnet(X, Y, alpha = 0, lambda = seq(0.01, 1, length = 100), intercept = FALSE)
  Y_matrix <- as.matrix(Y)
  predictions <- predict(model, X)
  residuals_squared <- sweep(predictions, 1, Y_matrix, FUN = "-")^2
  df <- model$df
  sigma2_hat <- colMeans(residuals_squared) / (1 - df / n)
  sure <- colSums(residuals_squared) + 2 * sigma2_hat * df - n * sigma2_hat
  best_lambda_index <- which.min(sure)
  best_lambda <- model$lambda[best_lambda_index]
  
  if (is.null(best_lambda)) {
    stop("Could not determine the optimal lambda. SURE calculation failed.")
  }
  
  gamma <- as.vector(coef(model, s = best_lambda)[-1])
  return(gamma)
}

# Simulation settings
set.seed(123)
n <- 100
p <- 50
sigma_noise <- 0.5

# Generate data
data <- generate_data(n, p, sigma_noise)
X <- data$X
Y <- data$Y

# Compute BEC and GBEC with error handling
gamma_bec <- tryCatch(
  compute_bec(X, Y),
  error = function(e) {
    warning(paste("Error computing BEC:", e$message))
    rep(NA, p) # Return a vector of NAs if there's an error
  }
)
gamma_gbec_ridge <- compute_gbec_ridge(X, Y, lambda = 0.1)
gamma_gbec_sure <- compute_gbec_sure(X, Y)

# Calculate true L1 and L2 norm of beta
true_l1_norm <- sum(abs(data$beta)) 
true_l2_norm <- norm(data$beta, "2")

# Calculate MSE for each method
mse_bec <- ifelse(all(is.na(gamma_bec)), NA, mean((gamma_bec - data$beta)^2))
mse_gbec_ridge <- mean((gamma_gbec_ridge - data$beta)^2)
mse_gbec_sure <- mean((gamma_gbec_sure - data$beta)^2)

# Comparison of methods (L1 Norm)
comparison_l1 <- data.frame(
  Method = c("BEC", "GBEC (Ridge)", "GBEC (SURE)", "True"),
  Norm = c(ifelse(all(is.na(gamma_bec)), NA, sum(abs(gamma_bec))), sum(abs(gamma_gbec_ridge)), sum(abs(gamma_gbec_sure)), true_l1_norm)
)

# Remove any NA values before plotting
comparison_l1 <- comparison_l1[complete.cases(comparison_l1), ]

# Comparison of methods (MSE)
comparison_mse <- data.frame(
  Method = c("BEC", "GBEC (Ridge)", "GBEC (SURE)"),
  MSE = c(mse_bec, mse_gbec_ridge, mse_gbec_sure)
)

# Remove NA rows (if BEC couldn't be calculated)
comparison_mse <- comparison_mse[complete.cases(comparison_mse), ]

# Create the bar plot with labels (L1 Norm)
ggplot(comparison_l1, aes(x = Method, y = Norm, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Norm)), vjust = -0.5) +
  geom_hline(yintercept = true_l1_norm, linetype = "dashed", color = "red") + 
  labs(title = "Comparison of BEC and GBEC Methods (L1 Norm)",
       x = "Method",
       y = "L1 Norm of Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create the bar plot with labels (MSE)
ggplot(comparison_mse, aes(x = Method, y = MSE, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.4f", MSE)), vjust = -0.5) +
  labs(title = "Comparison of BEC and GBEC Methods(Ridge and SURE) (MSE)",
       x = "Method",
       y = "Mean Squared Error") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#################################

# Simulation settings
set.seed(123)
n <- 100
p <- 50
sigma_noise <- 0.5

# Generate data
data <- generate_data(n, p, sigma_noise)
X <- data$X
Y <- data$Y

# Compute BEC and GBEC with error handling
gamma_bec <- tryCatch(
  compute_bec(X, Y),
  error = function(e) {
    warning(paste("Error computing BEC:", e$message))
    rep(NA, p) # Return a vector of NAs if there's an error
  }
)
gamma_gbec_ridge <- compute_gbec_ridge(X, Y, lambda = 0.1)
gamma_gbec_sure <- compute_gbec_sure(X, Y)

# Calculate true L1 and L2 norm of beta
true_l1_norm <- sum(abs(data$beta))
true_l2_norm <- norm(data$beta, "2") 

# Calculate MSE for each method
mse_bec <- ifelse(all(is.na(gamma_bec)), NA, mean((gamma_bec - data$beta)^2))
mse_gbec_ridge <- mean((gamma_gbec_ridge - data$beta)^2)
mse_gbec_sure <- mean((gamma_gbec_sure - data$beta)^2)

# Comparison of methods (L1 Norm)
comparison_l1 <- data.frame(
  Method = c("BEC", "GBEC (Ridge)", "GBEC (SURE)", "True"),
  Norm = c(ifelse(all(is.na(gamma_bec)), NA, sum(abs(gamma_bec))), sum(abs(gamma_gbec_ridge)), sum(abs(gamma_gbec_sure)), true_l1_norm)
)

# Remove any NA values before plotting
comparison_l1 <- comparison_l1[complete.cases(comparison_l1), ]

# Comparison of methods (L2 Norm)
comparison_l2 <- data.frame(
  Method = c("BEC", "GBEC (Ridge)", "GBEC (SURE)", "True"),
  Norm = c(ifelse(all(is.na(gamma_bec)), NA, norm(gamma_bec, "2")), norm(gamma_gbec_ridge, "2"), norm(gamma_gbec_sure, "2"), true_l2_norm)
)

# Remove any NA values before plotting
comparison_l2 <- comparison_l2[complete.cases(comparison_l2), ]

# Comparison of methods (MSE)
comparison_mse <- data.frame(
  Method = c("BEC", "GBEC (Ridge)", "GBEC (SURE)"),
  MSE = c(mse_bec, mse_gbec_ridge, mse_gbec_sure)
)

# Remove NA rows (if BEC couldn't be calculated)
comparison_mse <- comparison_mse[complete.cases(comparison_mse), ]

# Create the bar plot with labels (L1 Norm)
ggplot(comparison_l1, aes(x = Method, y = Norm, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Norm)), vjust = -0.5) +
  geom_hline(yintercept = true_l1_norm, linetype = "dashed", color = "black") + 
  labs(title = "Comparison of BEC and GBEC Methods (L1 Norm)",
       x = "Method",
       y = "L1 Norm of Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create the bar plot with labels (L2 Norm)
ggplot(comparison_l2, aes(x = Method, y = Norm, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Norm)), vjust = -0.5) +
  geom_hline(yintercept = true_l2_norm, linetype = "dashed", color = "black") + 
  labs(title = "Comparison of BEC and GBEC Methods (L2 Norm)",
       x = "Method",
       y = "L2 Norm of Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create the bar plot with labels (MSE)
ggplot(comparison_mse, aes(x = Method, y = MSE, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.4f", MSE)), vjust = -0.5) +
  labs(title = "Comparison of BEC and GBEC Methods (MSE)",
       x = "Method",
       y = "Mean Squared Error") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
