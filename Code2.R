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

# Function to compute BEC using pseudo-inverse (with error handling)
compute_bec <- function(X, Y) {
  Sigma <- cov(X)
  
  # Ensure Sigma is a matrix (not a dataframe)
  Sigma <- as.matrix(Sigma)
  
  # Ensure J is a column vector and has the correct dimensions
  J <- matrix(1, ncol(Sigma), 1)
  
  # Check if Sigma is invertible
  if (rcond(Sigma) < .Machine$double.eps) {
    warning("Covariance matrix is ill-conditioned, BEC results may be unreliable.")
  }
  
  Sigma_pseudo_inv <- ginv(Sigma)
  
  # Explicitly ensure dimensions match
  if (ncol(Sigma_pseudo_inv) != nrow(J)) {
    stop("Dimension mismatch between Sigma_pseudo_inv and J.")
  }
  
  c_opt <- 1 / (t(J) %*% Sigma_pseudo_inv %*% J)
  gamma <- c_opt * Sigma_pseudo_inv %*% J
  return(drop(gamma))
}

# Function to compute GBEC using ridge regression
compute_gbec_ridge <- function(X, Y, lambda) {
  model <- glmnet(X, Y, alpha = 0, lambda = lambda, intercept = FALSE)
  gamma <- as.vector(coef(model)[-1])
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

# Comparison of methods (including BEC even if there's an error)
comparison <- data.frame(
  Method = c("BEC", "GBEC (Ridge)", "GBEC (SURE)"),
  Norm = c(norm(gamma_bec, "2"), norm(gamma_gbec_ridge, "2"), norm(gamma_gbec_sure, "2"))
)

# Remove any NA values before plotting
comparison <- comparison[complete.cases(comparison),]

# Create the bar plot with labels
ggplot(comparison, aes(x = Method, y = Norm, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Norm)), vjust = -0.5) +  # Add labels above the bars
  labs(title = "Comparison of BEC and GBEC Methods",
       x = "Method",
       y = "L2 Norm of Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

