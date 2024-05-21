# Load necessary libraries
library(MASS)
library(glmnet)
library(ggplot2)

# ... (Functions generate_data, compute_bec, compute_gbec_ridge, compute_gbec_sure remain the same) ...

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

# Calculate MSE for each method
mse_bec <- mean((gamma_bec - data$beta)^2)
mse_gbec_ridge <- mean((gamma_gbec_ridge - data$beta)^2)
mse_gbec_sure <- mean((gamma_gbec_sure - data$beta)^2)

# Comparison of methods (MSE)
comparison_mse <- data.frame(
  Method = c("BEC", "GBEC (Ridge)", "GBEC (SURE)"),
  MSE = c(mse_bec, mse_gbec_ridge, mse_gbec_sure)
)

# Remove NA rows (if BEC couldn't be calculated)
comparison_mse <- comparison_mse[complete.cases(comparison_mse), ]


# Create the bar plot with labels
ggplot(comparison_mse, aes(x = Method, y = MSE, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.4f", MSE)), vjust = -0.5) +
  labs(title = "Comparison of GBEC (ridge) and GBEC(Sure) Methods (MSE)",
       x = "Method",
       y = "Mean Squared Error") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
