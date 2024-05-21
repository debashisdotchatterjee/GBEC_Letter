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

# Calculate true norm of beta
true_norm <- norm(data$beta, "2")

# Comparison of methods (including BEC and the true norm)
comparison <- data.frame(
  Method = c("BEC", "GBEC (Ridge)", "GBEC (SURE)", "True"),
  Norm = c(norm(gamma_bec, "2"), norm(gamma_gbec_ridge, "2"), norm(gamma_gbec_sure, "2"), true_norm)
)

# Remove any NA values before plotting
comparison <- comparison[complete.cases(comparison),]

# Create the bar plot with labels and the true norm as a horizontal line
ggplot(comparison, aes(x = Method, y = Norm, fill = Method)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f", Norm)), vjust = -0.5) +  # Add labels above the bars
  geom_hline(yintercept = true_norm, linetype = "dashed", color = "red") + # Add a dashed line for true norm
  labs(title = "Comparison of BEC and GBEC Methods (with True Norm)",
       x = "Method",
       y = "L2 Norm of Coefficients") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

