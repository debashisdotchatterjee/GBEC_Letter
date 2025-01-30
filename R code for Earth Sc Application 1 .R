###################################################################################################
# PACKAGES & SETUP
###################################################################################################
if(!require(sp)) install.packages("sp", dependencies=TRUE)
if(!require(gstat)) install.packages("gstat", dependencies=TRUE)
if(!require(sf)) install.packages("sf", dependencies=TRUE)
if(!require(stars)) install.packages("stars", dependencies=TRUE)
if(!require(spdep)) install.packages("spdep", dependencies=TRUE)
if(!require(Matrix)) install.packages("Matrix", dependencies=TRUE)
if(!require(ggplot2)) install.packages("ggplot2", dependencies=TRUE)
if(!require(dplyr)) install.packages("dplyr", dependencies=TRUE)
if(!require(tidyr)) install.packages("tidyr", dependencies=TRUE)
if(!require(reshape2)) install.packages("reshape2", dependencies=TRUE)
if(!require(knitr)) install.packages("knitr", dependencies=TRUE)
if(!require(kableExtra)) install.packages("kableExtra", dependencies=TRUE)
#install.packages("spdep", dependencies=TRUE)
#install.packages("gstat", dependencies=TRUE)

library(sp)
library(gstat)
library(sf)
library(stars)
library(spdep)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(knitr)
library(kableExtra)
library(spdep)
library(gstat)

# Create an output directory for saving plots and tables
output_dir <- "meuse_output"
dir.create(output_dir, showWarnings = FALSE)

###################################################################################################
# DATA LOADING & INSPECTION
###################################################################################################

# Load the meuse dataset
data(meuse)

# Convert the 'meuse' data to an sf (simple features) object for coordinate handling
meuse_sf <- st_as_sf(meuse, coords = c("x", "y"), crs = 28992)

cat("=== Summary of Meuse Data (sf object) ===\n")
print(summary(meuse_sf))

# We'll focus on heavy metal concentrations: cadmium, copper, lead, zinc
metals <- c("cadmium", "copper", "lead", "zinc")

###################################################################################################
# DATA PREPROCESSING
###################################################################################################

# Convert sf object to a data frame for manipulation
meuse_df <- as.data.frame(meuse_sf)

# Ensure these metal names are present
if(!all(metals %in% colnames(meuse_df))) {
  stop("Some metals are not present in the meuse dataset.")
}

# Subset the metals
metals_data <- meuse_df[, metals]

# Standardize the metal variables
metals_scaled <- scale(metals_data)  # returns a matrix
metals_scaled_df <- as.data.frame(metals_scaled)
colnames(metals_scaled_df) <- metals

###################################################################################################
# SPATIAL COVARIANCE MATRIX (DEMONSTRATION)
###################################################################################################

# Convert back to 'sp' style (for gstat variogram approach)
coordinates(meuse) <- ~x+y

# Compute an empirical variogram for cadmium
vgm_model <- variogram(cadmium ~ 1, data = meuse)
# Fit an exponential variogram model
fit_vgm <- fit.variogram(vgm_model, model = vgm(1, "Exp", 300, 1))

# Extract parameters from the fitted model
psill <- fit_vgm$psill[2]  # partial sill
range_ <- fit_vgm$range[2]
nugget <- fit_vgm$psill[1]

# Distances among sample points
coords_mat <- cbind(meuse$x, meuse$y)
n_points <- nrow(coords_mat)

# Build a simplified approximate covariance matrix:
# Cov(X_i, X_j) = psill * exp(-h/range_), plus nugget on diagonal for demonstration.
Sigma_mat <- matrix(0, nrow = n_points, ncol = n_points)
for(i in 1:n_points){
  for(j in 1:n_points){
    h_ij <- sqrt((coords_mat[i,1] - coords_mat[j,1])^2 + 
                   (coords_mat[i,2] - coords_mat[j,2])^2)
    if(i == j){
      # total sill
      Sigma_mat[i,j] <- psill + nugget
    } else {
      Sigma_mat[i,j] <- psill * exp(-h_ij / range_)
    }
  }
}

###################################################################################################
# SPATIAL PENALTY MATRIX (ADJACENCY-BASED)
###################################################################################################

# Construct a neighbor structure within 500m
nb_5 <- dnearneigh(coords_mat, 0, 500)
W_list <- nb2listw(nb_5, style = "B", zero.policy = TRUE)
W <- listw2mat(W_list)
# Scale adjacency matrix for demonstration
W_scaled <- W / max(W)

###################################################################################################
# GBEC-LIKE CONSTRUCTION ON METAL FEATURES
###################################################################################################
# The dimension here is p = 4 (cadmium, copper, lead, zinc). This is to show how to handle
# near-singularity in a small p scenario.

cov_metals <- cov(metals_scaled_df)  # 4x4
p <- ncol(cov_metals)
J_p <- rep(1, p)

# If covariance is nearly singular, use Moore-Penrose inverse; else standard inverse.
det_cov <- det(cov_metals)
if(det_cov < 1e-10) {
  cat("Covariance matrix nearly singular. Using Moore-Penrose inverse.\n")
  library(MASS)
  Sigma_inv <- ginv(cov_metals)
} else {
  Sigma_inv <- solve(cov_metals)
}

# Classical BEC
c_val <- 1 / as.numeric(t(J_p) %*% Sigma_inv %*% J_p)
gamma_bec <- c_val * Sigma_inv %*% J_p

cat("\n=== Classical BEC Coefficients ===\n")
print(gamma_bec)

# Ridge-based approach for stable GBEC
lambda_seq <- seq(0, 2, by = 0.1)

gbec_ridge <- function(Sigma, lambda){
  p_ <- ncol(Sigma)
  I_p <- diag(p_)
  denom <- t(J_p) %*% solve(Sigma + lambda * I_p) %*% J_p
  c_ridge <- 1 / as.numeric(denom)
  gamma_r <- c_ridge * solve(Sigma + lambda * I_p) %*% J_p
  return(gamma_r)
}

# We'll define an eq-cov measure: how far (Sigma*gamma) is from a constant vector.
# Minimizing sum((Sigma*gamma - mean(Sigma*gamma))^2) as a demonstration.
mse_vals <- numeric(length(lambda_seq))

for(i in seq_along(lambda_seq)){
  gamma_temp <- gbec_ridge(cov_metals, lambda_seq[i])
  eq_cov <- cov_metals %*% gamma_temp
  eq_mean <- mean(eq_cov)
  mse_vals[i] <- sum((eq_cov - eq_mean)^2)
}

optimal_lambda <- lambda_seq[which.min(mse_vals)]
cat("Optimal lambda based on eq-cov measure:", optimal_lambda, "\n")

gamma_optimal <- gbec_ridge(cov_metals, optimal_lambda)
cat("\n=== Ridge-based GBEC Coefficients (optimal) ===\n")
print(gamma_optimal)

###################################################################################################
# PLOTS & VISUALIZATION
###################################################################################################

# 1) Plot the fitted variogram for cadmium
vario_plot <- ggplot(vgm_model, aes(x = dist, y = gamma)) +
  geom_point(color = "blue") +
  xlab("Distance") +
  ylab("Semivariance") +
  ggtitle("Empirical Variogram for Cadmium") +
  theme_minimal()

vario_plot_file <- file.path(output_dir, "variogram_cadmium.png")
ggsave(vario_plot_file, plot = vario_plot, width = 6, height = 4)
print(vario_plot)

# 2) MSE-like measure vs. Lambda for ridge-based GBEC
mse_df <- data.frame(lambda = lambda_seq, MSE = mse_vals)
mse_plot <- ggplot(mse_df, aes(x = lambda, y = MSE)) +
  geom_line(color = "red", size = 1) +
  geom_point(color = "black") +
  theme_minimal() +
  ggtitle("MSE-like measure vs. Lambda for Ridge-GBEC") +
  xlab(expression(lambda)) +
  ylab("Eq-Cov Deviation")

mse_plot_file <- file.path(output_dir, "MSE_vs_Lambda.png")
ggsave(mse_plot_file, plot = mse_plot, width = 6, height = 4)
print(mse_plot)

# 3) Barplot comparing Classical BEC vs. Ridge-GBEC coefficients
coeff_df <- data.frame(
  Method = rep(c("Classical BEC", "Ridge-GBEC"), each = p),
  Metal = rep(metals, 2),
  Coeff = c(gamma_bec, gamma_optimal)
)

coeff_plot <- ggplot(coeff_df, aes(x = Metal, y = Coeff, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Comparison of BEC vs. Ridge-GBEC Coefficients") +
  theme_minimal()

coeff_plot_file <- file.path(output_dir, "Coefficients_Comparison.png")
ggsave(coeff_plot_file, plot = coeff_plot, width = 6, height = 4)
print(coeff_plot)

###################################################################################################
# TABLE OUTPUT
###################################################################################################

coef_table <- data.frame(
  Metal = metals,
  Classical_BEC = as.vector(gamma_bec),
  Ridge_GBEC = as.vector(gamma_optimal)
)

# Show in console
cat("\n=== Coefficient Comparison Table (Console) ===\n")
print(kable(coef_table, format = "pipe"))

# Also save table to CSV
write.csv(coef_table, file.path(output_dir, "CoefficientComparison.csv"), row.names = FALSE)

###################################################################################################
# CONCLUSION
###################################################################################################
cat("
===============================================================================
In this script, we:
1. Demonstrated an approximate spatial covariance matrix using a fitted variogram
   on the cadmium data from Meuse.
2. Implemented a classical BEC and a ridge-based GBEC on the covariance of
   standardized heavy metal concentrations to show how ridge shrinkage can help
   stabilize solutions.
3. Derived an 'eq-cov measure' as a toy example to pick an optimal lambda.
4. Produced multiple plots for visual inspection and saved them into the
   'meuse_output' directory.

We  can further refine:
- The spatial penalty W and incorporate it directly into the (Sigma + lambda*W)
  framework if your goal is to penalize correlations in a spatial sense.
- A SURE-based approach for dynamic selection of lambda.
- More advanced or domain-specific metrics.

All relevant plots and tables have been saved in the directory:
meuse_output
===============================================================================
")
