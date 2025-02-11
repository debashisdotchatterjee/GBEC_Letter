###################################################################################################
# PACKAGES & SETUP
###################################################################################################
if(!require(sp)) install.packages('sp', dependencies=TRUE)
if(!require(gstat)) install.packages('gstat', dependencies=TRUE)
if(!require(sf)) install.packages('sf', dependencies=TRUE)
if(!require(stars)) install.packages('stars', dependencies=TRUE)
if(!require(spdep)) install.packages('spdep', dependencies=TRUE)
if(!require(Matrix)) install.packages('Matrix', dependencies=TRUE)
if(!require(ggplot2)) install.packages('ggplot2', dependencies=TRUE)
if(!require(dplyr)) install.packages('dplyr', dependencies=TRUE)
if(!require(tidyr)) install.packages('tidyr', dependencies=TRUE)
if(!require(reshape2)) install.packages('reshape2', dependencies=TRUE)
if(!require(knitr)) install.packages('knitr', dependencies=TRUE)
if(!require(kableExtra)) install.packages('kableExtra', dependencies=TRUE)

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

###################################################################################################
# CREATE OUTPUT DIRECTORY
###################################################################################################
output_dir <- 'meuse_output'
dir.create(output_dir, showWarnings = FALSE)

###################################################################################################
# DATA LOADING & INSPECTION
###################################################################################################
data(meuse)

meuse_sf <- st_as_sf(meuse, coords = c('x', 'y'), crs = 28992)
cat('=== Summary of Meuse Data (sf object) ===\n')
print(summary(meuse_sf))

metals <- c('cadmium', 'copper', 'lead', 'zinc')

###################################################################################################
# DATA PREPROCESSING
###################################################################################################
meuse_df <- as.data.frame(meuse_sf)

if(!all(metals %in% colnames(meuse_df))) {
  stop('Some metals are not present in the meuse dataset.')
}

metals_data <- meuse_df[, metals]

########################################################################################
# FORCING NEAR-COLLINEARITY
########################################################################################
set.seed(123)
n_obs <- nrow(metals_data)

original_cadmium <- metals_data[, 'cadmium']

# tiny random noises
noise_copper <- rnorm(n_obs, mean=0, sd=0.001)
noise_lead <- rnorm(n_obs, mean=0, sd=0.002)
noise_zinc <- rnorm(n_obs, mean=0, sd=0.003)

# Force metals to be near-linear combos
metals_data[, 'copper'] <- original_cadmium + noise_copper
metals_data[, 'lead']   <- 2*original_cadmium + 3*metals_data[, 'copper'] + noise_lead
metals_data[, 'zinc']   <- metals_data[, 'lead'] + metals_data[, 'copper'] + noise_zinc

# Standardize
metals_scaled <- scale(metals_data)
metals_scaled_df <- as.data.frame(metals_scaled)
colnames(metals_scaled_df) <- metals

###################################################################################################
# SPATIAL COVARIANCE MATRIX (DEMO)
###################################################################################################
coordinates(meuse) <- ~x+y

vgm_model <- variogram(cadmium ~ 1, data = meuse)
fit_vgm <- fit.variogram(vgm_model, model = vgm(1, 'Exp', 300, 1))

psill <- fit_vgm$psill[2]
range_ <- fit_vgm$range[2]
nugget <- fit_vgm$psill[1]

coords_mat <- cbind(meuse$x, meuse$y)
n_points <- nrow(coords_mat)

Sigma_mat <- matrix(0, nrow=n_points, ncol=n_points)
for(i in 1:n_points){
  for(j in 1:n_points){
    h_ij <- sqrt((coords_mat[i,1] - coords_mat[j,1])^2 + 
                   (coords_mat[i,2] - coords_mat[j,2])^2)
    if(i == j){
      Sigma_mat[i,j] <- psill + nugget
    } else {
      Sigma_mat[i,j] <- psill * exp(-h_ij / range_)
    }
  }
}

###################################################################################################
# SPATIAL PENALTY MATRIX (ADJACENCY)
###################################################################################################
nb_5 <- dnearneigh(coords_mat, 0, 500)
W_list <- nb2listw(nb_5, style='B', zero.policy=TRUE)
W <- listw2mat(W_list)
W_scaled <- W / max(W)

###################################################################################################
# GBEC-LIKE CONSTRUCTION ON HEAVILY COLLINEAR METALS
###################################################################################################
cov_metals <- cov(metals_scaled_df)
p <- ncol(cov_metals)
J_p <- rep(1, p)

det_cov <- det(cov_metals)
cat(sprintf('Determinant of the (modified) metals covariance matrix: %.6e\n', det_cov))

if(det_cov < 1e-10) {
  cat('Covariance matrix nearly singular => Using Moore-Penrose inverse.\n')
  library(MASS)
  Sigma_inv <- ginv(cov_metals)
} else {
  Sigma_inv <- solve(cov_metals)
}

# Classical BEC
c_val <- 1 / as.numeric(t(J_p) %*% Sigma_inv %*% J_p)
gamma_bec <- c_val * Sigma_inv %*% J_p

cat('\n=== Classical BEC Coefficients ===\n')
print(gamma_bec)

# Ridge-based GBEC
lambda_seq <- seq(0, 2, by=0.1)

gbec_ridge <- function(Sigma, lambda){
  p_ <- ncol(Sigma)
  I_p <- diag(p_)
  denom <- t(J_p) %*% solve(Sigma + lambda*I_p) %*% J_p
  c_ridge <- 1 / as.numeric(denom)
  gamma_r <- c_ridge * solve(Sigma + lambda*I_p) %*% J_p
  gamma_r
}

# eq-cov measure => sum((Sigma*gamma - mean(Sigma*gamma))^2)
mse_vals <- numeric(length(lambda_seq))

for(i in seq_along(lambda_seq)){
  gamma_temp <- gbec_ridge(cov_metals, lambda_seq[i])
  eq_cov <- cov_metals %*% gamma_temp
  eq_mean <- mean(eq_cov)
  mse_vals[i] <- sum((eq_cov - eq_mean)^2)
}

optimal_lambda <- lambda_seq[which.min(mse_vals)]
cat('Optimal lambda based on eq-cov measure:', optimal_lambda, '\n')

gamma_optimal <- gbec_ridge(cov_metals, optimal_lambda)
cat('\n=== Ridge-based GBEC Coefficients (optimal) ===\n')
print(gamma_optimal)

###################################################################################################
# PLOTS & VISUALIZATION
###################################################################################################

# 1) Variogram plot for cadmium
vario_plot <- ggplot(vgm_model, aes(x=dist, y=gamma)) +
  geom_point(color='blue') +
  xlab('Distance') +
  ylab('Semivariance') +
  ggtitle('Empirical Variogram for Cadmium') +
  theme_minimal()

vario_plot_file <- file.path(output_dir, 'variogram_cadmium.png')
ggsave(vario_plot_file, plot=vario_plot, width=6, height=4)
print(vario_plot)

# 2) MSE vs. Lambda
mse_df <- data.frame(lambda=lambda_seq, MSE=mse_vals)
mse_plot <- ggplot(mse_df, aes(x=lambda, y=MSE)) +
  geom_line(color='red', linewidth=1) +
  geom_point(color='black') +
  theme_minimal() +
  ggtitle('MSE-like measure vs. Lambda for Ridge-GBEC') +
  xlab(expression(lambda)) +
  ylab('Eq-Cov Deviation')

mse_plot_file <- file.path(output_dir, 'MSE_vs_Lambda.png')
ggsave(mse_plot_file, plot=mse_plot, width=6, height=4)
print(mse_plot)

# 3) Barplot: Classical BEC vs. Ridge-based GBEC
coeff_df <- data.frame(
  Method = rep(c('Classical BEC', 'Ridge-GBEC'), each=p),
  Metal = rep(colnames(metals_scaled_df), 2),
  Coeff = c(gamma_bec, gamma_optimal)
)

coeff_plot <- ggplot(coeff_df, aes(x=Metal, y=Coeff, fill=Method)) +
  geom_bar(stat='identity', position='dodge') +
  ggtitle('Comparison: BEC vs. Ridge-GBEC (Collinearity)') +
  theme_minimal()

coeff_plot_file <- file.path(output_dir, 'Coefficients_Comparison.png')
ggsave(coeff_plot_file, plot=coeff_plot, width=6, height=4)
print(coeff_plot)

###################################################################################################
# TABLE OUTPUT
###################################################################################################
coef_table <- data.frame(
  Metal = colnames(metals_scaled_df),
  Classical_BEC = as.vector(gamma_bec),
  Ridge_GBEC = as.vector(gamma_optimal)
)

cat('\n=== Coefficient Comparison Table (Console) ===\n')
print(kable(coef_table, format='pipe'))

write.csv(coef_table, file.path(output_dir, 'CoefficientComparison.csv'), row.names=FALSE)

###################################################################################################
# CONCLUSION
###################################################################################################
cat('
===============================================================================
We artificially introduced near-collinearity among the four metals to make the 
covariance matrix nearly singular. This scenario forces Classical BEC to use 
the Moore-Penrose inverse, potentially giving unstable solutions. Ridge-based 
GBEC, however, finds an optimal lambda > 0 that yields a distinct solution, 
demonstrating the advantage of ridge shrinkage in ill-conditioned problems.

All relevant plots and tables are saved in:
meuse_output
===============================================================================
')
