# Load necessary libraries
library(readr)
library(dplyr)
library(ggplot2)
library(rjags)
library(rstan)
library(bayesplot)
library(coda)


# Read the data
library(readr)
#data_MZ_LF_iu <- read_csv("C:/Users/User/Desktop/Debashis 2024/Tropical Disease PLOS/Datasets/data-MZ-LF-iu.csv")
data_MZ_SCH_iu <- read_csv("C:/Users/User/Desktop/Debashis 2024/Tropical Disease PLOS/Datasets/data-MZ-SCH-iu.csv")
#data_MZ_STH_iu <- read_csv("C:/Users/User/Desktop/Debashis 2024/Tropical Disease PLOS/Datasets/data-MZ-STH-iu.csv")

#data <- data_MZ_LF_iu
data <- data_MZ_SCH_iu
#data <- data_MZ_STH_iu


# Display the first few rows of the data
head(data)

# Check the column names
colnames(data)

# Prepare data by filtering out rows with NA values in PopTot and PopTreat
data <- data %>%
  filter(!is.na(PopTot) & !is.na(PopTreat) & PopTot > 0 & PopTreat >= 0 & PopTreat <= PopTot) %>%
  mutate(prevalence = PopTreat / PopTot)

# Display summary statistics
summary(data)

# JAGS model with Jeffrey's prior
jeffreys_model <- "
model {
  for (i in 1:N) {
    cases[i] ~ dbin(theta[i], population[i])
    theta[i] ~ dbeta(0.5, 0.5)
  }
}
"

# Prepare data for JAGS
jags_data <- list(
  cases = data$PopTreat,
  population = data$PopTot,
  N = nrow(data)
)

# Parameters to monitor
parameters <- c("theta")

# Initial values
inits <- function() list(theta = runif(nrow(data)))

# Run the JAGS model
jags_fit <- jags.model(
  textConnection(jeffreys_model),
  data = jags_data,
  inits = inits,
  n.chains = 3,
  n.adapt = 1000
)

# Update and sample
update(jags_fit, 1000)
samples <- coda.samples(jags_fit, parameters, n.iter = 5000)

# Combine chains and convert to dataframe
combined_samples <- as.mcmc(do.call(rbind, samples))
samples_df <- as.data.frame(combined_samples)

# Summary of the samples
jags_summary <- summary(samples)
print(jags_summary)
xtable(jags_summary)
# Extract parameter names and print them
parameter_names <- varnames(samples)
print(parameter_names)

# Select the first few parameters for plotting
selected_parameters <- parameter_names[1:5]

# Plot posterior distributions for the first few parameters
jpeg("jags_posteriors.jpg")
mcmc_areas(samples, pars = selected_parameters, prob = 0.95) +
  ggtitle("Posterior distributions with 95% credible intervals for first few parameters")
dev.off()

# Trace plots for the first few parameters
jpeg("jags_trace.jpg")
mcmc_trace(samples, pars = selected_parameters) +
  ggtitle("Trace plots for first few parameters")
dev.off()

# Save JAGS summary to a CSV file
write.csv(as.data.frame(jags_summary$statistics), "jags_summary.csv")
xtable(head(as.data.frame(jags_summary$statistics)))
# Stan model with informative prior
informative_model <- "
data {
  int<lower=0> N;
  int<lower=0> cases[N];
  int<lower=0> population[N];
}
parameters {
  real<lower=0, upper=1> theta[N];
  real<lower=0> alpha;
  real<lower=0> beta;
}
model {
  alpha ~ gamma(2, 0.1);  // Prior for alpha
  beta ~ gamma(2, 0.1);   // Prior for beta
  for (i in 1:N) {
    cases[i] ~ binomial(population[i], theta[i]);
    theta[i] ~ beta(alpha, beta);
  }
}
"

# Prepare data for Stan
stan_data <- list(
  cases = data$PopTreat,
  population = data$PopTot,
  N = nrow(data)
)

# Fit the Stan model
fit <- stan(
  model_code = informative_model,
  data = stan_data,
  iter = 5000,
  chains = 3,
  seed = 123
)

# Summary of the fit
stan_summary <- summary(fit)
print(stan_summary)

# Extract samples
stan_samples <- extract(fit)

# Convert to array for plotting
stan_samples_array <- as.array(fit)

# Print the names of the parameters in the Stan model
stan_param_names <- dimnames(stan_samples_array)[[3]]
print(stan_param_names)
library(xtable)
# Select the first few parameters for plotting
selected_parameters <- grep("^theta\\[", stan_param_names, value = TRUE)[1:5]

# Plot posterior distributions for the first few parameters
jpeg("stan_posteriors.jpg")
mcmc_areas(stan_samples_array, pars = selected_parameters, prob = 0.95) +
  ggtitle("Posterior distributions with 95% credible intervals (Informative Prior)")
dev.off()

# Trace plots for the first few parameters
jpeg("stan_trace.jpg")
mcmc_trace(stan_samples_array, pars = selected_parameters) +
  ggtitle("Trace plots (Informative Prior)")
dev.off()

# Save Stan summary to a CSV file
stan_summary_df <- as.data.frame(stan_summary$summary)
write.csv(stan_summary_df, "stan_summary.csv")
xtable(stan_summary_df)
# Diagnostics for JAGS model
jpeg("jags_gelman.jpg")
gelman_diag <- gelman.diag(samples)
print(gelman_diag)
dev.off()

# Diagnostics for Stan model
jpeg("stan_gelman.jpg")
stan_trace(fit, pars = selected_parameters)
dev.off()

jpeg("stan_ac.jpg")
stan_ac(fit, pars = selected_parameters)
dev.off()

# Save Gelman-Rubin diagnostics for JAGS to a CSV file
write.csv(as.data.frame(gelman_diag$psrf), "jags_gelman.csv")

# Save effective sample size for JAGS and Stan
jpeg("jags_ess.jpg")
jags_ess <- effectiveSize(samples)
print(jags_ess)
dev.off()

####
# Convert Stan fit object to coda mcmc object
stan_coda <- As.mcmc.list(fit)

# Compute effective sample size
stan_ess <- effectiveSize(stan_coda)

# Save effective sample size to a CSV file
write.csv(as.data.frame(stan_ess), "stan_ess.csv")

# Optionally, plot the effective sample size
jpeg("stan_ess.jpg")
plot(stan_ess, main="Effective Sample Size for Stan Model")
dev.off()

########

# Save R-hat statistics for Stan to a CSV file
write.csv(stan_summary_df[, "Rhat"], "stan_rhat.csv")

