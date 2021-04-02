############################################################################################################
# This file contains the syntax for replicating the empirical example in
# Miocevic, M., & Golchi, S. (2021). Bayesian mediation analysis with power prior distributions. 
# Multivariate Behavior Research.
#
# The syntax calls the file med_wWeight_PROMIS.stan and uses the data set PROMISempex.csv.
# The example analysis contains three covariates:
# 1. gender (dummy variable coded so male is the reference category) 
# 2. race (dummy variable coded so white is the reference category) 
# 3. age (continuous covariate)
#
# If the example syntax is used for an analysis with a different number of covariates and/or
# covariates with different measurement scales, lines 33-35 and 55-57 in this syntax file 
# and 6-8, 12-17, and 33-34 in the .stan file will require modification.
# In the .stan file, the independent variable is denoted x, the mediator is denoted m, and the outcome is y
# The continuous covariate is denoted z, and w1 and w2 represent binary covariates.
# The current syntax file can be adapted to a new data set by changing variable names in lines 33-35 for the 
# analysis with no borrowing and lines 55-57 for the analysis with power priors.
############################################################################################################

# load required R packages
library(readr)
library(rstan)
library(dplyr)
library(ggplot2)

# Extracting observations from the current data set (site = Duke)

PROMISempex <- read_csv("LOCATION OF DATA SET/PROMISempex.csv")
df_current<-PROMISempex[PROMISempex$site=="Duke",]

## Create a data list for rstan
data_current = list(N = dim(df_current)[1], y = df_current$trouble_conc, 
                    a0 = df_current$weights, m = df_current$restless_sleep, x = df_current$pain,
                    z = df_current$age, w1 = df_current$female, w2 = df_current$nonwhite)

# Run the model with no borrowing from the historical data 
fit = stan(file="med_PROMIS.stan", 
           data=data_current, 
           warmup=1000, iter=5000, chains=3)

# Convergence diagnostics 
## Rhat
summary(fit)
## Trace plots
traceplot(fit)
# Visual examination of posteriors
stan_dens(fit, separate_chains = TRUE)
# Extract draws for the indirect effect
fit_ss = extract(fit)
med_eff_NP = fit_ss$med_eff

# Run the model with power priors based on previously computed a0 (weights) in the data set
## Create a data list for rstan
data_both  = list(N = dim(PROMISempex)[1], y = PROMISempex$trouble_conc, 
                    a0 = PROMISempex$weights, m = PROMISempex$restless_sleep, x = PROMISempex$pain,
                    z = PROMISempex$age, w1 = PROMISempex$female, w2 = PROMISempex$nonwhite)
# Run the model with no borrowing from the historical data 
fit_borrowing = stan(file="med_PROMIS.stan", 
           data=data_both, 
           warmup=1000, iter=5000, chains=3)

# Convergence diagnostics 
## Rhat
summary(fit_borrowing)
## Trace plots
traceplot(fit_borrowing)
# Visual examination of posteriors
stan_dens(fit_borrowing, separate_chains = TRUE)
# Extract draws for the indirect effect
fit_ssp = extract(fit_borrowing)
med_eff_PP = fit_ssp$med_eff

# Plot the indirect effects from the analyses with diffuse and power priors for comparison
df = data.frame(ME = c(med_eff_NP, med_eff_PP), 
                method = c(rep('NP', 12000), 
                           rep('PP', 12000)))

# Compute the posterior medians and 95% credibility intervals for ab
df_summary = df %>%
  group_by(method) %>%
  dplyr::summarize(ME_med = median(ME), MElow = quantile(ME, 0.025), MEupp = quantile(ME, 0.975))

# Plot the indirect effect from the analyses with diffuse priors (NP) and power priors (PP)
ggplot(df_summary, aes(x = ME_med, y = method, color = method)) + geom_point(size = 2) +
  geom_errorbarh(aes(xmax = MEupp, xmin = MElow), height = 0, size = 1) + 
  xlab('Indirect Effect') + 
  scale_color_brewer(palette = 'Set1') + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'none')
