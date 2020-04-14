library(rstan)
library(dplyr)
library(ggplot2)


df0 = read.csv('DANS_std.csv', header = T)
df = df0[1:200, c('z_critizied', 'z_impatient', 'z_ridicules', 'z_prettig', 'z_sympathiek')]
study = c(rep(1, 100), rep(0, 100))
Ns = length(study)
y = df$z_critizied
x = df$z_impatient
M = df$z_ridicules
W1 = df$z_prettig
z = df$z_sympathiek

y1 = y[study==1] 
x1 = x[study==1] 
M1 = M[study==1]
W11 = W1[study==1]
z1 = z[study==1]


mu = apply(df[1:100], 2, mean) # sample mean for current study
Sigma = cov(df[1:100]) # sample covariance matrix for current study
d = apply(df, 1, mahalanobis, center = mu, cov = ginv(Sigma), inverted = T) # Mahalanobis distance
w2 = 1- (d - min(d))/(max(d) - min(d)) # weights
eps = quantile(w2[study==1], p = .05) # truncation threshold
w2[w2<eps] =0 # truncating weights
w2[study==1] = 1 # current study weights are set equal to 1


data = list(N = length(y), y = y, w = w2, x = x, M = M, W1 = W1, z = z)
fit1 = stan(file='med_wWeight_DANS.stan', data=data, chains=0)

#############################################################################################
### No Prior

w_NP = rep(0, length(y))
w_NP[study==1] = 1
data = list(N = length(y), y = y, w = w_NP, x = x, M = M, W1 = W1, z = z)
fit = stan('fit'=fit1, 'data'=data, warmup=1000, iter=2000, chains=2)
fit_ss = extract(fit)
med_eff_NP = fit_ss$med_eff

#############################################################################################
### Individually Weighted Prior

data2 = list(N = length(y), y = y, w = w2, x = x, M = M, W1 = W1, z = z)
fit = stan('fit'=fit1, 'data'=data2, warmup=1000, iter=2000, chains=2)
fit_ss = extract(fit)
med_eff_IW = fit_ss$med_eff



df = data.frame(ME = c(med_eff_NP, med_eff_IW), 
                 method = c(rep('NP', 2000), 
                            rep('PP', 2000)))


df_summary = df %>%
  group_by(method) %>%
  dplyr::summarize(ME_med = median(ME), MElow = quantile(ME, 0.025), MEupp = quantile(ME, 0.975))

### plots ##################################################################################

ggplot(df_summary, aes(x = ME_med, y = method)) + geom_point(size = 2) +
  geom_errorbarh(aes(xmax = MEupp, xmin = MElow), height = 0, size = 1) + 
  xlab('Mediated Effect') + 
  scale_color_brewer(palette = 'Set1') + geom_vline(xintercept = 0.068) +
  theme(axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = 'bold'),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = 'bold'),
    legend.position = 'none')
ggsave(file = 'DANS_data.png')

dfw = data.frame(weight = w2[study==0])
ggplot(dfw, aes(x = weight)) + geom_histogram() + geom_vline(xintercept = eps) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'none')
ggsave(file = 'weights_DANS.png')

#mean(w2[study==0])
#sum(w2[study==0]>0)/length(w2[study==0])

