rm(list=ls())
set.seed(712)

library(R2jags)
library(RBesT)
library(ggplot2)

source("02-NormalSimSpec1.R")
mu_c <- -45

#derive the rMAP as it doesn't rely on current data
options(RBesT.MC.control=list(adapt_delta=0.999))
map_c_mcmc = gMAP(cbind(y, se_yh) ~ 1 | study, 
                  #weight = n,
                  data=dt,
                  family=gaussian,
                  beta.prior=cbind(0, se_mu_c),
                  tau.dist="HalfNormal",tau.prior=cbind(0,HNscale_c),
                  chains = n.chains)
#approximate the MAP
map_c_hat = automixfit(map_c_mcmc)
#robustify command needs a reference scale, although it won't be used
sigma(map_c_hat) = sigma

#non-informative priors for current treatment arm
prior_t = mixnorm(c(1,0,se_mu_t))

w_v = 0.5

#robustification
rmap_c = robustify(map_c_hat, weight = w_v, mean = robust_mean, sigma = robust_sd)

#simulate current control data
sample_c = rnorm(n_c, mu_c, sigma)
y_c = mean(sample_c)
se_yc = sd(sample_c)/sqrt(n_c)

#simulate current treatment data
sample_t = rnorm(n_t, mu_c + effsize, sigma)
y_t = mean(sample_t)
se_yt = sd(sample_t)/sqrt(n_t)

###############
#Alternative 1#
###############
jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, 
                         yh = yh, prec_yh = se_yh^-2,
                         prec_mu_t = se_mu_t^-2,                        
                         y_t = y_t, prec_yt = se_yt^-2,
                         w_v = w_v, robust_mean = robust_mean, robust_sd = robust_sd)
)
jags_obj = jags(model.file = "Model_NormalHNMix.bugs",
                data = jags_data,
                parameters.to.save = c("mu_c","mu_t"),
                n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                progress.bar = "none"
)
jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
post_c_alt1 = tt$mu_c
post_t_alt1 = tt$mu_t

#posterior difference
post_diff_alt1 = post_t_alt1 - post_c_alt1
decision_alt1 = mean(post_diff_alt1 <= Qcut)   


###############
#Alternative 2#
###############
#original MAP model
jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, 
                         yh = yh, prec_yh = se_yh^-2,
                         prec_mu_t = se_mu_t^-2,                        
                         y_t = y_t, prec_yt = se_yt^-2,
                         w_v = 0, robust_mean = robust_mean, robust_sd = robust_sd)
)
jags_obj = jags(model.file = "Model_NormalHNMix.bugs",
                data = jags_data,
                parameters.to.save = c("mu_c","mu_t"),
                n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                progress.bar = "none"
)

jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
post_c_alt2_b = tt$mu_c
post_t_alt2 = tt$mu_t

#robust model with the vague prior
jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, 
                         yh = yh, prec_yh = se_yh^-2,
                         prec_mu_t = se_mu_t^-2,                        
                         y_t = y_t, prec_yt = se_yt^-2,
                         w_v = 1, robust_mean = robust_mean, robust_sd = robust_sd)
)
jags_obj = jags(model.file = "Model_NormalHNMix.bugs",
                data = jags_data,
                parameters.to.save = c("mu_c","mu_t"),
                n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                progress.bar = "none"
)

jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
post_c_alt2_nb = data.frame(jags_auto$BUGSoutput$sims.matrix)$mu_c

#model averaging
nb <- rbinom(length(post_diff_alt1),1,w_v)
post_c_alt2 <- (1 - nb) * post_c_alt2_b + nb * post_c_alt2_nb

#posterior difference
post_diff_alt2 = post_t_alt2 - post_c_alt2
decision_alt2 = mean(post_diff_alt2 <= Qcut)

######################
#Genuine MAP and rMAP#
######################
#get the posterior mixtures and random sample
postmix_map_c = postmix(map_c_hat,m = y_c, se = se_yc); post_c_map = rmix(mix = postmix_map_c, n = length(post_diff_alt1))
postmix_rmap_c = postmix(rmap_c,m = y_c, se = se_yc); post_c_rmap = rmix(mix = postmix_rmap_c, n = length(post_diff_alt1))
postmix_map_t = postmix(prior_t, m = y_t, se = se_yt); post_t_map = rmix(mix = postmix_map_t, n = length(post_diff_alt1))

post_diff_map = post_t_map - post_c_map
decision_map = mean(post_diff_map <= Qcut)

post_diff_rmap = post_t_map - post_c_rmap
decision_rmap = mean(post_diff_rmap <= Qcut)


posteriordt <- rbind(
  tibble(Method = "MAP", Sample = post_c_map), 
  tibble(Method = "Method A", Sample = post_c_rmap), 
  tibble(Method = "Method B", Sample = post_c_alt1), 
  tibble(Method = "Method C", Sample = post_c_alt2)
)

posteriordt %>% 
  ggplot(aes(x = Sample, fill = Method)) + geom_density(alpha = 0.4) + labs(fill = "Method") + 
  scale_x_continuous(breaks = seq(-90,-10,10), name = "Current Control Mean mu_C") +
  theme(axis.title = element_text(face="bold",size=20),
        axis.text = element_text(size=16),
        legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.position = c(0.2, 0.5)
  )

round(c(median(post_c_map), median(post_c_rmap), median(post_c_alt1), median(post_c_alt2)), 2)
c(decision_map, decision_rmap, decision_alt1, decision_alt2)
