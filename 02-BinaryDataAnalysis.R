rm(list=ls())
set.seed(712)

library(R2jags)
library(RBesT)
library(ggplot2)

source("02-BinarySimSpec1.R")

#derive the rMAP as it doesn't rely on current data
options(RBesT.MC.control=list(adapt_delta=0.999))
map_c_mcmc = gMAP(cbind(r, n - r) ~ 1 | study, 
                  #weight = n,
                  data=dt,
                  family=binomial,
                  beta.prior=cbind(0, se_mu_c),
                  tau.dist="HalfNormal",tau.prior=cbind(0,HNscale_c),
                  chains = n.chains)
#approximate the MAP
map_c_hat = automixfit(map_c_mcmc)

#non-informative beta prior for current treatment arm
prior_t = mixbeta(c(1,a_t,b_t), param = "ab")

w_v <- 0.5
mu <- 0.25
#################################################################
#################################################################
#################################################################
#robustification, note that RBesT uses the mean/sample size-1 parametrization, see help of "robustify" for details
rmap_c = robustify(map_c_hat, weight = w_v, mean = a_c/(a_c+b_c), n = a_c+b_c- 1)

#simulate current control data
y_c = rbinom(1, n_c, mu)
#simulate current treatment data
y_t = rbinom(1, n_t, mu + effsize)

###############
#Alternative 1#
###############
jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, 
                         y_c = y_c, n_c = n_c,
                         y_h = dt$r, n_h = dt$n,
                         y_t = y_t, n_t = n_t, prec_t = se_mu_t^-2, 
                         w_v = w_v, robust_sd = robust_sd)
)
jags_obj = jags(model.file = "Model_BinaryHNMix.bugs",
                data = jags_data,
                parameters.to.save = c("p_c","p_t"),
                n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                progress.bar = "none"
)
jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
post_c_alt1 = tt$p_c
post_t_alt1 = tt$p_t

#posterior difference
post_diff_alt1 = post_t_alt1 - post_c_alt1
decision_alt1 = mean(post_diff_alt1 > Qcut)      

###############
#Alternative 2#
###############
#original MAP model
jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, 
                         y_c = y_c, n_c = n_c,
                         y_h = dt$r, n_h = dt$n,
                         y_t = y_t, n_t = n_t, prec_t = se_mu_t^-2, 
                         w_v = 0, robust_sd = robust_sd)
)
jags_obj = jags(model.file = "Model_BinaryHNMix.bugs",
                data = jags_data,
                parameters.to.save = c("p_c","p_t"),
                n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                progress.bar = "none"
)

jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
post_c_alt2_b = tt$p_c
post_t_alt2 = tt$p_t

#robust model with the vague prior
jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, 
                         y_c = y_c, n_c = n_c,
                         y_h = dt$r, n_h = dt$n,
                         y_t = y_t, n_t = n_t, prec_t = se_mu_t^-2, 
                         w_v = 1, robust_sd = robust_sd)
)
jags_obj = jags(model.file = "Model_BinaryHNMix.bugs",
                data = jags_data,
                parameters.to.save = c("p_c","p_t"),
                n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                progress.bar = "none"
)

jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
post_c_alt2_nb = data.frame(jags_auto$BUGSoutput$sims.matrix)$p_c

#model averaging
nb <- rbinom(length(post_diff_alt1),1,w_v)
post_c_alt2 <- (1 - nb) * post_c_alt2_b + nb * post_c_alt2_nb

#posterior difference
post_diff_alt2 = post_t_alt2 - post_c_alt2
decision_alt2 = mean(post_diff_alt2 > Qcut)

######################
#Genuine MAP and rMAP#
######################
#get the posterior mixtures and random sample
postmix_map_c = postmix(map_c_hat,r = y_c, n = n_c); post_c_map = rmix(mix = postmix_map_c, n = length(post_diff_alt1))
postmix_rmap_c = postmix(rmap_c,r = y_c, n = n_c); post_c_rmap = rmix(mix = postmix_rmap_c, n = length(post_diff_alt1))
postmix_map_t = postmix(prior_t, r = y_t, n = n_t); post_t_map = rmix(mix = postmix_map_t, n = length(post_diff_alt1))

post_diff_map = post_t_map - post_c_map
decision_map = mean(post_diff_map > Qcut)

post_diff_rmap = post_t_map - post_c_rmap
decision_rmap = mean(post_diff_rmap > Qcut)


posteriordt <- rbind(
  tibble(Method = "MAP", Sample = post_c_map), 
  tibble(Method = "Method A", Sample = post_c_rmap), 
  tibble(Method = "Method B", Sample = post_c_alt1), 
  tibble(Method = "Method C", Sample = post_c_alt2)
)

posteriordt %>% 
  ggplot(aes(x = Sample, linetype = Method, color = Method)) + geom_density(alpha = 0.4, size = 1.5) + labs(fill = "Method") + ylab("Kernel Density Estimate") +
  scale_x_continuous(breaks = seq(0, 1, 0.1), name = "Current Control Rate p_C") +
  theme(axis.title = element_text(face="bold",size=20),
        axis.text = element_text(size=16),
        legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.position = c(0.1, 0.5)
  )

round(c(median(post_c_map), median(post_c_rmap), median(post_c_alt1), median(post_c_alt2)), 4)
c(decision_map, decision_rmap, decision_alt1, decision_alt2)





#assemble the results
decisions <- c(decision_map, decision_rmap, decision_alt1, decision_alt2)
median_est <- c(median(post_c_map), median(post_c_rmap), median(post_c_alt1), median(post_c_alt2))

c(decisions, median_est, w_post)