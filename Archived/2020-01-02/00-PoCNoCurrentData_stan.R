rm(list=ls())

library(rstan)
library(RBesT)
library(ggplot2)
source("00-BhattacharyyaDistance.R")

dt = crohn
sigma = 88
dt$se_yh = sigma/sqrt(dt$n)
#historical control data
yh = dt$y
nh = dt$n
#control arm hyper-parameters
prec_mu_c = 1/10000; se_mu_c = sqrt(1/prec_mu_c)
HNscale_c = 44
n.chains = 3

##########################################################
##########################################################
##########################################################
stan_data = with(dt,list(H = length(study),
                         y = y,
                         n = n,
                         sigma = sigma,
                         HNscale_c = HNscale_c,
                         se_mu_c = se_mu_c)
)

stan_obj <- stan(
  file = "model_normal.stan",  # Stan program
  data = stan_data,    # named list of data
  chains = n.chains,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1,              # number of cores (could use one per chain)
  refresh = 0             # no progress shown
)


##########################################################
##########################################################
##########################################################

options(RBesT.MC.control=list(adapt_delta=0.999))
map_mcmc <- gMAP(cbind(y, se_yh) ~ 1 | study, 
                 weight = n,
                 data=dt,
                 family=gaussian,
                 beta.prior=cbind(0, sqrt(1/prec_mu_c)),
                 tau.dist="HalfNormal",tau.prior=cbind(0,HNscale_c),
                 chains = n.chains)

#approximate the MAP
map_hat <- mixfit(map_mcmc, Nc = 4)
# plot(map_hat)

postdt_rbest = rmix(mix = map_hat,n = length(postdt_jags))



#Bhattacharyya distance
bhatt.coeff(postdt_rbest,postdt_jags)

tt1 = data.frame(Method = "JAGS",Sample = postdt_jags)
tt2 = data.frame(Method = "RBesT",Sample = postdt_rbest)

plotdt = rbind(tt1,tt2)

ggplot(plotdt, aes(x=Sample, fill=Method)) + geom_density(alpha=.3)

summary(postdt_rbest)
summary(postdt_jags)
