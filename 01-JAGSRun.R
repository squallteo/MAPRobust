rm(list=ls())

library(R2jags)
library(RBesT)
library(ggplot2)
source("00-BhattacharyyaDistance.R")

#historical control data
yh = c(51,49,36,47)
se_yh = c(4,6,5,5); prec_yh = 1/se_yh^2
#current control data
y_c = 29
se_yc = 5; prec_yc = 1/se_yc^2
n_c = 20
#control arm hyper-parameters
prec_mu_c = 1/10000
HNscale_c = 44

#current treatment data


#MCMC control parameters
n.chains = 3

##########################################################
##########################################################
##########################################################
jags_data = list(prec_mu_c = prec_mu_c, prec_mu_t = 1/10000,
                 HNscale_c = HNscale_c, HNscale_t = 1/200,
                 y_c = y_c, prec_yc = prec_yc,
                 y_t = 76, prec_yt = 1/17,
                 yh = yh, prec_yh = prec_yh
                 )
jags_obj = jags(model.file = "model_normal_HN.txt",
                data = jags_data,
                parameters.to.save = c("theta_t","theta_c"),
                n.chains = n.chains, n.burnin = 2000, n.iter = 10000
                )

jags_auto = autojags(jags_obj,Rhat = 1.1,n.thin = 4,n.iter = 20000)
postdt_jags = data.frame(jags_auto$BUGSoutput$sims.matrix)
postdt_jags$post_diff = with(postdt_jags, theta_t - theta_c)
##########################################################
##########################################################
##########################################################
histdt = data.frame(yh = yh, se_yh = se_yh, study = 1:length(yh))

options(RBesT.MC.control=list(adapt_delta=0.999))
map_mcmc <- gMAP(cbind(yh, se_yh) ~ 1 | study, 
                 data=histdt,
                 family=gaussian,
                 beta.prior=cbind(0, sqrt(1/prec_mu_c)),
                 tau.dist="HalfNormal",tau.prior=cbind(0,HNscale_c),
                 chains = n.chains)

#approximate the MAP
map_hat <- mixfit(map_mcmc, Nc = 4)
#update approximated map with current trial data
post_c = postmix(map_hat, m=y_c, se=se_yc)

postdt_c_rbest = rmix(mix = post_c,n = nrow(post_c))

#Bhattacharyya distance
bhatt.coeff(postdt_c_rbest,postdt_jags$theta_c)


summary(postdt_c_rbest)
summary(postdt_jags$theta_c)
