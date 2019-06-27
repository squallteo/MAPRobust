rm(list=ls())

library(R2jags)
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
prec_mu_c = 1/10000
HNscale_c = 44
n.chains = 3

nsim = 10

distance = rep(NA,length=nsim)
for(s in 1:nsim){
  print(s)
  ##########################################################
  ##########################################################
  ##########################################################
  jags_data = list(prec_mu_c = prec_mu_c, HNscale_c = HNscale_c, yh = yh, nh = nh, sigma = sigma)
  jags_obj = jags(model.file = "Model_NormalHN_NoCurrent.bugs",
                  data = jags_data,
                  parameters.to.save = c("mu_c"),
                  n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                  progress.bar = "none"
  )
  
  jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10)
  postdt_jags = data.frame(jags_auto$BUGSoutput$sims.matrix)$mu_c
  ##########################################################
  ##########################################################
  ##########################################################
  
  options(RBesT.MC.control=list(adapt_delta=0.999))
  map_mcmc <- gMAP(cbind(y, se_yh) ~ 1 | study, 
                   #weight = n,
                   data=dt,
                   family=gaussian,
                   beta.prior=cbind(0, sqrt(1/prec_mu_c)),
                   tau.dist="HalfNormal",tau.prior=cbind(0,HNscale_c),
                   chains = n.chains)
  
  #approximate the MAP
  map_hat <- mixfit(map_mcmc, Nc = 4)
  postdt_rbest = rmix(mix = map_hat,n = length(postdt_jags))
  
  #Bhattacharyya distance
  distance[s] = bhatt.coeff(postdt_rbest,postdt_jags)
  
  tt1 = data.frame(Method = "JAGS",Sample = postdt_jags)
  tt2 = data.frame(Method = "RBesT",Sample = postdt_rbest)
  plotdt = rbind(tt1,tt2)
  ggplot(plotdt, aes(x=Sample, fill=Method)) + geom_density(alpha=.3)
}

