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
#current control data
y_c = -76.01
se_yc = 21.93
n_c = 20

#control arm hyper-parameters
se_mu_c = 100
HNscale_c = 44

#current treatment data/parameters
se_mu_t = 100
HNscale_t = 10000
y_t = -29.15
se_yt = 16.69
n_t = 39

#MCMC control parameters
n.chains = 3
nsim = 10

distance = rep(NA,length=nsim)
for(s in 1:nsim){
  print(s)
  ##########################################################
  ##########################################################
  ##########################################################
  jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, yh = yh, prec_yh = se_yh^-2, 
                           prec_mu_t = se_mu_t^-2, HNscale_t = HNscale_t, y_t = y_t, prec_yt = se_yt^-2)
  )
  jags_obj = jags(model.file = "Model_NormalHN.bugs",
                  data = jags_data,
                  parameters.to.save = c("mu_c"),
                  n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                  progress.bar = "none"
  )
  
  jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
  postdt_jags = data.frame(jags_auto$BUGSoutput$sims.matrix)$mu_c
  ##########################################################
  ##########################################################
  ##########################################################
  
  options(RBesT.MC.control=list(adapt_delta=0.999))
  map_mcmc <- gMAP(cbind(y, se_yh) ~ 1 | study, 
                   #weight = n,
                   data=dt,
                   family=gaussian,
                   beta.prior=cbind(0, se_mu_c),
                   tau.dist="HalfNormal",tau.prior=cbind(0,HNscale_c),
                   chains = n.chains)
  #approximate the MAP
  map_hat <- mixfit(map_mcmc, Nc = 4)
  #update with current control data
  post_mix = postmix(map_hat,m = y_c, se = se_yc)
  
  postdt_rbest = rmix(mix = post_mix,n = length(postdt_jags))
  
  #Bhattacharyya distance
  distance[s] = bhatt.coeff(postdt_rbest,postdt_jags)
  
  tt1 = data.frame(Method = "JAGS",Sample = postdt_jags)
  tt2 = data.frame(Method = "RBesT",Sample = postdt_rbest)
  plotdt = rbind(tt1,tt2)
  ggplot(plotdt, aes(x=Sample, fill=Method)) + geom_density(alpha=.3)
}

summary(distance)
