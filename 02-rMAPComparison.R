rm(list=ls())

library(R2jags)
library(RBesT)
library(ggplot2)
source("00-BhattacharyyaDistance.R")


#historical control data
dt = crohn
sigma = 88
dt$se_yh = sigma/sqrt(dt$n)
yh = dt$y
nh = dt$n

#meta analysis of historical data
# with(dt, meta::metamean(n = n, mean = y, sd = rep(sigma,length(study))))
#decide to consider true current mean from -60 to -40

#current control data
n_c = 20 #current control sample size
muvec_c = seq(-60, -40, by=2) #current control mean vector

#control arm hyper-parameters
se_mu_c = 100
HNscale_c = 44

#robustification parameter: the weight associated with the robust part
w_r = 0.2
robust_sd = 200

#current treatment data/parameters
se_mu_t = 100
HNscale_t = 10000
n_t = 40
effsize = 15

#MCMC control parameters
n.chains = 3
nsim = 10

#initialization at the control mean level


for(m in 1:length(muvec_c)){
  #simulate current control data
  sample_c = rnorm(n_c, muvec_c[m], sigma)
  y_c = mean(sample_c)
  se_yc = sd(sample_c)/sqrt(n_c)
  
  #simulate current treatment data
  sample_t = rnorm(n_t, muvec_c[m] + effsize, sigma)
  y_t = mean(sample_t)
  se_yt = sd(sample_t)/sqrt(n_t)
  
  #initialization at the simulation level
  distance = dur_jags = dur_rbest = rep(NA,length=nsim)
  
}


for(s in 1:nsim){
  print(s)
  ##########################################################
  ##########################################################
  ##########################################################
  start = Sys.time()
  #original MAP model
  jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, yh = yh, prec_yh = se_yh^-2, 
                           prec_mu_t = se_mu_t^-2, HNscale_t = HNscale_t, y_t = y_t, prec_yt = se_yt^-2)
  )
  jags_obj = jags(model.file = "Model_NormalHN.bugs",
                  data = jags_data,
                  parameters.to.save = c("mu_c","mu_t"),
                  n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                  progress.bar = "none"
  )
  
  jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
  tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
  post_c_jags_b = tt$mu_c
  post_t_jags = tt$mu_t
  
  #robust model with a fixed tau
  jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, prec_c = robust_sd^-2, y_c = y_c, prec_yc = se_yc^-2, yh = yh, prec_yh = se_yh^-2, 
                           prec_mu_t = se_mu_t^-2, HNscale_t = HNscale_t, y_t = y_t, prec_yt = se_yt^-2)
  )
  jags_obj = jags(model.file = "Model_NormalFixedTau.bugs",
                  data = jags_data,
                  parameters.to.save = c("mu_c"),
                  n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                  progress.bar = "none"
  )
  
  jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
  post_c_jags_nb = data.frame(jags_auto$BUGSoutput$sims.matrix)$mu_c
  
  #model averaging
  post_c_jags = (1 - w_r) * post_c_jags_b + w_r * post_c_jags_nb
  
  end = Sys.time()
  dur_jags[s] = as.numeric(end - start)
  
  #posterior difference
  post_diff_jags = post_c_jags - post_t_jags
  
  
  ##########################################################
  ##########################################################
  ##########################################################
  start = Sys.time()
  
  options(RBesT.MC.control=list(adapt_delta=0.999))
  map_mcmc = gMAP(cbind(y, se_yh) ~ 1 | study, 
                   #weight = n,
                   data=dt,
                   family=gaussian,
                   beta.prior=cbind(0, se_mu_c),
                   tau.dist="HalfNormal",tau.prior=cbind(0,HNscale_c),
                   chains = n.chains)
  #approximate the MAP
  map_hat = mixfit(map_mcmc, Nc = 3)
  #robustify command needs a reference scale, although it won't be used
  sigma(map_hat) = sigma
  #robustification
  rmap = robustify(map_hat, weight = w_r, mean = 0, sigma = robust_sd)
  
  #update with current control data
  post_mix = postmix(rmap,m = y_c, se = se_yc)
  
  postdt_rbest = rmix(mix = post_mix,n = length(postdt_jags))
  
  end = Sys.time()
  
  dur_rbest[s] = as.numeric(end - start)
  
  #Bhattacharyya distance
  distance[s] = bhatt.coeff(postdt_rbest,postdt_jags)
  
  # tt1 = data.frame(Method = "JAGS",Sample = postdt_jags)
  # tt2 = data.frame(Method = "RBesT",Sample = postdt_rbest)
  # plotdt = rbind(tt1,tt2)
  # p = ggplot(plotdt, aes(x=Sample, fill=Method)) + geom_density(alpha=.3)
  # print(p)
}

summary(distance)
summary(dur_jags)
summary(dur_rbest)



