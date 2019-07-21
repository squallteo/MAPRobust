rm(list=ls())

library(doParallel)
library(R2jags)
library(RBesT)
library(ggplot2)
source("00-BhattacharyyaDistance.R")
source("02-NormalSimSpec1.R")

ncores = parallel::detectCores()
cl = makeCluster(ncores)
registerDoParallel(cl)

m = 1
result = 
foreach(icount(100), .combine = rbind,.packages = c("RBesT","R2jags")) %dopar% {
  #simulate current control data
  sample_c = rnorm(n_c, muvec_c[m], sigma)
  y_c = mean(sample_c)
  se_yc = sd(sample_c)/sqrt(n_c)
  
  #simulate current treatment data
  sample_t = rnorm(n_t, muvec_c[m] + effsize, sigma)
  y_t = mean(sample_t)
  se_yt = sd(sample_t)/sqrt(n_t)
  
  
  jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, yh = yh, prec_yh = se_yh^-2, 
                           prec_mu_t = se_mu_t^-2, y_t = y_t, prec_yt = se_yt^-2)
  )
  jags_obj = jags(model.file = "Model_NormalHN.bugs",
                  data = jags_data,
                  parameters.to.save = c("mu_c","mu_t"),
                  n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                  progress.bar = "none"
  )
  
  jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
  tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
  post_t_jags = tt$mu_t
  
  #calculate OCs
  #non-informative priors for current treatment arm
  prior_t = mixnorm(c(1,0,se_yt), sigma=sigma, param = 'mn')
  post_t_mix = postmix(prior_t,m = y_t, se = se_yt)
  post_t_rbest = rmix(mix = post_t_mix,n = length(post_t_jags))
  
  c(var(post_t_rbest), var(post_t_jags), bhatt.coeff(post_t_rbest,post_t_jags))
}




