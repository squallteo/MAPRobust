rm(list=ls())

library(R2jags)
library(RBesT)
library(ggplot2)
library(doParallel)
source("00-BhattacharyyaDistance.R")
source("02-NormalSimSpec1.R")

nsim = 10
ncores = parallel::detectCores()
cl = makeCluster(ncores)
registerDoParallel(cl)

#initialization


for(m in 1:length(muvec_c)){

  #parallel computing at this level
  oc_rbest = foreach(icount(nsim), .combine = rbind,.packages = "RBesT") %dopar% {
    #simulate current control data
    sample_c = rnorm(n_c, muvec_c[m], sigma)
    y_c = mean(sample_c)
    se_yc = sd(sample_c)/sqrt(n_c)
    
    #simulate current treatment data
    sample_t = rnorm(n_t, muvec_c[m] + effsize, sigma)
    y_t = mean(sample_t)
    se_yt = sd(sample_t)/sqrt(n_t)
    
    options(RBesT.MC.control=list(adapt_delta=0.999))
    map_c_mcmc = gMAP(cbind(y, se_yh) ~ 1 | study, 
                      #weight = n,
                      data=dt,
                      family=gaussian,
                      beta.prior=cbind(0, se_mu_c),
                      tau.dist="HalfNormal",tau.prior=cbind(0,HNscale_c),
                      chains = n.chains)
    #approximate the MAP
    map_c_hat = mixfit(map_c_mcmc, Nc = 3)
    #robustify command needs a reference scale, although it won't be used
    sigma(map_c_hat) = sigma
    #robustification
    rmap_c = robustify(map_c_hat, weight = w_r, mean = 0, sigma = robust_sd)
    
    #calculate OCs
    #non-informative priors for current treatment arm
    prior_t = mixnorm(c(1,0,se_mu_t), sigma=sigma)
    design_rbest =   oc2S(prior1 = prior_t, prior2 = rmap_c, n1 = n_t, n2 = n_c, 
                          decision = success_rule, sigma1 = sigma, sigma2 = sigma)
    c(design_rbest(muvec_c[m],muvec_c[m]), design_rbest(muvec_c[m] + effsize,muvec_c[m]))
  }
  
  power_jags = foreach(icount(100), .combine = c,.packages = "R2jags") %dopar% {
    #simulate current control data
    sample_c = rnorm(n_c, muvec_c[m], sigma)
    y_c = mean(sample_c)
    se_yc = sd(sample_c)/sqrt(n_c)
    
    #simulate current treatment data
    sample_t = rnorm(n_t, muvec_c[m] + effsize, sigma)
    y_t = mean(sample_t)
    se_yt = sd(sample_t)/sqrt(n_t)
    
    #original MAP model
    jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, yh = yh, prec_yh = se_yh^-2,
                             prec_mu_t = se_mu_t^-2,                        y_t = y_t, prec_yt = se_yt^-2)
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
                             prec_mu_t = se_mu_t^-2,                        y_t = y_t, prec_yt = se_yt^-2)
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
    
    #posterior difference
    post_diff_jags = post_t_jags - post_c_jags
    
    1*(mean(post_diff_jags <= Qcut) > Pcut)
  }
  
  
  


  
  
}






