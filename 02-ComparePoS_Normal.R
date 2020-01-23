rm(list=ls())

library(R2jags)
library(RBesT)
library(ggplot2)
library(doParallel)
source("00-BhattacharyyaDistance.R")
source("02-NormalSimSpec1.R")

nsim = 1000
ncores = min(parallel::detectCores(), 40)
cl = makeCluster(ncores)
registerDoParallel(cl)

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
#robustification
rmap_c = robustify(map_c_hat, weight = w_v, mean = 0, sigma = robust_sd)

#non-informative priors for current treatment arm
prior_t = mixnorm(c(1,0,se_mu_t))


for(m in 1:length(muvec_c)){
  print(m)
  ##################################
  #parallel computing at this level#
  ##################################
  results = 
    foreach(icount(nsim), .combine = rbind,.packages = c("RBesT","R2jags")) %dopar% {
      #simulate current control data
      sample_c = rnorm(n_c, muvec_c[m], sigma)
      y_c = mean(sample_c)
      se_yc = sd(sample_c)/sqrt(n_c)
      
      #simulate current treatment data
      sample_t = rnorm(n_t, muvec_c[m] + effsize, sigma)
      y_t = mean(sample_t)
      se_yt = sd(sample_t)/sqrt(n_t)
      
      ###################
      #Mixture Prior MAC#
      ###################
      jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, 
                               yh = yh, prec_yh = se_yh^-2,
                               prec_mu_t = se_mu_t^-2,                        
                               y_t = y_t, prec_yt = se_yt^-2,
                               w_v = w_v, robust_sd = robust_sd)
      )
      jags_obj = jags(model.file = "Model_NormalHNMix.bugs",
                      data = jags_data,
                      parameters.to.save = c("mu_c","mu_t"),
                      n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                      progress.bar = "none"
      )
      jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
      tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
      post_c_mac = tt$mu_c
      post_t_mac = tt$mu_t
      
      #posterior difference
      post_diff_mac = post_t_mac - post_c_mac
      decision_mac = (mean(post_diff_mac <= Qcut) > Pcut)      

      #############################
      #mixture prior full Bayesian#
      #############################
      jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, 
                               yh = yh, prec_yh = se_yh^-2,
                               prec_mu_t = se_mu_t^-2,                        
                               y_t = y_t, prec_yt = se_yt^-2, robust_sd = robust_sd,
                               beta_a = w_v, beta_b = 1 - w_v)
      )
      jags_obj = jags(model.file = "Model_NormalHNMix_FB.bugs",
                      data = jags_data,
                      parameters.to.save = c("mu_c","mu_t"),
                      n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                      progress.bar = "none"
      )
      jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
      tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
      post_c_fb = tt$mu_c
      post_t_fb = tt$mu_t
      
      #posterior difference
      post_diff_fb = post_t_fb - post_c_fb
      decision_fb = (mean(post_diff_fb <= Qcut) > Pcut)      
      
      
      ######################
      #Robust mixture prior#
      ######################
      #original MAP model
      jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, 
                               yh = yh, prec_yh = se_yh^-2,
                               prec_mu_t = se_mu_t^-2,                        
                               y_t = y_t, prec_yt = se_yt^-2,
                               w_v = 0, robust_sd = robust_sd)
      )
      jags_obj = jags(model.file = "Model_NormalHNMix.bugs",
                      data = jags_data,
                      parameters.to.save = c("mu_c","mu_t"),
                      n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                      progress.bar = "none"
      )
      
      jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
      tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
      post_c_rmix_b = tt$mu_c
      post_t_rmix = tt$mu_t
      
      #robust model with a fixed tau
      jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, 
                               yh = yh, prec_yh = se_yh^-2,
                               prec_mu_t = se_mu_t^-2,                        
                               y_t = y_t, prec_yt = se_yt^-2,
                               w_v = 1, robust_sd = robust_sd)
      )
      jags_obj = jags(model.file = "Model_NormalHNMix.bugs",
                      data = jags_data,
                      parameters.to.save = c("mu_c","mu_t"),
                      n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                      progress.bar = "none"
      )
      
      jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
      post_c_rmix_nb = data.frame(jags_auto$BUGSoutput$sims.matrix)$mu_c
      
      #model averaging
      nb <- rbinom(length(post_diff_fb),1,w_v)
      post_c_rmix <- (1 - nb) * post_c_rmix_b + nb * post_c_rmix_nb
      # post_c_rmix = (1 - w_v) * post_c_rmix_b + w_v * post_c_rmix_nb
      
      #posterior difference
      post_diff_rmix = post_t_rmix - post_c_rmix
      decision_rmix = (mean(post_diff_rmix <= Qcut) > Pcut)
      
      ######
      #rMAP#
      ######
      #get the posterior mixtures
      postmix_c = postmix(rmap_c,m = y_c, se = se_yc)
      postmix_t = postmix(prior_t, m = y_t, se = se_yt)
      
      post_c_rmap = rmix(mix = postmix_c, n = length(post_diff_fb))
      post_t_rmap = rmix(mix = postmix_t, n = length(post_diff_fb))
      
      post_diff_rmap = post_t_rmap - post_c_rmap
      decision_rmap = (mean(post_diff_rmap <= Qcut) > Pcut)
      
      
      c(decision_mac, decision_rmix, decision_fb, decision_rmap)
    }
  
  ttt = data.frame(Prop = colMeans(results), Method = c("MAC","rMix", "FullBayes", "rMAP"), TrueCtrlMean = muvec_c[m])
  if(m==1) outdt = ttt else outdt = rbind(outdt,ttt)
  
}

stopCluster(cl)


if(effsize != 0){
  ggplot(data=outdt,aes(x=TrueCtrlMean,y=Prop,color=Method)) + geom_line(size = 1.5) + 
    scale_x_continuous(breaks = muvec_c, name = "True Control Mean") + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), name = "Probability of Success") + 
    geom_vline(xintercept = -49.9, linetype=2, size=1) + 
    theme(axis.title = element_text(face="bold",size=15),
          axis.text = element_text(size=12),
          legend.title=element_text(size=15,face="bold"),
          legend.text=element_text(size=12)
    )
}

if(effsize == 0){
  ggplot(data=outdt,aes(x=TrueCtrlMean,y=Prop,color=Method)) + geom_line(size = 1.5) + 
    scale_x_continuous(breaks = muvec_c, name = "True Control Mean") + 
    scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.01), name = "Probability of Success") + 
    geom_vline(xintercept = -49.9, linetype=2, size=1) + 
    geom_hline(yintercept = 0.025, linetype=2, size=1) + 
    theme(axis.title = element_text(face="bold",size=15),
          axis.text = element_text(size=12),
          legend.title=element_text(size=15,face="bold"),
          legend.text=element_text(size=12)
    )
}
