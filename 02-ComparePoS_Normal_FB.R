rm(list=ls())

library(R2jags)
library(RBesT)
library(ggplot2)
library(doParallel)
source("00-BhattacharyyaDistance.R")
source("02-NormalSimSpec1.R")

nsim = 500
ncores = parallel::detectCores()
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
map_c_hat = mixfit(map_c_mcmc, Nc = 3)
#robustify command needs a reference scale, although it won't be used
sigma(map_c_hat) = sigma
#robustification
rmap_c = robustify(map_c_hat, weight = w_r, mean = 0, sigma = robust_sd)

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
      
      ######################
      #Robust mixture prior#
      ######################
      jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, yh = yh, prec_yh = se_yh^-2,
                               prec_mu_t = se_mu_t^-2,                        y_t = y_t, prec_yt = se_yt^-2,
                               w_r = w_r, robust_sd = robust_sd)
      )
      jags_obj = jags(model.file = "Model_NormalHNMix.bugs",
                      data = jags_data,
                      parameters.to.save = c("mu_c","mu_t"),
                      n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                      progress.bar = "none"
      )
      jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
      tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
      post_c_jags = tt$mu_c
      post_t_jags = tt$mu_t
      
      #posterior difference
      post_diff_jags = post_t_jags - post_c_jags
      
      decision_jags = (mean(post_diff_jags <= Qcut) > Pcut)
      
      #get the posterior mixtures
      postmix_c = postmix(rmap_c,m = y_c, se = se_yc)
      postmix_t = postmix(prior_t, m = y_t, se = se_yt)
      
      post_c_rbest = rmix(mix = postmix_c, n = length(post_diff_jags))
      post_t_rbest = rmix(mix = postmix_t, n = length(post_diff_jags))
      
      post_diff_rbest = post_t_rbest - post_c_rbest
      decision_rbest = (mean(post_diff_rbest <= Qcut) > Pcut)
      
      c(decision_jags,decision_rbest)
    }
  
  ttt = data.frame(Prop = colMeans(results), Method = c("JAGS", "RBesT"), TrueCtrlMean = muvec_c[m])
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
