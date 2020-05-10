rm(list=ls())

library(R2jags)
library(RBesT)
library(ggplot2)
library(doParallel)
library(tidyverse)
#source("02-NormalSimSpec0.R") #effect size 0
source("02-NormalSimSpec1.R") #effect size -20
muvec_c = c(-60, -55, -50, -45, -40)

nsim = 5000
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
map_c_hat = mixfit(map_c_mcmc,Nc = 2)
#robustify command needs a reference scale, although it won't be used
sigma(map_c_hat) = sigma

#non-informative priors for current treatment arm
prior_t = mixnorm(c(1,0,se_mu_t))

#################################################################
#################################################################
#################################################################

outlist <- NULL

for(w in 1:length(w_vec)){
  w_v = w_vec[w]
  
  #robustification
  rmap_c = robustify(map_c_hat, weight = w_v, mean = 0, sigma = robust_sd)
  
  results_lst <- NULL
  
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
        
        ##############
        #Genuine rMAP#
        ##############
        #get the posterior mixtures
        postmix_c = postmix(rmap_c,m = y_c, se = se_yc)
        postmix_t = postmix(prior_t, m = y_t, se = se_yt)
        
        post_c_rmap = rmix(mix = postmix_c, n = 10000)
        post_t_rmap = rmix(mix = postmix_t, n = 10000)
        
        post_diff_rmap = post_t_rmap - post_c_rmap
        decision_rmap = (mean(post_diff_rmap <= Qcut) > Pcut)
        
        #assemble the results
        w_post <- tail(as.matrix(postmix_c)[1,],1)
        decisions <- c(decision_rmap)
        median_est <- c(median(post_c_rmap))
        
        c(decisions, median_est, w_post)
      }
    
    results_lst[[m]] <- results
    
  }
  
  for(m in 1:length(muvec_c)){
    results = results_lst[[m]]
    
    tt1 <- tibble(TrueCtrlMean = muvec_c[m], Prop = mean(results[,1]), Method = "rMAP")
    tt2 <- as_tibble(results[,2]) %>% mutate(TrueCtrlMean = muvec_c[m])
    tt3 <- tibble(results[,3], TrueCtrlMean = muvec_c[m])
    
    if(m==1){
      decision <- tt1
      median_est <- tt2
      w_post <- tt3
    }
    else{
      decision <- rbind(decision, tt1)
      median_est <- rbind(median_est, tt2)
      w_post <- rbind(w_post, tt3)
    }
  }
  
  names(median_est) <- c("rMAP", "TrueCtrlMean")
  names(w_post) <- c("wPost", "TrueCtrlMean")
  
  w_summary <-
    w_post %>% group_by(TrueCtrlMean) %>% 
    summarize(mean = mean(wPost), stddev = sd(wPost), min = min(wPost), 
              q1 = quantile(wPost, 0.25), median = median(wPost),  q3 = quantile(wPost, 0.75), max = max(wPost))
  
  est_summary <-
    median_est %>% group_by(TrueCtrlMean) %>% 
    summarize_all(list(mean=mean, sd=sd, median=median)) %>%
    mutate(rMAP_MSE = (mean - TrueCtrlMean)^2 + sd^2)
  
  outlist[[w]] = cbind(decision[,1:2], est_summary[,c("mean", "rMAP_MSE")])
  
}


stopCluster(cl)
