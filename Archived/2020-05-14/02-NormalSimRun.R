rm(list=ls())

library(R2jags)
library(RBesT)
library(ggplot2)
library(doParallel)
source("00-BhattacharyyaDistance.R")
#source("02-NormalSimSpec0.R") #effect size 0
source("02-NormalSimSpec1.R") #effect size -20

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
map_c_hat = automixfit(map_c_mcmc)
#robustify command needs a reference scale, although it won't be used
sigma(map_c_hat) = sigma

#non-informative priors for current treatment arm
prior_t = mixnorm(c(1,0,se_mu_t))

#################################################################
#################################################################
#################################################################

t0 <- Sys.time()

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
        
        ###############
        #Alternative 1#
        ###############
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
        post_c_alt1 = tt$mu_c
        post_t_alt1 = tt$mu_t
        
        #posterior difference
        post_diff_alt1 = post_t_alt1 - post_c_alt1
        decision_alt1 = (mean(post_diff_alt1 <= Qcut) > Pcut)      
        
        
        ###############
        #Alternative 2#
        ###############
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
        post_c_alt2_b = tt$mu_c
        post_t_alt2 = tt$mu_t
        
        #robust model with the vague prior
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
        post_c_alt2_nb = data.frame(jags_auto$BUGSoutput$sims.matrix)$mu_c
        
        #model averaging
        nb <- rbinom(length(post_diff_alt1),1,w_v)
        post_c_alt2 <- (1 - nb) * post_c_alt2_b + nb * post_c_alt2_nb
        
        #posterior difference
        post_diff_alt2 = post_t_alt2 - post_c_alt2
        decision_alt2 = (mean(post_diff_alt2 <= Qcut) > Pcut)
        
        
        ###########################
        #Alternative 3: full Bayes#
        ###########################
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
        post_c_alt3 = tt$mu_c
        post_t_alt3 = tt$mu_t
        
        #posterior difference
        post_diff_alt3 = post_t_alt3 - post_c_alt3
        decision_alt3 = (mean(post_diff_alt3 <= Qcut) > Pcut)      
        
        
        ##############
        #Genuine rMAP#
        ##############
        #get the posterior mixtures
        postmix_c = postmix(rmap_c,m = y_c, se = se_yc)
        postmix_t = postmix(prior_t, m = y_t, se = se_yt)
        
        post_c_rmap = rmix(mix = postmix_c, n = length(post_diff_alt1))
        post_t_rmap = rmix(mix = postmix_t, n = length(post_diff_alt1))
        
        post_diff_rmap = post_t_rmap - post_c_rmap
        decision_rmap = (mean(post_diff_rmap <= Qcut) > Pcut)
        
        #assemble the results
        w_post <- tail(as.matrix(postmix_c)[1,],1)
        decisions <- c(decision_rmap, decision_alt1, decision_alt2, decision_alt3)
        median_est <- c(median(post_c_rmap), median(post_c_alt1), median(post_c_alt2), median(post_c_alt3))
        
        c(decisions, median_est, w_post)
      }
    
    results_lst[[m]] <- results
    
  }
  
  rdname <- paste("es_",effsize,"_wv_",w_v,".RData",sep="")
  save(results_lst,file = rdname)
}

stopCluster(cl)

t1 <- Sys.time()

t1 - t0
