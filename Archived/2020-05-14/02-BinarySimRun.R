rm(list=ls())

library(R2jags)
library(RBesT)
library(ggplot2)
library(doParallel)
source("00-BhattacharyyaDistance.R")
#source("02-BinarySimSpec0.R") #effect size 0
source("02-BinarySimSpec1.R") #effect size 0.2

nsim = 500
ncores = min(parallel::detectCores(), 40)
cl = makeCluster(ncores)
registerDoParallel(cl)

#derive the rMAP as it doesn't rely on current data
options(RBesT.MC.control=list(adapt_delta=0.999))
map_c_mcmc = gMAP(cbind(r, n - r) ~ 1 | study, 
                  #weight = n,
                  data=dt,
                  family=binomial,
                  beta.prior=cbind(0, se_mu_c),
                  tau.dist="HalfNormal",tau.prior=cbind(0,HNscale_c),
                  chains = n.chains)
#approximate the MAP
map_c_hat = automixfit(map_c_mcmc)

#non-informative beta prior for current treatment arm
prior_t = mixbeta(c(1,a_t,b_t))

#################################################################
#################################################################
#################################################################

t0 <- Sys.time()

for(w in 1:length(w_vec)){
  w_v = w_vec[w]
  
  #robustification, note RBesT uses the mean/sample size-1 parametrization, see help of "robustify" for details
  rmap_c = robustify(map_c_hat, weight = w_v, mean = a_c/(a_c+b_c), n = a_c+b_c- 1)
  
  results_lst <- NULL
  
  for(m in 1:length(muvec_c)){
    print(m)
    ##################################
    #parallel computing at this level#
    ##################################
    results = 
      foreach(icount(nsim), .combine = rbind,.packages = c("RBesT","R2jags")) %dopar% {
        #simulate current control data
        y_c = rbinom(1, n_c, muvec_c[m])
        #simulate current treatment data
        y_t = rbinom(1, n_t, muvec_c[m] + effsize)
  
        ###############
        #Alternative 1#
        ###############
        jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, 
                                 y_c = y_c, n_c = n_c,
                                 y_h = dt$r, n_h = dt$n,
                                 y_t = y_t, n_t = n_t, prec_t = se_mu_t^-2, 
                                 w_v = w_v, robust_sd = robust_sd)
        )
        jags_obj = jags(model.file = "Model_BinaryHNMix.bugs",
                        data = jags_data,
                        parameters.to.save = c("p_c","p_t"),
                        n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                        progress.bar = "none"
        )
        jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
        tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
        post_c_alt1 = tt$p_c
        post_t_alt1 = tt$p_t
        
        #posterior difference
        post_diff_alt1 = post_t_alt1 - post_c_alt1
        decision_alt1 = (mean(post_diff_alt1 > Qcut) > Pcut)      
        
        
        ###############
        #Alternative 2#
        ###############
        #original MAP model
        jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, 
                                 y_c = y_c, n_c = n_c,
                                 y_h = dt$r, n_h = dt$n,
                                 y_t = y_t, n_t = n_t, prec_t = se_mu_t^-2, 
                                 w_v = 0, robust_sd = robust_sd)
        )
        jags_obj = jags(model.file = "Model_BinaryHNMix.bugs",
                        data = jags_data,
                        parameters.to.save = c("p_c","p_t"),
                        n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                        progress.bar = "none"
        )
        
        jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
        tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
        post_c_alt2_b = tt$p_c
        post_t_alt2 = tt$p_t
        
        #robust model with the vague prior
        jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, 
                                 y_c = y_c, n_c = n_c,
                                 y_h = dt$r, n_h = dt$n,
                                 y_t = y_t, n_t = n_t, prec_t = se_mu_t^-2, 
                                 w_v = 1, robust_sd = robust_sd)
        )
        jags_obj = jags(model.file = "Model_BinaryHNMix.bugs",
                        data = jags_data,
                        parameters.to.save = c("p_c","p_t"),
                        n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                        progress.bar = "none"
        )

        jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
        post_c_alt2_nb = data.frame(jags_auto$BUGSoutput$sims.matrix)$p_c
        
        #model averaging
        nb <- rbinom(length(post_diff_alt1),1,w_v)
        post_c_alt2 <- (1 - nb) * post_c_alt2_b + nb * post_c_alt2_nb
        
        #posterior difference
        post_diff_alt2 = post_t_alt2 - post_c_alt2
        decision_alt2 = (mean(post_diff_alt2 > Qcut) > Pcut)
        
        
        ###########################
        #Alternative 3: full Bayes#
        ###########################
        # jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, 
        #                          y_c = y_c, n_c = n_c,
        #                          y_h = dt$r, n_h = dt$n,
        #                          y_t = y_t, n_t = n_t, prec_t = se_mu_t^-2, 
        #                          w_v = w_v, robust_sd = robust_sd,
        #                          beta_a = w_v, beta_b = 1 - w_v)
        # )
        # jags_obj = jags(model.file = "Model_BinaryHNMix_FB.bugs",
        #                 data = jags_data,
        #                 parameters.to.save = c("p_c","p_t"),
        #                 n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
        #                 progress.bar = "none"
        # )
        # jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
        # tt = data.frame(jags_auto$BUGSoutput$sims.matrix)
        # post_c_alt3 = tt$p_c
        # post_t_alt3 = tt$p_t
        # 
        # #posterior difference
        # post_diff_alt3 = post_t_alt3 - post_c_alt3
        # decision_alt3 = (mean(post_diff_alt3 > Qcut) > Pcut)      
        
        post_c_alt3 = 1
        decision_alt3 = 0
        
        
        ##############
        #Genuine rMAP#
        ##############
        #get the posterior mixtures
        postmix_c = postmix(rmap_c,r = y_c, n = n_c)
        postmix_t = postmix(prior_t, r = y_t, n = n_t)
        
        post_c_rmap = rmix(mix = postmix_c, n = length(post_diff_alt1))
        post_t_rmap = rmix(mix = postmix_t, n = length(post_diff_alt1))
        
        post_diff_rmap = post_t_rmap - post_c_rmap
        decision_rmap = (mean(post_diff_rmap > Qcut) > Pcut)
        
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
