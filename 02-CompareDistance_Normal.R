rm(list=ls())

library(doParallel)
library(R2jags)
library(RBesT)
library(ggplot2)
source("00-BhattacharyyaDistance.R")
source("02-NormalSimSpec1.R")

nsim = 500
ncores = detectCores()
cl = makeCluster(ncores)
registerDoParallel(cl)



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

for(m in 1:length(muvec_c)){

  results = 
  foreach(icount(nsim), .combine = rbind, .packages = c("R2jags","RBesT")) %dopar% {
    #simulate current control data
    sample_c = rnorm(n_c, muvec_c[m], sigma)
    y_c = mean(sample_c)
    se_yc = sd(sample_c)/sqrt(n_c)
    
    #simulate current treatment data
    sample_t = rnorm(n_t, muvec_c[m] + effsize, sigma)
    y_t = mean(sample_t)
    se_yt = sd(sample_t)/sqrt(n_t)
    
    jags_data = with(dt,list(prec_mu_c = se_mu_c^-2, HNscale_c = HNscale_c, y_c = y_c, prec_yc = se_yc^-2, yh = yh, prec_yh = se_yh^-2, 
                             prec_mu_t = se_mu_t^-2,                        y_t = y_t, prec_yt = se_yt^-2)
    )
    jags_obj = jags(model.file = "Model_NormalHN.bugs",
                    data = jags_data,
                    parameters.to.save = c("mu_c"),
                    n.chains = n.chains, n.burnin = 2000, n.iter = 10000,
                    progress.bar = "none"
    )
    
    jags_auto = autojags(jags_obj, Rhat = 1.1, n.thin = 4, n.iter = 40000, n.update = 10, progress.bar = "none")
    post_c_jags_b = data.frame(jags_auto$BUGSoutput$sims.matrix)$mu_c
    
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
    
    ##########################################################
    ##########################################################
    ##########################################################
    #update with current control data
    post_c_mix = postmix(rmap_c,m = y_c, se = se_yc)
    post_c_rbest = rmix(mix = post_c_mix,n = length(post_c_jags))
    
    #Bhattacharyya distance
    dist = bhatt.coeff(post_c_jags,post_c_rbest)
    #Point Estimates: median and mean for JAGS and RBesT respectively
    PtEst = c(
      median(post_c_jags), mean(post_c_jags),
      median(post_c_rbest), mean(post_c_rbest)
    )
      
    #quantity to be returned in each iteration
    c(dist,PtEst)
    
  }
  
  #MSE based on MEDIAN point estimators
  MyMSE = function(sample,truth){
    p1 = var(sample)
    p2 = mean(sample - truth)^2
    data.frame(Val = c(p1 + p2, p2/(p1+p2)), Type = c("MSE", "BiasProp"))
  }
  
  tt1 = data.frame(MyMSE(results[,2],muvec_c[m]), Method = "JAGS", TrueCtrlMean = muvec_c[m])
  tt2 = data.frame(MyMSE(results[,4],muvec_c[m]), Method = "RBesT", TrueCtrlMean = muvec_c[m])
  tt3 = data.frame(Val = mean(results[,1]), Type = "Distance", Method = "", TrueCtrlMean = muvec_c[m])
  ttt = rbind(tt1,tt2,tt3)
  
  if(m==1) outdt = ttt else outdt = rbind(outdt,ttt)
  
}

stopCluster(cl)

plotdt = subset(outdt,Type=="MSE")
ggplot(data=plotdt,aes(x=TrueCtrlMean,y=Val,color=Method)) + geom_line(size = 1.5) + 
  scale_x_continuous(breaks = muvec_c, name = "True Control Mean") + 
  ylab("Mean Square Error") + geom_vline(xintercept = -49.9, linetype=2, size=1) + 
  theme(axis.title = element_text(face="bold",size=15),
        axis.text = element_text(size=12),
        legend.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=12)
        )

plotdt = subset(outdt,Type=="Distance")
ggplot(data=plotdt,aes(x=TrueCtrlMean,y=Val)) + geom_line(size = 1.5) + 
  scale_x_continuous(breaks = muvec_c, name = "True Control Mean") + 
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), name = "Bhattacharyya Coefficient") + 
  geom_hline(yintercept = 0.9, linetype=2, size=1) + 
  theme(axis.title = element_text(face="bold",size=15),
        axis.text = element_text(size=12),
        legend.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=12)
  )

plotdt = subset(outdt,Type=="BiasProp")
ggplot(data=plotdt,aes(x=TrueCtrlMean,y=1-Val,color=Method)) + geom_line(size = 1.5) + 
  scale_x_continuous(breaks = muvec_c, name = "True Control Mean") + 
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), name = "Proportion of Variance in MSE") + 
  geom_vline(xintercept = -49.9, linetype=2, size=1) + 
  theme(axis.title = element_text(face="bold",size=15),
        axis.text = element_text(size=12),
        legend.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=12)
  )
                                     