rm(list=ls())

library(RBesT)
library(tidyverse)
source("01-NormalSimSpec0.R") #effect size 0
# source("01-NormalSimSpec1.R") #effect size -20

mu_c_tab = c(-60, -55, -50, -45, -40)

for(w in 1:length(w_vec)){
  w_v <- w_vec[w]
  rdname <- paste("es_",effsize,"_wv_",w_v,".RData",sep="")
  load(rdname)
  
  for(m in 1:length(muvec_c)){
    results = results_lst[[m]]
    
    tt1 <- tibble(TrueCtrlMean = muvec_c[m], Prop = colMeans(results[,1:4]), Method = c("MAP", "rMAP", "Alt1", "Alt2"))
    tt2 <- as_tibble(results[,5:8]) %>% mutate(TrueCtrlMean = muvec_c[m])
    tt3 <- tibble(results[,9], TrueCtrlMean = muvec_c[m])
    
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
  
  names(median_est) <- c("MAP", "rMAP","Alt1", "Alt2", "TrueCtrlMean")
  names(w_post) <- c("wPost", "TrueCtrlMean")
  
  w_summary <-
    w_post %>% group_by(TrueCtrlMean) %>% 
    summarize(mean = mean(wPost), stddev = sd(wPost), min = min(wPost), 
              q1 = quantile(wPost, 0.25), median = median(wPost),  q3 = quantile(wPost, 0.75), max = max(wPost))
  
  est_summary <-
    median_est %>% group_by(TrueCtrlMean) %>% 
    summarize_all(list(mean=mean, sd=sd, median=median)) %>%
    mutate(MAP_MSE = (MAP_mean - TrueCtrlMean)^2 + MAP_sd^2, 
           rMAP_MSE = (rMAP_mean - TrueCtrlMean)^2 + rMAP_sd^2,
           Alt1_MSE = (Alt1_mean - TrueCtrlMean)^2 + Alt1_sd^2,
           Alt2_MSE = (Alt2_mean - TrueCtrlMean)^2 + Alt2_sd^2
           )
  
  
  #assemble results
  for(t in 1:length(mu_c_tab)){
    part1 <- decision %>% filter(Method %in% c("MAP","rMAP", "Alt1", "Alt2") & TrueCtrlMean == mu_c_tab[t]) %>% select(Prop)
    part1 <- round(part1, 3)
    
    tt1 <- median_est %>% filter(TrueCtrlMean == mu_c_tab[t]) %>% group_by(TrueCtrlMean) %>% 
      summarize_all(list(mean=mean, sd=sd)) %>%
      mutate(MAP_MSE = (MAP_mean - TrueCtrlMean)^2 + MAP_sd^2, 
             rMAP_MSE = (rMAP_mean - TrueCtrlMean)^2 + rMAP_sd^2,
             Alt1_MSE = (Alt1_mean - TrueCtrlMean)^2 + Alt1_sd^2,
             Alt2_MSE = (Alt2_mean - TrueCtrlMean)^2 + Alt2_sd^2
      )
      
    part2 <- cbind(t(tt1[,c("MAP_mean","rMAP_mean", "Alt1_mean", "Alt2_mean")]),
                   t(tt1[,c("MAP_MSE", "rMAP_MSE", "Alt1_MSE", "Alt2_MSE")])
    )
    part2 <- round(part2,1)
    tt2 <- cbind(part1, part2)
    
    if(t==1)  out <- tt2 else out <- cbind(out, tt2)
  }
  
  ttt <- data.frame(w_v, Method = c("MAP", "rMAP","Alt1", "Alt2"), out)
  rownames(ttt) <- NULL
  
  if(w==1) outdt <- ttt else outdt <- rbind(outdt, ttt)
  
}


# write.csv(outdt,"out.csv")

