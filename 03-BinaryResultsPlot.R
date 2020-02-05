rm(list=ls())

library(ggplot2)
library(RBesT)
library(tidyverse)
# source("02-NormalSimSpec0.R") #effect size 0
source("02-NormalSimSpec1.R") #effect size -20

w_v <- 0.5
rdname <- paste("Results/Normal/es_",effsize,"_wv_",w_v,".RData",sep="")
load(rdname)

for(m in 1:length(muvec_c)){
  results = results_lst[[m]]
  
  tt1 <- tibble(TrueCtrlMean = muvec_c[m], Prop = colMeans(results[,1:4]), Method = c("rMAP","Alt1", "Alt2", "Alt3"))
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

names(w_post) <- c("wPost", "TrueCtrlMean")

w_post %>% 
  ggplot(aes(x=TrueCtrlMean, y=wPost, group=TrueCtrlMean)) + geom_boxplot(outlier.size = 1) +
  scale_x_continuous(breaks = muvec_c, name = "True Control Mean") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), name = expression("Posterior w"[V])) +
  geom_hline(yintercept = w_v, linetype=2, color = "red", size=1)
  
# 
# median_est %>% ggplot(aes(x=TrueCtrlMean, y=Alt3, group=TrueCtrlMean)) + geom_boxplot()
# 
# if(effsize != 0){
#   decision %>% filter(Method != "Alt3") %>%
#   ggplot(aes(x=TrueCtrlMean,y=Prop,color=Method)) + geom_line(size = 1.5) +
#     scale_x_continuous(breaks = muvec_c, name = "True Control Mean") +
#     scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), name = "Probability of Success") +
#     geom_vline(xintercept = -49.9, linetype=2, size=1) +
#     theme(axis.title = element_text(face="bold",size=15),
#           axis.text = element_text(size=12),
#           legend.title=element_text(size=15,face="bold"),
#           legend.text=element_text(size=12)
#     )
# }
# 
# if(effsize == 0){
#   decision %>% filter(Method != "Alt3") %>%
#   ggplot(aes(x=TrueCtrlMean,y=Prop,color=Method)) + geom_line(size = 1.5) +
#     scale_x_continuous(breaks = muvec_c, name = "True Control Mean") +
#     scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.01), name = "Probability of Success") +
#     geom_vline(xintercept = -49.9, linetype=2, size=1) +
#     geom_hline(yintercept = 0.025, linetype=2, size=1) +
#     theme(axis.title = element_text(face="bold",size=15),
#           axis.text = element_text(size=12),
#           legend.title=element_text(size=15,face="bold"),
#           legend.text=element_text(size=12)
#     )
# }
