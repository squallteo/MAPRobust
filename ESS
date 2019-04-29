rm(list=ls())

library(RBesT)
library(ggplot2)
set.seed(712)
n_ctrl = 20
histdt = data.frame(study=1:4, n=c(56,63,121,123), r=c(6,9,18,7))

map_mcmc = gMAP(cbind(r, n-r) ~ 1 | study,
                data=histdt,
                tau.dist="HalfNormal",
                tau.prior=0.5,
                beta.prior=2,
                family=binomial)

map_approx = automixfit(map_mcmc)

ess_morita = ess_moment = rep(NA,n_ctrl+1)

for(r in 0:n_ctrl){
  post_ctrl = postmix(map_approx,r = r, n = n_ctrl)
  ess_morita[r+1] = ess(post_ctrl,method = "morita")
  ess_moment[r+1] = ess(post_ctrl)
}

ess_Hobbs = c(-17,15,35,41,31,20,9,4,2,0,-1,-1,-3,-2,-3,-3,-5,-5,-8,-11,-19)

tt1 = data.frame(y = 0:n_ctrl, ESS = ess_Hobbs, Method = "Hobbs")
tt2 = data.frame(y = 0:n_ctrl, ESS = ess_morita, Method = "Morita")
tt3 = data.frame(y = 0:n_ctrl, ESS = ess_moment, Method = "Moment")
plotdt = rbind(tt1,tt2,tt3)

ggplot(plotdt,aes(x=y,y=ESS,color=Method,label=ESS)) + geom_path(size=1.5) + geom_point(size=6) + 
  geom_text(color="black",size=7) +
  scale_x_continuous(breaks=0:n_ctrl,name="Observed number of events in current control arm (y/20)") +
  geom_hline(yintercept = 0,size = 1) +
  theme(legend.text=element_text(size=20),
        legend.key.size=unit(2,"cm"),
        axis.title.x = element_text(face="bold",size=15),
        axis.text.x  = element_text(vjust=0.5, size=12,face="bold"),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y  = element_text(size=12,face="bold")
  )
