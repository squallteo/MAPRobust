library(RBesT)
library(ggplot2)

histdt = data.frame(yh = c(51,49,36,47), se_yh = sqrt(c(10,6,5,20)), study = 1:4)

map_mcmc <- gMAP(cbind(yh, se_yh) ~ 1 | study, 
                 data=histdt,
                 family=gaussian,
                 beta.prior=cbind(0, 100),
                 tau.dist="HalfNormal",tau.prior=cbind(0,44))

#approximation of MAP
# map <- automixfit(map_mcmc)
# print(map)
# plot(map)$mix

map_hat <- mixfit(map_mcmc, Nc = 3)

#update approximated map with current trial data
y_c = 29
se_yc = 22
n_c = 20

post_c = postmix(map_hat, m=y_c, se=se_yc)


# p1 <- pmixdiff(post_act, post_pbo, 0); print(p1)
# p2 <- pmixdiff(post_act, post_pbo, -50); print(p2)
# print(p1>0.95 & p2>0.5)
# poc(post_act, post_pbo)
# 
# #robustification
# ## add a 20% non-informative mixture component
map_robust <- robustify(map_hat, weight=0.2, mean=-50)
# 
# #decision rules
# ## dual decision criteria
# ## pay attention to "lower.tail" argument and the order of active and pbo
# poc <- decision2S(pc=c(0.95,0.5), qc=c(0,-50), lower.tail=TRUE)
# print(poc)
# 
# #design options
# ## set up prior for active group
# weak_prior <- mixnorm(c(1,-50,1), sigma=crohn_sigma, param = 'mn')
# n_act <- 40
# n_pbo <- 20
# 
# ## four designs
# ## "b" means a balanced design, 1:1
# ## "ub" means 40 in active and 20 in placebo
# design_noprior_b  <- oc2S(weak_prior, weak_prior, n_act, n_act, poc,
#                           sigma1=crohn_sigma, sigma2=crohn_sigma)
# design_noprior_ub <- oc2S(weak_prior, weak_prior, n_act, n_pbo, poc,
#                           sigma1=crohn_sigma, sigma2=crohn_sigma)
# design_nonrob_ub  <- oc2S(weak_prior, map, n_act, n_pbo, poc,
#                           sigma1=crohn_sigma, sigma2=crohn_sigma)
# design_rob_ub     <- oc2S(weak_prior, map_robust, n_act, n_pbo, poc,
#                           sigma1=crohn_sigma, sigma2=crohn_sigma)
# 
# #type I error
# cfb_truth <- seq(-120, -40, by=1)
# 
# typeI1 <- design_noprior_b(cfb_truth, cfb_truth)
# typeI2 <- design_noprior_ub(cfb_truth, cfb_truth)
# typeI3 <- design_nonrob_ub(cfb_truth, cfb_truth)
# typeI4 <- design_rob_ub(cfb_truth, cfb_truth)
# 
# ocI <- rbind(data.frame(cfb_truth=cfb_truth, typeI=typeI1,
#                         design="40:40 with non-informative priors"),
#              data.frame(cfb_truth=cfb_truth, typeI=typeI2,
#                         design="40:20 with non-informative priors"),
#              data.frame(cfb_truth=cfb_truth, typeI=typeI3,
#                         design="40:20 with non-robust prior for placebo"),
#              data.frame(cfb_truth=cfb_truth, typeI=typeI4,
#                         design="40:20 with robust prior for placebo")
# )
# 
# qplot(cfb_truth, typeI, data=ocI, colour=design, geom="line", main="Type I Error") +
#   xlab(expression(paste('True value of change from baseline ', mu[act] == mu[pbo]))) +
#   ylab('Type I error') +
#   coord_cartesian(ylim=c(0,0.2)) +
#   theme(legend.justification=c(1,1),legend.position=c(0.95,0.85))
# 
# #power
# delta <- seq(-80,0,by=1)
# m <- summary(map)["mean"]
# cfb_truth1 <- m + delta   # active for 1
# cfb_truth2 <- m + 0*delta # pbo for 2
# 
# power1 <- design_noprior_b(cfb_truth1, cfb_truth2)
# power2 <- design_noprior_ub(cfb_truth1, cfb_truth2)
# power3 <- design_nonrob_ub(cfb_truth1, cfb_truth2)
# power4 <- design_rob_ub(cfb_truth1, cfb_truth2)
# 
# ocP <- rbind(data.frame(cfb_truth1=cfb_truth1, cfb_truth2=cfb_truth2,
#                         delta=delta, power=power1,
#                         design="40:40 with non-informative priors"),
#              data.frame(cfb_truth1=cfb_truth1, cfb_truth2=cfb_truth2,
#                         delta=delta, power=power2,
#                         design="40:20 with non-informative priors"),
#              data.frame(cfb_truth1=cfb_truth1, cfb_truth2=cfb_truth2,
#                         delta=delta, power=power3,
#                         design="40:20 with non-robust prior for placebo"),
#              data.frame(cfb_truth1=cfb_truth1, cfb_truth2=cfb_truth2,
#                         delta=delta, power=power4,
#                         design="40:20 with robust prior for placebo")
# )
# 
# qplot(delta, power, data=ocP, colour=design, geom="line", main="Power") +
#   xlab('True value of difference (act - pbo)')+ ylab('Power') +
#   scale_y_continuous(breaks=c(seq(0,1,0.2),0.9)) +
#   scale_x_continuous(breaks=c(seq(-80,0,20),-70)) +
#   geom_hline(yintercept=0.9,linetype=2) + 
#   geom_vline(xintercept=-70,linetype=2) +
#   theme(legend.justification=c(1,1),legend.position=c(0.95,0.85))
# 
# #final analysis with trial data
# ## one can either use summary data or individual data. See ?postmix.
# y.act <- -29.15
# y.act.se <- 16.69
# n.act <- 39
# 
