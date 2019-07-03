dat <- crohn
crohn_sigma <- 88
dat$y.se <- crohn_sigma/sqrt(dat$n)

library(RBesT)
set.seed(1234)
map_mcmc <- gMAP(cbind(y, y.se) ~ 1 | study, 
                 weights=n,data=dat,
                 family=gaussian,
                 beta.prior=cbind(0, crohn_sigma),
                 tau.dist="HalfNormal",tau.prior=cbind(0,crohn_sigma/2))
print(map_mcmc)

map <- automixfit(map_mcmc)

map_robust <- robustify(map, weight=0.2, mean=-50)

poc <- decision2S(pc=c(0.95,0.5), qc=c(0,-50), lower.tail=TRUE)
print(poc)

## set up prior for active group
weak_prior <- mixnorm(c(1,-50,1), sigma=crohn_sigma, param = 'mn')
n_act <- 40
n_pbo <- 20

## four designs
## "b" means a balanced design, 1:1
## "ub" means 40 in active and 20 in placebo
design_noprior_b  <- oc2S(weak_prior, weak_prior, n_act, n_act, poc,
                          sigma1=crohn_sigma, sigma2=crohn_sigma)
design_noprior_ub <- oc2S(weak_prior, weak_prior, n_act, n_pbo, poc,
                          sigma1=crohn_sigma, sigma2=crohn_sigma)
design_nonrob_ub  <- oc2S(weak_prior, map, n_act, n_pbo, poc,
                          sigma1=crohn_sigma, sigma2=crohn_sigma)
design_rob_ub     <- oc2S(weak_prior, map_robust, n_act, n_pbo, poc,
                          sigma1=crohn_sigma, sigma2=crohn_sigma)