library(RBesT)

dat <- crohn
crohn_sigma <- 88
dat$y.se <- crohn_sigma/sqrt(dat$n)

map_mcmc <- gMAP(cbind(y, y.se) ~ 1 | study, 
                 weights=n,data=dat,
                 family=gaussian,
                 beta.prior=cbind(0, crohn_sigma),
                 tau.dist="HalfNormal",tau.prior=cbind(0,crohn_sigma/2))

#autofix
map_a <- automixfit(map_mcmc)
#option 2: specify the number of mixture components
map_2 <- mixfit(map_mcmc, Nc = 2)

map_hat = map_2
#robustify
map_robust <- robustify(map_hat, weight=0.2, mean=-50)

#update with current control data
y.c <- -76.01
y.c.se <- 21.93
n.c <- 20

post_c <- postmix(map_robust, m=y.c, se=y.c)

#posterior sample from the posterior distribution
tt = rmix(post_c,100000)

