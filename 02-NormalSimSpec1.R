#historical control data
dt = crohn
sigma = 88
dt$se_yh = sigma/sqrt(dt$n)
yh = dt$y
nh = dt$n

#meta analysis of historical data
# with(dt, meta::metamean(n = n, mean = y, sd = rep(sigma,length(study))))
#point est is -50. Decide to consider true current mean from -60 to -40

#current control data
n_c = 20 #current control sample size
muvec_c = seq(-60, -40, by=5) #current control mean vector

#control arm hyper-parameters
se_mu_c = 100
HNscale_c = 44

#robustification parameter: the weight associated with the robust part
w_vec = c(0.2, 0.5, 0.8)
robust_mean = -50
robust_sd = 88

#current treatment data/parameters
se_mu_t = 100
n_t = 40
effsize = -60

#decision rule to claim trial success
#pr(mu_t - mu_c < Qcut) > Pcut
Qcut = 0
Pcut = 0.975
success_rule = decision2S(pc = Pcut, qc = Qcut, lower.tail = T, link = "identity")

#MCMC control parameters
n.chains = 3


#simulation to compute PoS for frequentist no borrowing
#comment out in simulation
# nsim = 10000
# 
# mean_c <- -50
# mean_t <- mean_c + effsize
# common_sd <- 88
# 
# rej <- rep(NA, nsim)
# 
# for(s in 1:nsim){
#   s_c <- rnorm(n_c, mean_c, common_sd)
#   s_t <- rnorm(n_t, mean_t, common_sd)
# 
#   tt <- t.test(s_t, s_c, alternative = "less", conf.level = 0.975)
#   rej[s] <- (tt$p.value < 0.025)
# }
# 
# 
# mean(rej)
