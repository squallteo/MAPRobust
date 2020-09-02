#historical control data
dt = AS

#meta analysis of historical data
# with(dt, meta::metaprop(event = r, n = n, method = "Inverse"))
#point est is 0.25 (0.202, 0.312). Decide to consider true current rate from 0.2 to 0.32

#current control data
n_c = 20 #current control sample size
# muvec_c = seq(0.2, 0.32, by=0.02) #current control mean vector
muvec_c = c(0.20, 0.22, 0.25, 0.28, 0.32)

#control arm hyper-parameters
se_mu_c = 2
HNscale_c = 0.5

#robustification parameter: the weight associated with the robust part
w_vec = c(0.2, 0.5, 0.8)
a_c = 0.5; b_c = 0.5 #vague beta prior for p(c), rMAP only
robust_sd = 2 #vague normal prior sd of logit p(c), alternatives only

#current treatment data/parameters
se_mu_t = 2 #sd of logit p(t), alternatives only
a_t = 0.5; b_t = 0.5 #beta prior for p(t), rMAP only
n_t = 40
effsize = 0.3

#decision rule to claim trial success
#pr(mu_t - mu_c < Qcut) > Pcut
Qcut = 0
Pcut = 0.95
success_rule = decision2S(pc = Pcut, qc = Qcut, lower.tail = F, link = "identity")

#MCMC control parameters
n.chains = 3

#simulation to compute PoS for frequentist no borrowing
#comment out in simulation
# nsim = 10000
# 
# mean_c <- 0.25
# mean_t <- mean_c + effsize
# 
# rej <- rep(NA, nsim)
# 
# for(s in 1:nsim){
#   s_c <- rbinom(1, n_c, mean_c)
#   s_t <- rbinom(1, n_t, mean_t)
#   
#   tt <- prop.test(c(s_t, s_c), c(n_t, n_c), alternative = "greater", conf.level = 0.95)
#   rej[s] <- (tt$p.value < 0.05)
# }
# 
# 
# mean(rej)