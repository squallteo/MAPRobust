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
muvec_c = seq(-60, -40, by=2) #current control mean vector

#control arm hyper-parameters
se_mu_c = 100
HNscale_c = 44

#robustification parameter: the weight associated with the robust part
w_r = 0.2
robust_sd = 200

#current treatment data/parameters
se_mu_t = 100
n_t = 40
effsize = -15

#decision rule to claim trial success
#pr(mu_t - mu_c < Qcut) > Pcut
Qcut = 0
Pcut = 0.975
success_rule = decision2S(pc = Pcut, qc = Qcut, lower.tail = T, link = "identity")

#MCMC control parameters
n.chains = 3