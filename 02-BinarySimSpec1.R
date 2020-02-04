#historical control data
dt = AS
sigma = 88

#meta analysis of historical data
# with(dt, meta::metaprop(event = r, n = n, method = "Inverse"))
#point est is 0.25 (0.202, 0.312). Decide to consider true current rate from 0.2 to 0.32

#current control data
n_c = 20 #current control sample size
muvec_c = seq(0.2, 0.32, by=0.01) #current control mean vector

#control arm hyper-parameters
se_mu_c = 2
HNscale_c = 0.5

#robustification parameter: the weight associated with the robust part
w_vec = c(0.2, 0.5, 0.8)
a_c = 1; b_c = 1 #vague beta prior for p(c), rMAP only
robust_sd = 10 #vague normal prior sd of logit p(c), alternatives only

#current treatment data/parameters
se_mu_t = 2 #sd of logit p(t), alternatives only
a_t = 0.5; b_t = 0.5 #beta prior for p(t), rMAP only
n_t = 40
effsize = 0.2

#decision rule to claim trial success
#pr(mu_t - mu_c < Qcut) > Pcut
Qcut = 0
Pcut = 0.95
success_rule = decision2S(pc = Pcut, qc = Qcut, lower.tail = F, link = "identity")

#MCMC control parameters
n.chains = 3