rm(list=ls())
library(RBesT)
library(tidyverse)
library(ggplot2)
library(pwr)

n_t <- 60
n_c <- 20

ctrl_range <- seq(0.05, 0.35, 0.01)
es <- 0.3

#frequentist power calculation
#20 current control
FPwr20 <- pwr.2p2n.test(h = ES.h(ctrl_range, ctrl_range + es), n1 = 20, n2 = n_t, sig.level = 0.05, alternative = "two.sided")
#40 current control
FPwr40 <- pwr.2p2n.test(h = ES.h(ctrl_range, ctrl_range + es), n1 = 40, n2 = n_t, sig.level = 0.05, alternative = "two.sided")

# tibble(Ctrl = ctrl_range, FPwr20 = FPwr20$power, FPwr40 = FPwr40$power) %>% filter(Ctrl %in% seq(0.05, 0.35, 0.05)) %>% View()

se_mu_c <- 2
se_mu_t <- 2

#historical data
histdt <- tibble(study = 1:4, n = c(221, 651, 20, 214), r = c(33, 98, 3, 36))
# with(histdt, meta::metaprop(event = r, n = n, method = "Inverse", prediction = T))
histdt_dummy <- tibble(study = 1, n = 1, r = 0)

#wrapper function to generate operating characteristics
oc_RBesT <- function(rule, histdt, label, ...){
  #MAP prior for controls
  options(RBesT.MC.control=list(adapt_delta=0.999))
  map_c_mcmc = gMAP(cbind(r, n - r) ~ 1 | study, 
                    data=histdt,
                    family=binomial,
                    beta.prior=cbind(0, se_mu_c),
                    tau.dist="HalfNormal",tau.prior=cbind(0,HNscale_c),
                    chains = 3)
  
  
  map_c_hat = automixfit(map_c_mcmc)
  
  #non-informative prior for current treatment
  prior_t = mixbeta(c(1,0,se_mu_t))
  
  #compute type I error rate and power. no need to update with current trial data
  design <- oc2S(prior_t, map_c_hat, n_t, n_c, rule)
  
  #tibble for type I error
  tt1 <- tibble(ctrl = ctrl_range, es = 0, prob = design(ctrl_range, ctrl_range), Type = "Type I error", Method = label)
  #tibble for power
  tt2 <- tibble(ctrl = ctrl_range, es = es, prob = design(ctrl_range + es, ctrl_range), Type = "Power", Method = label)
  
  rbind(tt1, tt2)
}


#generate a plot to compare Half-Normal priors: scale = 2 and 0.2
priordt <- rbind(tibble(Sigma = 2, sample = extraDistr::rhnorm(10000, sigma = 2)), tibble(Sigma = 0.2, sample = extraDistr::rhnorm(10000, sigma = 0.2)))
priordt %>% 
  ggplot(aes(x = sample, fill = factor(Sigma))) + geom_density(alpha = 0.5) + labs(fill = "HN\nScale") + 
  scale_x_continuous(breaks = seq(0, 6,  1), name = "") +
  theme(axis.title = element_text(face="bold",size=20),
        axis.text = element_text(size=20),
        legend.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.position = c(0.2, 0.5)
  )

#####################################################
#####################################################
#####################################################

#single criterion decision rule
rule_single <- decision2S(pc = c(0.975), qc = c(0), lower.tail = F)
HNscale_c <- 2

oc_hn2 <- oc_RBesT(rule = rule_single, histdt = histdt, label = "Borrowing", n_t, n_c, ctrl_range, es, se_mu_c, se_mu_t, HNscale_c)
oc_hn2_nb <- oc_RBesT(rule = rule_single, histdt = histdt_dummy, label = "No Borrowing", n_t, n_c, ctrl_range, es, se_mu_c, se_mu_t, HNscale_c)

plotdt <- rbind(oc_hn2, oc_hn2_nb, 
                tibble(ctrl = ctrl_range, es = es, prob = FPwr20$power, Type = "Power", Method = "Frequentist (20)")
                )
plotdt %>%
  ggplot(aes(x=ctrl,y=prob, color=Type, linetype=Method)) + geom_line(size = 1.5) +
  scale_x_continuous(breaks = seq(0.05, 0.35, 0.02), name = "True Current Control Rate") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), name = "Probability") +
  geom_vline(xintercept = 0.145, linetype=2, size=1) +
  geom_hline(yintercept = 0.025, linetype=2, size=1) +
  scale_color_manual(values = c("blue", "red")) +
  scale_linetype_manual(values = c("solid", "twodash", "dotted")) +
  theme(axis.title = element_text(face="bold",size=15),
        axis.text = element_text(size=12),
        legend.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=20), 
        legend.key.size = unit(1, "cm")
  )

#dual criteria decision rule
rule_dual <- decision2S(pc = c(0.975, 0.6), qc = c(0, 0.2), lower.tail = F)

oc_hn2 <- oc_RBesT(rule = rule_dual, histdt = histdt, label = "Borrowing", n_t, n_c, ctrl_range, es, se_mu_c, se_mu_t, HNscale_c)
oc_hn2_nb <- oc_RBesT(rule = rule_dual, histdt = histdt_dummy, label = "No Borrowing", n_t, n_c, ctrl_range, es, se_mu_c, se_mu_t, HNscale_c)

plotdt <- rbind(oc_hn2, oc_hn2_nb)

plotdt %>%
  ggplot(aes(x=ctrl,y=prob, color=Type, linetype=Method)) + geom_line(size = 1.5) +
  scale_x_continuous(breaks = seq(0.05, 0.35, 0.02), name = "True Current Control Rate") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), name = "Probability") +
  geom_vline(xintercept = 0.145, linetype=2, size=1) +
  geom_hline(yintercept = 0.025, linetype=2, size=1) +
  scale_color_manual(values = c("blue", "red")) +
  scale_linetype_manual(values = c("solid", "twodash")) +
  theme(axis.title = element_text(face="bold",size=15),
        axis.text = element_text(size=12),
        legend.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=20), 
        legend.key.size = unit(1, "cm")
  )

#####################################################################
#####################################################################
#####################################################################

HNscale_c <- 0.2

oc_hn0.2 <- oc_RBesT(rule = rule_dual, histdt = histdt, label = "Borrowing", n_t, n_c, ctrl_range, es, se_mu_c, se_mu_t, HNscale_c)
oc_hn0.2_nb <- oc_RBesT(rule = rule_dual, histdt = histdt_dummy, label = "No Borrowing", n_t, n_c, ctrl_range, es, se_mu_c, se_mu_t, HNscale_c)

plotdt <- rbind(oc_hn0.2, oc_hn0.2_nb)

plotdt %>%
  ggplot(aes(x=ctrl,y=prob, color=Type, linetype=Method)) + geom_line(size = 1.5) +
  scale_x_continuous(breaks = seq(0.05, 0.35, 0.02), name = "True Current Control Rate") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), name = "Probability") +
  geom_vline(xintercept = 0.145, linetype=2, size=1) +
  geom_hline(yintercept = 0.025, linetype=2, size=1) +
  scale_color_manual(values = c("blue", "red")) +
  scale_linetype_manual(values = c("solid", "twodash")) +
  theme(axis.title = element_text(face="bold",size=15),
        axis.text = element_text(size=12),
        legend.title=element_text(size=15,face="bold"),
        legend.text=element_text(size=20), 
        legend.key.size = unit(1, "cm")
  )
