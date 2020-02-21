tt = rbeta(100000,0.5,0.5)
plot(density(tt))

quantile(tt,c(0.025, 0.5, 0.975))

ilogit = function(x) exp(x)/(1+exp(x))
lgt = rnorm(100000,0,2)
plot(density(ilogit(lgt)))


quantile(ilogit(lgt),c(0.025, 0.5, 0.975))
