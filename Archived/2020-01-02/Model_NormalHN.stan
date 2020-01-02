data {
  int<lower=0> H;          // number of historical trials 
  real y[H];              // historical mean
  int<lower=0> n[H];     //historical sample size
  real<lower=0> sigma;     // common historical standard deviation
  real<lower=0> HNscale_c;  // scale of half-normal hyper-prior on control
  real<lower=0> se_mu_c;  // scale of half-normal hyper-prior on control
}
parameters {
  real mu_c; 
  real<lower=0> tau;
  vector[H] epsilon;
}
transformed parameters {
  vector[H] theta;
  vector[H] se_yh;
  theta = mu_c + epsilon;
  se_yh = sigma/sqrt(n);
}
model {
  target += normal_lpdf(y | theta, se_yh);
  target += normal_lpdf(mu_c | 0, se_mu_c);
  target += normal_lpdf(epsilon | 0, tau);
  target += normal_lpdf(tau | 0, HNscale_c);
}
