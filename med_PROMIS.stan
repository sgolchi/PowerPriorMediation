data{
  int N;
  real y[N];
  real m[N];
  real x[N];
  real z[N];
  real w1[N]; 
  real w2[N];
  real<lower = 0, upper = 1> a0[N];   
}
parameters {
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  real beta5;
  real beta6;
  real i2;
  real i3;
  real a; 
  real b;
  real cp;
  real<lower = 0> sig2;
  real<lower = 0> sig3;
}
model {
  a ~ normal(0, 1000);
  b ~ normal(0, 1000);
  cp ~ normal(0, 1000);
  sig2 ~ cauchy(0, 1000);
  sig3 ~ cauchy(0, 1000);
    for (n in 1:N) {
      target += a0[n] * normal_lpdf(m[n] | i2 + a*x[n] + beta1*z[n] + beta2*w1[n] + beta3*w2[n], sig2);
      target += a0[n] * normal_lpdf(y[n] | i3 + cp*x[n] + b*m[n] + beta4*z[n] + beta5*w1[n] + beta6*w2[n], sig3);
    }
}
generated quantities{
  real med_eff = a*b;
}
  
