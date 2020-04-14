data{
  int N;
  real y[N];
  real M[N];
  real x[N];
  real z[N];
  real W1[N]; 
  real<lower = 0, upper = 1> w[N];   
}
parameters {
  real beta1;
  real beta2;
  real beta3;
  real beta4;
  real i2;
  real i3;
  real a; 
  real b;
  real c;
  real<lower = 0> sig2;
  real<lower = 0> sig3;
}
model {
  a ~ normal(0, 1000);
  b ~ normal(0, 1000);
  c ~ normal(0, 1000);
  sig2 ~ cauchy(0, 1000);
  sig3 ~ cauchy(0, 1000);
    for (n in 1:N) {
      target += w[n] * normal_lpdf(M[n] | i2 + a*x[n] + beta1*z[n] + beta2*W1[n] , sig2);
      target += w[n] * normal_lpdf(y[n] | i3 + c*x[n] + b*M[n] + beta3*z[n] + beta4*W1[n] , sig3);
    }
}
generated quantities{
  real med_eff = a*b;
}
  
