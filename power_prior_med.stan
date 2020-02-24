data{
  int N;                                          // number of data rows
  int K;                                          // number of covariates other than the exposure and mediator
  real y[N];                                      // outcome
  real M[N];                                      // mediator
  real x[N];                                      // exposure
  matrix[N, K] Z;                                 // covariates
  real<lower = 0, upper = 1> w[N];                // weights
}
parameters {
  vector[K] beta1;                                // covariate coeeficient vector in the mediator model
  vector[K] beta2;                                // covariate coefficient vector in teh outcome model
  real i2;                                        // intercept of the mediator model
  real i3;                                        // intercept of the outcome model
  real a;                                         // a-path parameter
  real b;                                         // b-path parameter
  real c;                                         // c-path parameter
  real<lower = 0> sig2;                           // variance parameter of the mediator model
  real<lower = 0> sig3;                           // variance parameter of teh outcome model
}
model {
  a ~ normal(0, 1000);
  b ~ normal(0, 1000);
  c ~ normal(0, 1000);
  sig2 ~ cauchy(0, 1000);
  sig3 ~ cauchy(0, 1000);
    for (n in 1:N) {
      target += w[n] * normal_lpdf(M[n] | i2 + a*x[n] + Z[n]*beta1, sig2);
      target += w[n] * normal_lpdf(y[n] | i3 + c*x[n] + b*M[n] + Z[n]*beta2, sig3);
    }
}
generated quantities{
  real med_eff = a*b;
}
  
