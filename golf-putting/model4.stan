// J observations of n (trials), y (successes) for each distance value x
data {
  int<lower=0> J;
  int n[J];
  vector[J] f;
  vector[J] x;
  real r;
  real R;
}

// Parametrised by sigma_ang and sigma_dist
parameters {
  real<lower=0> sigma_angle;
  real<lower=0> sigma_distance;
  real<lower=0> sigma_y;  
}

model {
  vector[J] p;
  vector[J] s;
  
  sigma_angle ~ normal(0,1);    // half-normal prior (because parameter value is constrained above with lower=0)
  sigma_distance ~ normal(0,1); // half-normal prior (because parameter value is constrained above with lower=0)
  sigma_y ~ normal(0,1); // half-normal prior (because parameter value is constrained above with lower=0)
  
  for (j in 1:J) {
    p[j] = (2*Phi(asin((R-r)/x[j]) / sigma_angle) - 1) *
      (Phi( 2 / ((x[j]+1) * sigma_distance)) -
       Phi(-1 / ((x[j]+1) * sigma_distance)));
    
    s[j] = sqrt(p[j]*(1-p[j]) / n[j] + sigma_y^2);
    
    target += normal_lpdf(f[j] | p[j], s[j]);
  }

}
