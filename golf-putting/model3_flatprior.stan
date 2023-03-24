// J observations of n (trials), y (successes) for each distance value x
data {
  int<lower=0> J;
  int n[J];
  int y[J];
  vector[J] x;
  real r;
  real R;
}

// Parametrised by sigma_ang and sigma_dist
parameters {
  real sigma_angle;
  real sigma_distance;
}

model {
  vector[J] p;
  
  for (j in 1:J) {
    p[j] = (2*Phi(asin((R-r)/x[j]) / sigma_angle) - 1) * (Phi(2 / ((x[j]+1) * sigma_distance)) - Phi(-1 / ((x[j]+1) * sigma_distance)));
  }
  
  y ~ binomial(n, p);
}

generated quantities {
  real sigma_degrees;
  //vector[J] y_rep;
  //vector[J] p;
  
  sigma_degrees = (180/pi())*sigma_angle;
  
  // for (j in 1:J) {
  //   p[j] = 2*Phi(asin((R-r)/x[j])/sigma) - 1;
  //   y_rep[j] = binomial_rng(n[j], p[j]);
  // }
}