data {
  int N;
  int y[N];
}

parameters {
  real lambda;
}

model {
  y ~ poisson(lambda);
}
