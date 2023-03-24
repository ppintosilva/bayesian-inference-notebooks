library(tidyverse)
library(rstan)
library(bayesplot)

golf <- tibble::tribble(
  ~x, ~n ,~y,
  2 ,1443, 1346,
  3 ,694 ,577,
  4 ,455 ,337,
  5 ,353 ,208,
  6 ,272 ,149,
  7 ,256 ,136,
  8 ,240 ,111,
  9 ,217 ,69,
  10, 200, 67,
  11, 237, 75,
  12, 202, 52,
  13, 192, 46,
  14, 174, 54,
  15, 167, 28,
  16, 201, 27,
  17, 195, 31,
  18, 191, 33,
  19, 147, 20,
  20, 152, 24
)

model1 <- "
// J observations of n (trials), y (successes) for each distance value x
data {
  int<lower=0> J;
  int n[J];
  int y[J];
  vector[J] x;
}

// Parametrised by a and b
parameters {
  real a;
  real b;
}


model {
  vector[J] p;
  
  for (j in 1:J) {
    p[j] = inv_logit(a + b*x[j]);
  }
  
  y ~ binomial(n, p);
}

generated quantities {
  vector[J] y_rep;
  vector[J] p;
  
  for (j in 1:J) {
    p[j] = inv_logit(a + b*x[j]);
    y_rep[j] = binomial_rng(n[j], p[j]);
  }
}
"

fit1 <- stan(
  model_code = model1,
  chains = 4,
  iter = 2000,
  data = list(
    J = nrow(golf), x = golf$x, y = golf$y, n = golf$n
  )
)

fit1_yrep <- extract(fit1, pars = "y_rep")[[1]]

fit1_draws <- posterior::as_draws(fit1_yrep)

fit1_summary <- posterior::summarise_draws(fit1_draws) %>%
  select(mean, sd) %>%
  rename(ypred_mean = mean, ypred_sd = sd)


fit1_params <- extract(fit1, pars = c("a", "b")) %>% as_tibble()

fit1_params_summary <- fit1_params %>%
  summarise(
    across(everything(), mean)
  )

logistic <- function(x, a, b) 1/(1 + exp(-a - b*x))

bind_cols(
  golf,
  fit1_summary
) %>%
  mutate(
    across(starts_with("y"), ~ .x/n)
  ) %>%
  ggplot() +
  geom_point(
    aes(x = x, y = y, color = "data"),
    size = 4
  ) +
  geom_function(
    fun = logistic,
    args = list(a = fit1_params_summary$a, b = fit1_params_summary$b),
    aes(color = "LR model")
  ) +
  geom_pointrange(
    aes(x = x, 
        y = ypred_mean, 
        ymin = ypred_mean - ypred_sd,
        ymax = ypred_mean + ypred_sd, 
        color = "LR model"
    ),
    size = .50
  ) +
  xlab("Distance from hole (feet)") +
  ylab("Probability of success") +
  xlim(0,22) +
  ylim(0,1) +
  scale_color_manual(
    values = c(
      "data" = "blue",
      "LR model" = "red"
    )
  )
  theme_bw()


# model 2 - geometrical model

model2 <- "
// J observations of n (trials), y (successes) for each distance value x
data {
  int<lower=0> J;
  int n[J];
  int y[J];
  vector[J] x;
  real r;
  real R;
}

// Parametrised by a and b
parameters {
  real sigma;
}

model {
  vector[J] p;
  
  for (j in 1:J) {
    p[j] = 2*Phi(asin((R-r)/x[j]) / sigma) - 1;
  }
  
  y ~ binomial(n, p);
}

generated quantities {
  real sigma_degrees;
  vector[J] y_rep;
  vector[J] p;
  
  sigma_degrees = (180/pi())*sigma;
  
  for (j in 1:J) {
    p[j] = 2*Phi(asin((R-r)/x[j])/sigma) - 1;
    y_rep[j] = binomial_rng(n[j], p[j]);
  }
}
"

r_feet <- 1.68/2 * 0.0833 # inches to feet
R_feet <- 4.25/2 * 0.0833 # inches to feet

fit2 <- stan(
  model_code = model2,
  chains = 4,
  iter = 2000,
  data = list(
    J = nrow(golf), x = golf$x, y = golf$y, n = golf$n, r = r_feet, R = R_feet # to feet
  )
)

fit2_quant <- extract(fit2)

fit2_draws <- posterior::as_draws(fit2_quant$y_rep)

fit2_summary <- posterior::summarise_draws(fit2_draws) %>%
  select(mean, sd) %>%
  rename(ypred_mean = mean, ypred_sd = sd)

fit2_params <- bind_cols(
  fit2_quant$sigma %>% tibble(sigma = .),
  fit2_quant$sigma_degrees %>% tibble(sigma_degrees = .)
)

fit2_params_summary <- fit2_params %>%
  summarise(
    across(everything(), mean)
  )

fit2_params_summary

pr_angle <- function(x, sigma, r = r_feet, R = R_feet) 2*pnorm(asin((R-r)/x)/sigma) - 1

bind_cols(
  golf
) %>%
  mutate(
    across(starts_with("y"), ~ .x/n)
  ) %>%
  ggplot() +
  geom_point(
    aes(x = x, y = y, color = "data"),
    size = 4
  ) +
  # geom_function(
  #   fun = logistic,
  #   args = list(a = fit1_params_summary$a, b = fit1_params_summary$b),
  #   aes(color = "LR model")
  # ) +
  geom_function(
    fun = pr_angle,
    args = list(sigma = fit2_params_summary$sigma),
    aes(color = "Geom model")
  ) +
  xlab("Distance from hole (feet)") +
  ylab("Probability of success") +
  xlim(0,22) +
  ylim(0,1) +
  scale_color_manual(
    values = c(
      "data" = "blue",
      "LR model" = "red",
      "Geom model" = "black"
    )
  ) +
theme_bw()


# model 3 - new data

newgolf <- tibble::tribble(
  ~x, ~n, ~y,
  0.28 ,45198 , 45183,
  0.97 ,183020, 182899,
  1.93 ,169503, 168594,
  2.92 ,113094, 108953,
  3.93 ,73855 , 64740,
  4.94 ,53659 , 41106,
  5.94 ,42991 , 28205,
  6.95 ,37050 , 21334,
  7.95 ,33275 , 16615,
  8.95 ,30836 , 13503,
  9.95 ,28637 , 11060,
  10.95, 26239, 9032,
  11.95, 24636, 7687,
  12.95, 22876, 6432,
  14.43, 41267, 9813,
  16.43, 35712, 7196,
  18.44, 31573, 5290,
  20.44, 28280, 4086,
  21.95, 13238, 1642,
  24.39, 46570, 4767,
  28.40, 38422, 2980,
  32.39, 31641, 1996,
  36.39, 25604, 1327,
  40.37, 20366, 834,
  44.38, 15977, 559,
  48.37, 11770, 311,
  52.36, 8708 , 231,
  57.25, 8878 , 204,
  63.23, 5492 , 103,
  69.18, 3087 , 35,
  75.19, 1742 , 24
)


# mod3 <- cmdstanr::cmdstan_model("model3.stan", compile = FALSE)
# mod3$check_syntax()

fit3 <- stan(
  file = "golf-putting/model3_flatprior.stan",
  chains = 4,
  iter = 2000,
  data = list(
    J = nrow(newgolf),
    x = newgolf$x,
    y = newgolf$y,
    n = newgolf$n,
    r = r_feet,
    R = R_feet # to feet
  )
)

np <- nuts_params(fit3)

head(log_posterior(fit3))

color_scheme_set("darkgray")

mcmc_parcoord(as.array(fit3) %>% head(), np = )

mcmc_trace(as.array(fit3), pars = c("sigma_angle", "sigma_distance"), np = np) +
  xlab("Post-warmup iteration")

fit3_weak <- stan(
  file = "golf-putting/model3_weakprior.stan",
  chains = 4,
  iter = 2000,
  data = list(
    J = nrow(newgolf),
    x = newgolf$x,
    y = newgolf$y,
    n = newgolf$n,
    r = r_feet,
    R = R_feet # to feet
  )
)

np <- nuts_params(fit3_weak)

mcmc_trace(as.array(fit3_weak), pars = c("sigma_angle", "sigma_distance"), np = np) +
  xlab("Post-warmup iteration")

fit3_quant <- extract(fit3_weak)

fit3_params <- bind_cols(
  fit3_quant$sigma_angle %>% tibble(sigma_angle = .),
  fit3_quant$sigma_distance %>% tibble(sigma_distance = .)
)

fit3_params_mean <- 
  fit3_params %>%
  mutate(chain = ((row_number() - 1) %/% 1000) + 1) %>%
  group_by(chain) %>%
  summarise(
    across(everything(), .fns = list(mean = mean, sd  = sd))
  ) %>%
  pivot_longer(cols = 2:5, names_to = c("parameter", "statistic"), values_to = "value", names_pattern = "sigma_([[:alnum:]]+)_([[:alnum:]]+)") %>%
  pivot_wider(names_from = "statistic", values_from = "value") %>%
  filter(sd > 0.00001) %>%
  filter(chain == 1) %>%
  pull(mean) %>%
  set_names(c("angle", "distance"))

fit3_params_mean2 <- 
  fit3_params %>%
  summarise(
    across(everything(), .fns = list(mean = mean))
  )
  
pr_model3 <- function(x, sigma_angle, sigma_distance,
                      r = r_feet, R = R_feet) {
  # term 1
  (2*pnorm(asin((R-r)/x)/sigma_angle) - 1) * 
    # term 2
    (pnorm(2/((x + 1)*sigma_distance)) - pnorm(-1/((x + 1)*sigma_distance)))
}

newgolf %>%
  mutate(
    across(starts_with("y"), ~ .x/n)
  ) %>%
  ggplot() +
  geom_point(
    aes(x = x, y = y),
    size = 2,
    color = "red"
  ) +
  geom_function(
    fun = pr_model3,
    args = list(sigma_angle = fit3_params_mean2$sigma_angle_mean,# fit3_params_mean["angle"],
                sigma_distance = fit3_params_mean2$sigma_distance_mean), # fit3_params_mean["distance"]),
    color = "red"    
  ) +
  xlab("Distance from hole (feet)") +
  ylab("Probability of success") +
  xlim(0,80) +
  ylim(0,1) +
  theme_bw()


mod4 <- cmdstanr::cmdstan_model("golf-putting/model4.stan", compile = FALSE)
mod4$check_syntax()

fit4 <- stan(
  file = "golf-putting/model4.stan",
  chains = 4,
  iter = 2000,
  data = list(
    J = nrow(newgolf),
    x = newgolf$x,
    f = newgolf$y/newgolf$n,
    n = newgolf$n,
    r = r_feet,
    R = R_feet # to feet
  )
)


fit4_ext <- extract(fit4)

fit4_params <- bind_cols(
  fit4_ext$sigma_angle %>% tibble(sigma_angle = .),
  fit4_ext$sigma_distance %>% tibble(sigma_distance = .),
  fit4_ext$sigma_y %>% tibble(sigma_y = .)
)

fit4_means <- 
  fit4_params %>%
  summarise(
    across(everything(), .fns = list(mean = mean))
  )

traceplot(fit4)

stan_rhat(fit4, bins = 12)
stan_ess(fit4, bins = 12)
stan_mcse(fit4, bins = 12)
stan_par(fit4, par = "sigma_distance")
stan_par(fit4, par = "sigma_angle")
stan_par(fit4, par = "sigma_y")

newgolf %>%
  mutate(
    across(starts_with("y"), ~ .x/n)
  ) %>%
  ggplot() +
  geom_point(
    aes(x = x, y = y),
    size = 2,
    color = "red"
  ) +
  geom_function(
    fun = pr_model3,
    args = list(sigma_angle = fit4_means$sigma_angle_mean,
                sigma_distance = fit4_means$sigma_distance_mean),
    color = "red"    
  ) +
  xlab("Distance from hole (feet)") +
  ylab("Probability of success") +
  xlim(0,80) +
  ylim(0,1) +
  theme_bw()
