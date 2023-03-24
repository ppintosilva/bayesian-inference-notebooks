library(tidyverse)
library(rstan)

discoveries <- read_csv(
  file = "discoveries/evaluation_discoveries.csv",
  col_types = list(
    time = "i",
    discoveries = "i"
  )
)

discoveries %>%
  ggplot() +
  geom_line(
    aes(x = time, y = discoveries)
  )

set.seed(1)

fit <- rstan::stan(
  file = "discoveries/discoveries.stan",
  chains = 4,
  iter = 1000,
  warmup = 500,
  data = list(
    N = nrow(discoveries),
    y = discoveries$discoveries
  )
)

lambda_draws <- extract(fit, pars = "lambda")[[1]] %>%
  matrix(ncol = 4, nrow = 500)

Rhat(lambda_draws)
ess_bulk(lambda_draws)/4
ess_tail(lambda_draws)/4

lambda_draws %>% as_tibble() %>% set_names(c("c1", "c2", "c3", "c4"))
 
# 16.1.7 - Model converged fine

