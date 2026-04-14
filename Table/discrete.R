source("derivation.R")

library(assertthat)
library(tidyverse)
library(logitnorm)

# helper
auc_to_d <- function(auc) {
  qnorm(auc) * sqrt(2)
}

# n samples
n <- 1e5

# AI
theta_m <- 0.5
d_m <- 5

# Human
theta_h <- 0.7
d_h <- 1

bins <- 3


# Monte Carlo samples
sample_mc <- function(n, theta, d) {
  y <- runif(n) < theta
  p <- ifelse(y,
    rlogitnorm(n, d / 2, 1),
    rlogitnorm(n, -d / 2, 1)
  )
  list(
    p = p,
    y = y
  )
}

# get logit posterior per bin
get_llr <- function(p, y, bins) {
  cf <- tibble(p, y) |>
    mutate(bin = cut(p, breaks = bins, include.lowest = TRUE)) |>
    group_by(bin) |>
    mutate(pos = mean(y)) |>
    ungroup()

  cf$pos
}

# simulate combined accuracy
sim_com_acc <- function(n, theta_m, theta_h, d_m, d_h, bins) {
  sm <- sample_mc(n, theta_m, d_m)
  sh <- sample_mc(n, theta_h, d_h)

  ym <- sm$y
  yh <- sh$y

  pm <- sm$p
  ph <- sh$p

  am <- get_llr(pm, ym, bins)
  ah <- get_llr(ph, yh, bins)

  obs <- mean(ifelse(am > ah, ym, yh))
  obs
}


# get mean absolute error when simulation uses discrete confidence binds
get_mae <- function(n, theta_m, theta_h, d_m, d_h, bins) {
  # theoretical combined accuracy
  pre <- combined_accuracy_d(theta_m = theta_m, theta_h = theta_h, d_m = d_m, d_h = d_h)

  # simulate combined accuracy
  obs <- sim_com_acc(n, theta_m, theta_h, d_m, d_h, bins)

  list(pre = pre, obs = obs, mae = abs(obs - pre), pct = abs(obs - pre) / pre * 100)
}

get_mae(n, theta_m, theta_h, d_m, d_h, bins)

df <- tribble(
  ~n, ~theta_m, ~theta_h, ~d_m, ~d_h,
  n, 0.6, 0.9, 2, 1,
  n, 0.6, 0.9, 1, 2,
  n, 0.3, 0.6, 2, 1,
  n, 0.3, 0.6, 1, 2,
  n, 0.2, 0.8, 2, 1,
  n, 0.2, 0.8, 1, 2,
  n, 0.2, 0.4, 2, 1,
  n, 0.2, 0.4, 1, 2,
  n, 0.6, 0.9, 1, 3,
  n, 0.6, 0.9, 3, 1,
  n, 0.3, 0.6, 1, 3,
  n, 0.3, 0.6, 3, 1,
  n, 0.2, 0.8, 3, 1,
  n, 0.2, 0.8, 1, 3,
  n, 0.2, 0.4, 3, 1,
  n, 0.2, 0.4, 1, 3,
  n, 0.6, 0.9, 1, .5,
  n, 0.6, 0.9, .5, 1,
  n, 0.3, 0.6, 1, .5,
  n, 0.3, 0.6, .5, 1,
  n, 0.2, 0.8, .5, 1,
  n, 0.2, 0.8, 1, .5,
  n, 0.2, 0.4, .5, 1,
  n, 0.2, 0.4, 1, .5,
)

for (b in c(3, 5, 10, 50)) {
  ba <- df |>
    mutate(bins = b) |>
    rowwise() |>
    mutate(res = list(get_mae(n, theta_m, theta_h, d_m, d_h, bins))) |>
    unnest_wider(res) |>
    arrange(desc(mae)) |>
    glimpse()
}