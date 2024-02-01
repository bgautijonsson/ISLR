library(tidyverse)
library(purrr)
library(geomtextpath)
library(scales)
library(patchwork)
library(ggh4x)
theme_set(theme_classic())

#### Define function for calculating in and out-of-sample deviance ####
get_logliks <- function(i) {
  n_obs <- 200
  set.seed(i)
  n_pars <- sample(1:12, size = 1)
  beta <- rnorm(n_pars) 
  
  X_train <- matrix(
    rnorm(n_obs * n_pars),
    ncol = n_pars,
    nrow = n_obs
  )
  X_test <- matrix(
    rnorm(n_obs * n_pars),
    ncol = n_pars,
    nrow = n_obs
  )
  y_train <- X_train %*% beta + rnorm(n_obs)
  y_test <- X_test %*% beta + rnorm(n_obs)
  
  m <- lm(y_train ~ X_train)
  
  y_hat_train <- X_train %*% coef(m)[-1] + coef(m)[1]
  y_hat_test <- X_test %*% coef(m)[-1] + coef(m)[1]
  
  sigma <- sqrt(summary(m)$sigma)
  
  
  loglik_train <- dnorm(y_train, mean = y_hat_train, sd = sigma, log = TRUE) |> 
    sum()
  
  loglik_test <- dnorm(y_test, mean = y_hat_test, sd = sigma, log = TRUE) |> 
    sum()
  
  tibble(
    iter = i,
    n_pars = n_pars,
    train = loglik_train,
    test = loglik_test
  )
}



d <- map_dfr(1:5000, get_logliks)



p1 <- d |> 
  group_by(n_pars) |> 
  summarise_at(
    vars(train, test),
    \(x) mean(-2 * x)
  ) |> 
  ggplot(aes(n_pars + 2, train)) +
  geom_point(size = 3) +
  geom_labelline(aes(label = "Training"), linewidth = 1, size = 5) +
  geom_point(aes(y = test), size = 3) +
  geom_labelline(aes(y = test, label = "Testing"), linewidth = 1, size = 5) +
  geom_textsegment(
    aes(xend = n_pars + 2, yend = test, label = round(test - train)),
    lty = 2,
    linewidth = 0.2
  ) +
  scale_x_continuous(
    breaks = 2:15,
    guide = guide_axis_truncated()
  ) +
  scale_y_continuous(
    limits = c(550, 590),
    guide = guide_axis_truncated()
  ) +
  labs(
    x = "Number of parameters",
    y = "Deviance"
  )

p2 <- d |> 
  mutate(diff = -2 * (test - train)) |> 
  summarise(
    mean = mean(diff),
    .by = n_pars
  ) |> 
  ggplot(aes(n_pars + 2, mean)) +
  geom_textabline(
    label = "2 x (number of parameters)",
    intercept = 0, 
    slope = 2,
    hjust = 0.13
    ) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = 0:15,
    limits = c(0, NA),
    # expand = expansion(),
    guide = guide_axis_truncated()
  ) +
  scale_y_continuous(
    limits = c(0, 31),
    breaks = breaks_extended(7),
    # expand = expansion(),
    guide = guide_axis_truncated()
  ) +
  labs(
    x = "Number of parameters",
    y = expression(D[test] - D[train])
  )

p <- p1 + p2 +
  plot_annotation(
    title = "Simulating how well the AIC predicts the out-of-sample deviance for linear models with Gaussian residuals",
    subtitle = "The AIC assumes that out-of-sample deviance is equal to in-sample deviance plus twice the number of parameters",
    caption = "Plots show averages from 5.000 simulated datasets with 200 observations in the training and testing sets",
    theme = theme(
      plot.title = element_text(face = "bold")
    )
  )

p

ggsave(
  plot = p,
  filename = "Figures/AIC.png",
  width = 8, height = 0.5 * 8, scale = 1.6
)

