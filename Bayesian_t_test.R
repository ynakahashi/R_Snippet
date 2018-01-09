library(tidyverse)
library(ggplot2)

N_01   <- 100
mu_01  <- 100
sig_01 <- 10

N_02   <- 100
mu_02  <- 103
sig_02 <- 10

set.seed(123)
Pop_01 <- rnorm(N_01, mu_01, sig_01)
Pop_02 <- rnorm(N_02, mu_02, sig_02)

Pop_01 %>%
   as_data_frame() %>%
   mutate("Group" = "01") %>%
   bind_rows(as_data_frame(Pop_02) %>% mutate("Group" = "02")) %>% 
   ggplot(aes(value, fill = Group)) +
   geom_histogram(alpha = 0.5, position = "identity") +
   theme_bw()

t.test(Pop_01, Pop_02)

likelihood <- function(pars) {
   mu_01  <- pars[1]
   sig_01 <- pars[2]
   mu_02  <- pars[3]
   sig_02 <- pars[4]
   sum(
      dnorm(Pop_01, mu_01, sig_01, log = TRUE),
      dnorm(Pop_02, mu_02, sig_02, log = TRUE)
   )
}


pars <- c(100, 10, 100, 20)
optim(par = optim(par = pars, fn = likelihood, control=list(fnscale = -1))$par,
      fn = likelihood, control=list(fnscale = -1))

prior <- function(pars) {
   mu_01  <- pars[1]
   sig_01 <- pars[2]
   mu_02  <- pars[3]
   sig_02 <- pars[4]
   sum(
      dunif(mu_01, -10000, 10000, log = TRUE),
      dunif(mu_02, -10000, 10000, log = TRUE),
      # dnorm(mu_01, mean(Pop), 1000*sd(Pop), log = TRUE),
      # dnorm(mu_02, mean(Pop), 1000*sd(Pop), log = TRUE),
      dunif(sig_01, -10000, 10000, log = TRUE),
      dunif(sig_02, -10000, 10000, log = TRUE)
      # dexp(sig_01, rate = 0.1, log = TRUE),
      # dexp(sig_02, rate = 0.1, log = TRUE)
   )
}


posterior <- function(pars) {
   likelihood(pars) + prior(pars)
}


N_iter <- 30000
pars <- c(50, 30, 20, 50)
results <- matrix(0, nrow = N_iter, ncol = 8)
results[1, ] <- c(pars, 0, 0, 0, 0)


candidate <- pars
system.time(
for (i in 2:N_iter) {
   for (j in 1:4) { # Gibbs Sampling
      # candidate <- pars
      proposal <- candidate
      # candidate <- candidate + rnorm(4, sd = 0.5) # Metropolis
      proposal[j] <- candidate[j] + rnorm(1, sd = 0.5) # Gibbs Sampling
      ratio <- max(0, exp(posterior(proposal) - posterior(candidate)))
      t <- runif(1)
      if (t < ratio) { candidate <- proposal }
   } # Gibbs Sampling
   results[i, ] <- c(candidate, ratio, t, posterior(proposal), posterior(candidate))
}
)



par(mfrow = c(2, 2))
plot(results[, 1])
plot(results[, 2])
plot(results[, 3])
plot(results[, 4])


Res_MCMC <- results[seq(15000, N_iter, 10), ]
summary(Res_MCMC[, c(1:4)])
apply(Res_MCMC, 2, mean)
# sum(Res_MCMC[, 5] > Res_MCMC[, 6])


hist(Res_MCMC[, 1])
hist(Res_MCMC[, 2])
hist(Res_MCMC[, 3])
hist(Res_MCMC[, 4])

par(mfrow = c(1, 1))
hist(Res_MCMC[, 3] - Res_MCMC[, 1])
quantile(Res_MCMC[, 1] - Res_MCMC[, 3], seq(0, 1, 0.025))
mean(Res_MCMC[, 3] - Res_MCMC[, 1] < 0)

