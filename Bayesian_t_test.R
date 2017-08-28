library(tidyverse)
library(ggplot2)

N_01   <- 1000
mu_01  <- 100
sig_01 <- 10

N_02   <- 1000
mu_02  <- 120
sig_02 <- 10

set.seed(123)
Pop_01 <- rnorm(N_01, mu_01, sig_01)
Pop_02 <- rnorm(N_02, mu_02, sig_02)

Dat_Sample <- Pop_01 %>%
   as_data_frame() %>%
   mutate("Group" = "01") %>%
   bind_rows(as_data_frame(Pop_02) %>% mutate("Group" = "02"))

# ggplot(Dat_Sample, aes(value, fill = Group)) +
#    geom_histogram(alpha = 0.5, position = "identity") +
#    # facet_wrap(~ Group, ncol = 1) +
#    theme_bw()

likelihood <- function(pars) {
   mu_01  <- pars[1]
   sig_01 <- pars[2]
   mu_02  <- pars[3]
   sig_02 <- pars[4]
   sum(
      log(dnorm(Pop_01, mu_01, sig_01)),
      log(dnorm(Pop_02, mu_02, sig_02))
   )
   # prod(
   #    dnorm(Pop_01, mu_01, sig_01),
   #    dnorm(Pop_02, mu_02, sig_02)
   # )
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
      log(dunif(mu_01, -10000, 10000)),
      log(dunif(mu_02, -10000, 10000)),
      log(dunif(sig_01, -10000, 10000)),
      log(dunif(sig_02, -10000, 10000))
   )
   # sum(
   #    log(dnorm(mu_01, mean(Dat_Sample$value), 1000 * sd(Dat_Sample$value))),
   #    log(dnorm(mu_02, mean(Dat_Sample$value), 1000 * sd(Dat_Sample$value))),
   #    log(dexp(sig_01, rate = 0.1)),
   #    log(dexp(sig_02, rate = 0.1))
   # )
   # sum(
   #    log(dnorm(mu_01, mean(Dat_Sample$value), 1000 * sd(Dat_Sample$value))),
   #    log(dnorm(mu_02, mean(Dat_Sample$value), 1000 * sd(Dat_Sample$value))),
   #    log(dinvgamma(sig_01^2, 10)),
   #    log(dinvgamma(sig_02^2, 10))
   # )
   # prod(
   #    dnorm(mu_01, mean(Dat_Sample$value), 1000 * sd(Dat_Sample$value)),
   #    dnorm(mu_02, mean(Dat_Sample$value), 1000 * sd(Dat_Sample$value)),
   #    dexp(sig_01, rate = 0.1),
   #    dexp(sig_02, rate = 0.1)
   # )
   # dinvgamma(sig_01^2, 10, 1) *
      # dinvgamma(sig_02^2, 10, 1)
      # log(dunif(sig_01^2, -10000, 10000)),
      # log(dunif(sig_02^2, -10000, 10000)))
}

# posterior <- function(pars) {
#    likelihood(pars) * prior(pars)
# }

posterior <- function(pars) {
   likelihood(pars) + prior(pars)
}


N_iter <- 2000
pars <- c(mean(Pop_01), sd(Pop_01), mean(Pop_02), sd(Pop_02))
results <- matrix(0, nrow = N_iter, ncol = 8)
results[1, ] <- c(pars, 0, 0, 0, 0)

for (i in 2:N_iter) {
   for (j in 1:4) {
      candidate <- pars
      # candidate[c(1, 3)] <- candidate[c(1, 3)] + rnorm(2, sd = 0.01)
      # candidate[c(2, 4)] <- candidate[c(2, 4)] + rnorm(2, sd = 0.01)
      candidate[j] <- candidate[j] + rnorm(1, sd = 0.5)
      # ratio <- posterior(candidate) / posterior(pars)
      # if (ratio > 1) {
      #    pars <- candidate
      # }
      ratio <- posterior(pars) / posterior(candidate)
      t <- runif(1, 0.5, 1.5)
      if (t < ratio) { pars <- candidate }
   }
   results[i, ] <- c(pars, ratio, t, posterior(pars), posterior(candidate))
}

apply(results, 2, mean)
sum(results[, 5] > results[, 6])

par(mfrow = c(2, 2))
plot(results[seq(500, N_iter, 10), 1])
plot(results[seq(500, N_iter, 10), 2])
plot(results[seq(500, N_iter, 10), 3])
plot(results[seq(500, N_iter, 10), 4])

hist(results[seq(500, N_iter, 10), 1])
hist(results[seq(500, N_iter, 10), 2])
hist(results[seq(500, N_iter, 10), 3])
hist(results[seq(500, N_iter, 10), 4])


plot(results[seq(500, N_iter, 5), 1])
plot(results[501:N_iter, 2])
plot(results[501:N_iter, 3])
plot(results[501:N_iter, 4])
plot(results[501:N_iter, 5])
plot(results[501:N_iter, 6])
plot(results[501:N_iter, 7])
plot(results[501:N_iter, 8])

apply(results, 2, mean)


Res_MCMC <- results[501:N_iter, ]




plot(Res_MCMC[1:(i-1), 1])
plot(Res_MCMC[1:(i-1), 2])
plot(Res_MCMC[1:(i-1), 3])
plot(Res_MCMC[1:(i-1), 4])

plot(results[1:(i-1), 1])
plot(results[1:(i-1), 2])
plot(results[1:(i-1), 3])
plot(results[1:(i-1), 4])

