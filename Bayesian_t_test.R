library(tidyverse)
library(ggplot2)

N_01   <- 100
mu_01  <- 20
sig_01 <- 10

N_02   <- 100
mu_02  <- 50
sig_02 <- 10

# set.seed(123)
Pop_01 <- rnorm(N_01, mu_01, sig_01)
Pop_02 <- rnorm(N_02, mu_02, sig_02)

# Dat_Sample <- Pop_01 %>%
#    as_data_frame() %>% 
#    mutate("Group" = "01") %>% 
#    bind_rows(as_data_frame(Pop_02) %>% mutate("Group" = "02"))
# 
# ggplot(Dat_Sample, aes(value, fill = Group)) +
#    geom_histogram(alpha = 0.5, position = "identity") +
#    # facet_wrap(~ Group, ncol = 1) +
#    theme_bw()

# likelihood <- function(pars) {
#    mu_01  <- pars[1] 
#    sig_01 <- pars[2]
#    mu_02  <- pars[3] 
#    sig_02 <- pars[4]
#    -sum(
#       log(dnorm(Pop_01, mu_01, sig_01)),
#       log(dnorm(Pop_02, mu_02, sig_02)))
# }
# 
# 
# pars <- c(0, 10, 10, 20)
# optim(par = optim(par = pars, fn = likelihood)$par, fn = likelihood)

prior <- function(pars) {
   mu_01  <- pars[1]
   sig_01 <- pars[2]
   mu_02  <- pars[3]
   sig_02 <- pars[4]
   dnorm(mu_01, mean(Dat_Sample$value), sd(Dat_Sample$value)) *
      dnorm(mu_02, mean(Dat_Sample$value), sd(Dat_Sample$value)) *
      dexp(sig_01, rate = 0.1) *
      dexp(sig_02, rate = 0.1)
}
# 
# posterior <- function(pars) {
#    likelihood(pars) * prior(pars)
# }


pars <- c(mean(Pop_01), sd(Pop_01), mean(Pop_02), sd(Pop_02))
N_iter <- 10000
results <- matrix(0, nrow = N_iter, ncol = 4)
results[1, ] <- pars
for (i in 2:N_iter){
   candidate <- pars
   candidate[c(1, 3)] <- candidate[c(1, 3)] + rnorm(2, sd = 0.1)
   candidate[c(2, 4)] <- candidate[c(2, 4)] + rnorm(2, sd = 0.1)
   ratio <- posterior(candidate) / posterior(pars)
   if (runif(1) < ratio) {
      pars <- candidate
   }s
   results[i, ] <- pars
}

Res_MCMC <- results[501:N_iter, ]
par(mfrow = c(2, 2))
plot(Res_MCMC[1:(i-1), 1])
plot(Res_MCMC[1:(i-1), 2])
plot(Res_MCMC[1:(i-1), 3])
plot(Res_MCMC[1:(i-1), 4])

plot(results[1:(i-1), 1])
plot(results[1:(i-1), 2])
plot(results[1:(i-1), 3])
plot(results[1:(i-1), 4])

