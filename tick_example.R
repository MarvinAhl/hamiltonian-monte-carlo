source("hmc.R")

prob <- function(param) dbinom(6, size=100, prob=param%%1) * dbeta(param, 1, 1)
param <- c(runif(1, min=0, max=1))  # Just some random initial value
M <- diag(1)  # Mass Matrix, when in doubt, just set to diag
dt <- 0.001  # Delta t of one integration time step
L <- 20  # Integration time steps
N <- 1000  # Number of samples to generate

samples <- hmc(prob, param, M, dt, L, N)
samples <- samples %% 1  # Dont forget periodic boundaries

# Multiple plots
par(mfrow=c(2, 1))
x <- seq(0, 1, 0.005)
plot(x, dbeta(x, 7, 95), main="Analytic Posterior", xlab="x", ylab="beta", type="l", col="Black", lwd=2)
hist(samples, main="HMC posterior samples", breaks=seq(0, 1, 0.02), xlim=c(0, 1))