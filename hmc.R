# Author: Marvin Ahlborn
# 
# Implementation of the Hamiltonian Monte Carlo Algorithm based
# on its Wikipedia article and "Bayesian Statistics" by Ben Lambert.

library(MASS)

grad <- function(f, x0, eps=1e-4)
# Computes gradient of function f numerically via
# finite difference at point x0.
{
    n <- length(x0)
    f0 <- f(x0)
    del_f <- vector(mode="numeric", length=n)

    for (i in 1:n)
    {
        xeps <- x0
        xeps[i] <- xeps[i] + eps
        del_f[i] <- (f(xeps) - f(x0)) / eps
    }
    return(del_f)
}

hmc <- function(prob, x0, M, dt, L, N)
# Hamiltonian Monte-Carlo sampling algorithm.
# p is the unnormalized probability density to sample from.
# x0 is the starting point for sampling.
# M is the symmetric and positive definite mass matrix.
# dt is Leapfrog integration time interval.
# L is the number of integration steps per sample.
# N is the number of samples to sample.
{
    len <- length(x0)
    M_inv <- solve(M)
    U <- function(x) min(1e6, -log(prob(x)))  # Potential (capped for numerical stability reasons)
    H <- function(x, p) U(x) + 0.5 * p %*% (M_inv %*% p)  # Hamiltonian
    samples <- matrix(0, nrow=N, ncol=len, byrow=TRUE)
    samples[1,] <- x0

    x <- x0
    # Sampling
    for (n in 2:N)
    {
        x_init <- x
        p_init <- mvrnorm(1, vector(mode="numeric", length=len), M)  # Random initial momentum
        p <- p_init

        # Leapfrog integration
        for (t in 1:L)
        {
            grad_U <- grad(U, x)
            p_temp <- p - dt/2 * grad_U
            x <- x + dt * (M_inv %*% p_temp)

            grad_U <- grad(U, x)
            p <- p_temp - dt/2 * grad_U
        }

        # Accept-Reject step
        H_next <- H(x, p)
        H_last <- H(x_init, p_init)
        a <- min(1, exp(-H_next) / exp(-H_last))

        if (!rbinom(1, size=1, prob=a)) {
            x <- x_init  # Reject step
        }

        samples[n,] <- x
    }
    return(samples)
}