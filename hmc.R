library(MASS)

grad <- function(f, x0, eps=1e-5)
# Computes gradient of function f numerically via
# finite difference at point x0.
{
    n <- length(x0)
    f0 <- f(x0)
    grad <- vector(mode="numeric", length=n)

    for (i in 1:n)
    {
        xeps <- x0
        xeps[i] <- xeps[i] + eps
        grad[i] <- (f(xeps) - f(x0)) / eps
    }
    return(grad)
}

hmc <- function(p, x0, M, dt, L, N)
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
    U <- function(x) -log(p(x))  # Potential
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
        for (i in 1:L)
        {
            # TODO
            grad_U <- grad(U, x)
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