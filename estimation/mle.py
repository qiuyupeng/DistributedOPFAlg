"""Maximum likelihood estimation for any parameterized pdf."""

import numpy as np
from scipy.optimize import minimize


def mle(pdf, data, theta0, bounds=None):
    """Return the theta maximizing sum(log pdf(x_i, *theta)) over the data."""
    data = np.asarray(data, dtype=float)

    def neg_log_likelihood(theta):
        return -np.sum(np.log(pdf(data, *theta)))

    return minimize(neg_log_likelihood, theta0, bounds=bounds).x


def sample(pdf, theta, n, low, high, grid_size=100000):
    """Draw n samples from pdf(x, *theta) on [low, high] by inverse transform sampling."""
    x = np.linspace(low, high, grid_size)
    cdf = np.cumsum(pdf(x, *theta))
    cdf /= cdf[-1]
    return np.interp(np.random.rand(n), cdf, x)


if __name__ == "__main__":
    def gaussian_pdf(x, mu, sigma):
        return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))

    samples = sample(gaussian_pdf, theta=(3.0, 1.5), n=5000, low=-20, high=20)
    mu, sigma = mle(gaussian_pdf, samples, theta0=[0.0, 1.0],
                    bounds=[(None, None), (1e-9, None)])
    print("mu = %.4f, sigma = %.4f" % (mu, sigma))
