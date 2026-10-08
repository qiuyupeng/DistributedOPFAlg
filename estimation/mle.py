"""Maximum likelihood estimation for any parameterized pdf."""

import numpy as np
from scipy.optimize import minimize


def mle(pdf, data, theta0, bounds=None):
    """Return the theta maximizing sum(log pdf(x_i, *theta)) over the data."""
    data = np.asarray(data, dtype=float)

    def neg_log_likelihood(theta):
        return -np.sum(np.log(pdf(data, *theta)))

    return minimize(neg_log_likelihood, theta0, bounds=bounds).x


if __name__ == "__main__":
    def gaussian_pdf(x, mu, sigma):
        return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))

    samples = np.random.normal(3.0, 1.5, size=5000)
    mu, sigma = mle(gaussian_pdf, samples, theta0=[0.0, 1.0],
                    bounds=[(None, None), (1e-9, None)])
    print("mu = %.4f, sigma = %.4f" % (mu, sigma))
