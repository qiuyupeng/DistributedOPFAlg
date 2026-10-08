"""Generic maximum likelihood estimation for any parameterized pdf.

The estimator only needs a callable ``pdf(x, *theta)`` (or ``logpdf``) that
evaluates the density of the data ``x`` under the parameters ``theta``. Nothing
about the distribution's family is assumed: the negative log-likelihood is
minimized numerically and standard errors come from the observed Fisher
information (the numerical Hessian of the negative log-likelihood).
"""

import numpy as np
from scipy.optimize import minimize


class MLEResult(object):
    def __init__(self, params, std_errors, cov, log_likelihood, success, message, n_iter):
        self.params = params
        self.std_errors = std_errors
        self.cov = cov
        self.log_likelihood = log_likelihood
        self.success = success
        self.message = message
        self.n_iter = n_iter

    def __repr__(self):
        return ("MLEResult(params=%s, std_errors=%s, log_likelihood=%.6g, success=%s)"
                % (np.array2string(self.params, precision=6),
                   np.array2string(self.std_errors, precision=6),
                   self.log_likelihood, self.success))


def _numerical_hessian(f, theta, rel_step=1e-4):
    """Central finite-difference Hessian of a scalar function f at theta."""
    theta = np.asarray(theta, dtype=float)
    k = theta.size
    h = rel_step * np.maximum(np.abs(theta), 1.0)
    hess = np.empty((k, k))
    f0 = f(theta)
    for i in range(k):
        ei = np.zeros(k)
        ei[i] = h[i]
        hess[i, i] = (f(theta + ei) - 2.0 * f0 + f(theta - ei)) / h[i] ** 2
        for j in range(i + 1, k):
            ej = np.zeros(k)
            ej[j] = h[j]
            hess[i, j] = hess[j, i] = (f(theta + ei + ej) - f(theta + ei - ej)
                                       - f(theta - ei + ej) + f(theta - ei - ej)) / (4.0 * h[i] * h[j])
    return hess


class MaxLikelihoodEstimator(object):
    """Maximum likelihood estimator for an arbitrary parameterized density.

    Parameters
    ----------
    pdf : callable, optional
        ``pdf(x, *theta)`` returning the density of each data point (vectorized
        over ``x``). Exactly one of ``pdf`` / ``logpdf`` must be given.
    logpdf : callable, optional
        ``logpdf(x, *theta)``; preferred when available because it is more
        numerically stable than ``log(pdf)``.
    bounds : sequence of (low, high), optional
        Box constraints per parameter, ``None`` for unbounded, e.g.
        ``[(None, None), (1e-9, None)]`` for a location and a positive scale.
    """

    def __init__(self, pdf=None, logpdf=None, bounds=None):
        if (pdf is None) == (logpdf is None):
            raise ValueError("Provide exactly one of pdf or logpdf.")
        if logpdf is None:
            def logpdf(x, *theta):
                with np.errstate(divide="ignore", invalid="ignore"):
                    return np.log(pdf(x, *theta))
        self.logpdf = logpdf
        self.bounds = bounds

    def neg_log_likelihood(self, theta, data):
        with np.errstate(all="ignore"):
            ll = np.sum(self.logpdf(data, *theta))
        # Parameters outside the support (or producing zero density) are
        # infeasible; a large finite value keeps the optimizer well-behaved.
        if not np.isfinite(ll):
            return 1e300
        return -ll

    def fit(self, data, theta0, method=None, n_restarts=0, restart_scale=1.0, seed=None,
            options=None):
        """Estimate the parameters by minimizing the negative log-likelihood.

        Parameters
        ----------
        data : array_like
            Observed samples (shape ``(n,)`` or ``(n, d)`` for multivariate pdfs).
        theta0 : array_like
            Initial guess for the parameters.
        method : str, optional
            ``scipy.optimize.minimize`` method. Defaults to L-BFGS-B when bounds
            are given, otherwise Nelder-Mead (derivative-free, robust).
        n_restarts : int
            Extra random restarts around ``theta0`` to guard against local optima.
        restart_scale : float
            Std of the Gaussian perturbation used for restarts.
        """
        data = np.asarray(data, dtype=float)
        theta0 = np.atleast_1d(np.asarray(theta0, dtype=float))
        if method is None:
            method = "L-BFGS-B" if self.bounds is not None else "Nelder-Mead"
        if options is None and method == "Nelder-Mead":
            options = {"xatol": 1e-10, "fatol": 1e-12, "maxiter": 20000 * theta0.size}
        rng = np.random.default_rng(seed)

        nll = lambda theta: self.neg_log_likelihood(theta, data)

        starts = [theta0] + [theta0 + restart_scale * rng.standard_normal(theta0.size)
                             for _ in range(n_restarts)]
        best = None
        for start in starts:
            start = self._clip(start)
            res = minimize(nll, start, method=method, bounds=self.bounds, options=options)
            if best is None or res.fun < best.fun:
                best = res

        theta_hat = best.x
        cov, std_errors = self._covariance(nll, theta_hat)
        return MLEResult(params=theta_hat, std_errors=std_errors, cov=cov,
                         log_likelihood=-best.fun, success=bool(best.success),
                         message=best.message, n_iter=best.nit)

    def _clip(self, theta):
        if self.bounds is None:
            return theta
        lo = np.array([-np.inf if b[0] is None else b[0] for b in self.bounds])
        hi = np.array([np.inf if b[1] is None else b[1] for b in self.bounds])
        return np.clip(theta, lo, hi)

    @staticmethod
    def _covariance(nll, theta_hat):
        """Inverse observed Fisher information = asymptotic covariance of theta_hat."""
        try:
            hess = _numerical_hessian(nll, theta_hat)
            cov = np.linalg.inv(hess)
            std_errors = np.sqrt(np.diag(cov))
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            k = theta_hat.size
            cov = np.full((k, k), np.nan)
            std_errors = np.full(k, np.nan)
        return cov, std_errors


def gaussian_pdf(x, mu, sigma):
    if sigma <= 0:
        return np.zeros_like(x)
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))


if __name__ == "__main__":
    rng = np.random.default_rng(0)

    # Gaussian example: only the pdf is supplied, the estimator knows nothing else.
    true_mu, true_sigma = 3.0, 1.5
    samples = rng.normal(true_mu, true_sigma, size=5000)
    mle = MaxLikelihoodEstimator(pdf=gaussian_pdf, bounds=[(None, None), (1e-9, None)])
    result = mle.fit(samples, theta0=[0.0, 1.0])
    print("Gaussian  true (mu, sigma) = (%.3f, %.3f)" % (true_mu, true_sigma))
    print("          %s" % result)
    print("          closed form      = (%.6f, %.6f)" % (samples.mean(), samples.std()))

    # Same procedure on a different family: Gamma(shape k, scale theta).
    from scipy.special import gammaln

    def gamma_logpdf(x, k, theta):
        return (k - 1) * np.log(x) - x / theta - gammaln(k) - k * np.log(theta)

    true_k, true_theta = 2.5, 0.8
    samples = rng.gamma(true_k, true_theta, size=5000)
    mle = MaxLikelihoodEstimator(logpdf=gamma_logpdf, bounds=[(1e-9, None), (1e-9, None)])
    result = mle.fit(samples, theta0=[1.0, 1.0])
    print("Gamma     true (k, theta)  = (%.3f, %.3f)" % (true_k, true_theta))
    print("          %s" % result)
