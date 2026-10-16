"""Standard normal distribution helpers (no scipy, so the code ports to MATLAB 1:1)."""
import math

import numpy as np

from src.bisection import bisection

_erf = np.vectorize(math.erf)


def norm_pdf(x):
    """Standard normal density phi(x)."""
    x = np.asarray(x, dtype=float)
    return np.exp(-0.5 * x**2) / math.sqrt(2.0 * math.pi)


def norm_cdf(x):
    """Standard normal CDF  Phi(x) = (1 + erf(x / sqrt(2))) / 2."""
    x = np.asarray(x, dtype=float)
    return 0.5 * (1.0 + _erf(x / math.sqrt(2.0)))


def norm_ppf(q):
    """Inverse CDF: the z with Phi(z) = q.

    Phi is continuous and increasing, so Phi(z) - q changes sign on [-10, 10]
    and the IVT guarantees a solution, which bisection finds.
    """
    if not 0.0 < q < 1.0:
        raise ValueError("q must be in (0, 1)")
    z, _ = bisection(lambda x: norm_cdf(x) - q, -10.0, 10.0, tol=1e-13)
    return z
