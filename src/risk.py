"""Module 4 - Risk optimization with limits and continuity.

Put a fraction f of wealth in a stock (GBM with drift mu, volatility sigma) and
1 - f in a risk-free asset that compounds continuously at rate r (Module 1).
Rebalancing continuously - the dt -> 0 limit of Module 3 - wealth after T years is
    W_T = W_0 exp( g(f) T + f sigma sqrt(T) Z ),     Z ~ N(0, 1),
    g(f) = r + f (mu - r) - f^2 sigma^2 / 2          (long-run growth rate).

1. Growth: g is continuous and concave, g'(f) = (mu - r) - f sigma^2, so the best
   growth is at f* = (mu - r) / sigma^2 (the Kelly / Merton fraction). We also find
   it numerically: g'(f) as the limit of difference quotients, then bisection on g'.
2. Risk: the Value-at-Risk VaR(f) is continuous in f. For a risk limit L,
   h(f) = VaR(f) - L has h(0) < 0 < h(f_max), so by the Intermediate Value
   Theorem some f_VaR gives VaR(f_VaR) = L exactly. Bisection finds it.
3. Decision: f_opt = min(f*, f_VaR, f_max).
"""
import math

import numpy as np

from src.bisection import bisection
from src.normal import norm_ppf


def growth_rate(f, mu, sigma, r):
    """Long-run (log) growth rate g(f) of a portfolio with fraction f in the stock."""
    return r + f * (mu - r) - 0.5 * f**2 * sigma**2


def kelly_fraction(mu, sigma, r):
    """Analytic maximizer of g:  f* = (mu - r) / sigma^2."""
    return (mu - r) / sigma**2


def difference_quotient(fun, x, h):
    """(fun(x + h) - fun(x)) / h. Its limit as h -> 0 is the derivative fun'(x)."""
    return (fun(x + h) - fun(x)) / h


def derivative(fun, x, h=1e-6):
    """Numerical derivative (central difference quotient)."""
    return (fun(x + h) - fun(x - h)) / (2.0 * h)


def numerical_kelly_fraction(mu, sigma, r, f_lo=0.0, f_hi=10.0):
    """Maximize g without the formula: solve g'(f) = 0 by bisection (g' is continuous)."""
    g = lambda f: growth_rate(f, mu, sigma, r)
    return bisection(lambda f: derivative(g, f), f_lo, f_hi, tol=1e-10)


def value_at_risk(f, mu, sigma, r, T, alpha=0.95):
    """VaR at level alpha: the loss (fraction of starting wealth) exceeded with prob. 1 - alpha.

    The (1 - alpha) quantile of W_T / W_0 is exp(g(f) T - z_alpha |f| sigma sqrt(T)).
    A negative VaR means even the bad case is a gain.
    """
    z = norm_ppf(alpha)
    f = np.asarray(f, dtype=float)
    return 1.0 - np.exp(growth_rate(f, mu, sigma, r) * T - z * np.abs(f) * sigma * math.sqrt(T))


def max_fraction_under_var(L, mu, sigma, r, T, alpha=0.95, f_max=2.0):
    """Largest f in [0, f_max] with VaR(f) <= L, found by bisection on h(f) = VaR(f) - L.

    Returns (f_var, history). If even f_max respects the limit the constraint is
    not binding and f_max is returned.
    """
    h = lambda f: value_at_risk(f, mu, sigma, r, T, alpha) - L
    if h(0.0) > 0:
        raise ValueError("risk limit L is below the VaR of the risk-free portfolio")
    if h(f_max) <= 0:
        return f_max, np.empty((0, 4))
    return bisection(h, 0.0, f_max)


def optimal_fraction(mu, sigma, r, T, L, alpha=0.95, f_max=2.0):
    """Best growth subject to the VaR limit:  f_opt = min(f*, f_VaR, f_max), at least 0."""
    f_star = kelly_fraction(mu, sigma, r)
    f_var, _ = max_fraction_under_var(L, mu, sigma, r, T, alpha, f_max)
    return max(0.0, min(f_star, f_var, f_max))


def simulate_wealth(f, mu, sigma, r, T, n_steps, n_paths, rng, W0=1.0):
    """Monte Carlo wealth paths, rebalanced to fraction f every dt = T / n_steps.

    Returns an (n_paths, n_steps+1) array. As dt -> 0 this converges to the
    continuous-rebalancing formula above.
    """
    dt = T / n_steps
    rf_return = math.exp(r * dt) - 1.0
    W = np.empty((n_paths, n_steps + 1))
    W[:, 0] = W0
    for k in range(n_steps):
        Z = rng.standard_normal(n_paths)
        stock_return = np.exp((mu - 0.5 * sigma**2) * dt + sigma * math.sqrt(dt) * Z) - 1.0
        W[:, k + 1] = W[:, k] * (1.0 + f * stock_return + (1.0 - f) * rf_return)
    return W


def empirical_var(final_wealth, W0=1.0, alpha=0.95):
    """VaR estimated from simulated final wealth."""
    return 1.0 - np.quantile(np.asarray(final_wealth) / W0, 1.0 - alpha)
