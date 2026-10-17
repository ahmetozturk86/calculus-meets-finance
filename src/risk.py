"""Module 4 - Risk optimization with limits and continuity.

Put a fraction f of wealth in a stock (GBM with drift mu, volatility sigma) and
1 - f in a risk-free asset that compounds continuously at rate r (Module 1).
Rebalancing continuously - the dt -> 0 limit of Module 3 - wealth after T years is
    W_T = W_0 exp( g(f) T + f sigma sqrt(T) Z ),     Z ~ N(0, 1),
    g(f) = r + f (mu - r) - f^2 sigma^2 / 2          (long-run growth rate).
All rates are continuously compounded (Module 1): mu = 10% is an average yearly gain of e^0.1 - 1 = 10.5%.

1. Growth: g is continuous and concave, g'(f) = (mu - r) - f sigma^2, so the best
   growth is at f* = (mu - r) / sigma^2 (the Kelly / Merton fraction). We also find
   it numerically: g'(f) as the limit of difference quotients, then bisection on g'.
2. Risk: the Value-at-Risk VaR(f) is continuous in f. For a risk limit L,
   h(f) = VaR(f) - L has h(0) < 0 < h(f_max), so by the Intermediate Value
   Theorem some f_VaR gives VaR(f_VaR) = L exactly. Bisection finds it. It is the
   only root for any horizon T: VaR(f) = L is a quadratic equation in f whose two
   roots have opposite signs (see var_root_formula).
3. Decision: f_opt = min(f*, f_VaR, f_max).
"""
import math

import numpy as np

from src.bisection import bisection
from src.normal import norm_cdf, norm_ppf


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
    not binding and f_max is returned. The root is unique (var_root_formula), so every
    f below it keeps the limit, even for long horizons where VaR(f) first falls.
    """
    h = lambda f: value_at_risk(f, mu, sigma, r, T, alpha) - L
    if h(0.0) > 0:
        raise ValueError("risk limit L is below the VaR of the risk-free portfolio")
    if h(f_max) <= 0:
        return f_max, np.empty((0, 4))
    return bisection(h, 0.0, f_max)


def var_root_formula(L, mu, sigma, r, T, alpha=0.95):
    """The root of VaR(f) = L from the quadratic formula, a check on the bisection.

    VaR(f) = L  <=>  (sigma^2 T / 2) f^2 + (z sigma sqrt(T) - (mu - r) T) f + (ln(1 - L) - r T) = 0.
    The last coefficient is negative exactly when h(0) < 0, so the product of the two roots is
    negative: one root is negative and one positive, for every horizon T. The positive one is
    the only f >= 0 with VaR(f) = L. Works on arrays of mu and sigma (one fund per entry).
    """
    if 1.0 - math.exp(r * T) - L > 0:
        raise ValueError("risk limit L is below the VaR of the risk-free portfolio")
    sigma = np.asarray(sigma, dtype=float)
    a = 0.5 * sigma**2 * T
    b = norm_ppf(alpha) * sigma * math.sqrt(T) - (np.asarray(mu, dtype=float) - r) * T
    if L >= 1.0:                                     # no limit: a loss of 100% never has 5% chance
        root = np.full(np.broadcast(a, b).shape, np.inf)
    else:
        c = math.log(1.0 - L) - r * T
        disc = np.sqrt(b**2 - 4.0 * a * c)
        with np.errstate(divide="ignore", invalid="ignore"):   # two equal forms; each avoids cancellation
            root = np.where(b >= 0, -2.0 * c / (b + disc), (disc - b) / (2.0 * a))
    return float(root) if root.ndim == 0 else root


# --- Fat tails: the bell curve plus rare crashes (Merton jump-diffusion) -------
# A crash multiplies the stock by (1 + J) (e.g. J = -0.25) and happens on average
# lam times a year. The drift is compensated (mu - lam J), so the stock's average
# return is still mu: only the shape changes, with a fatter left tail. With
# N ~ Poisson(lam T) crashes and continuous rebalancing,
#     log(W_T / W_0) = (g(f) - lam J f) T + N ln(1 + f J) + f sigma sqrt(T) Z.

def poisson_pmf(k, m):
    """P(N = k) for N ~ Poisson(m)."""
    if m == 0:
        return 1.0 if k == 0 else 0.0
    return math.exp(-m + k * math.log(m) - math.lgamma(k + 1))


def binomial_pmf(k, n, p):
    """P(B = k) for B ~ Binomial(n, p). As n -> inf with p = lam / n it tends to poisson_pmf(k, lam)."""
    if p == 0:
        return 1.0 if k == 0 else 0.0
    if k > n:
        return 0.0
    return math.exp(math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
                    + k * math.log(p) + (n - k) * math.log1p(-p))


def jump_log_wealth_cdf(x, f, mu, sigma, r, T, lam, J, k_max=30):
    """P(log(W_T / W_0) <= x): a mixture of bell curves, one for each number of crashes k."""
    m = (growth_rate(f, mu, sigma, r) - lam * J * f) * T
    sd = abs(f) * sigma * math.sqrt(T)
    jump = math.log(1.0 + f * J)
    k = np.arange(k_max + 1)
    weights = np.array([poisson_pmf(i, lam * T) for i in k])
    if sd == 0.0:                                   # no bell-curve noise (f = 0): steps at m + k jump
        return float(np.sum(weights * (x >= m + k * jump)))
    return float(np.sum(weights * norm_cdf((x - m - k * jump) / sd)))


def value_at_risk_jump(f, mu, sigma, r, T, lam, J, alpha=0.95, k_max=30):
    """VaR at level alpha when rare crashes can happen. With lam = 0 it equals value_at_risk.

    The mixture CDF has no inverse formula, but it is continuous and increasing, so by
    the IVT it crosses 1 - alpha exactly once, and bisection finds that quantile.
    """
    if f == 0:
        return 1.0 - math.exp(r * T)
    if 1.0 + f * J <= 0:
        raise ValueError("a crash would wipe out the portfolio: need 1 + f J > 0")
    m = (growth_rate(f, mu, sigma, r) - lam * J * f) * T
    sd = abs(f) * sigma * math.sqrt(T)
    lo = m + k_max * min(0.0, math.log(1.0 + f * J)) - 12.0 * sd
    hi = m + max(0.0, k_max * math.log(1.0 + f * J)) + 12.0 * sd
    x, _ = bisection(lambda x: jump_log_wealth_cdf(x, f, mu, sigma, r, T, lam, J, k_max) - (1.0 - alpha),
                     lo, hi, tol=1e-12)
    return 1.0 - math.exp(x)


def max_fraction_under_var_jump(L, mu, sigma, r, T, lam, J, alpha=0.95, f_max=2.0):
    """Largest f in [0, f_max] with VaR_jump(f) <= L: bisection on h(f) = VaR_jump(f) - L (IVT)."""
    h = lambda f: value_at_risk_jump(f, mu, sigma, r, T, lam, J, alpha) - L
    if h(0.0) > 0:
        raise ValueError("risk limit L is below the VaR of the risk-free portfolio")
    if h(f_max) <= 0:
        return f_max, np.empty((0, 4))
    return bisection(h, 0.0, f_max)


def growth_rate_jump(f, mu, sigma, r, lam, J):
    """Long-run growth with crashes, E[log(W_T / W_0)] / T = g(f) - lam J f + lam ln(1 + f J).

    The average return is the same as without crashes, but ln(1 + x) < x for x != 0, so crashes
    lower the growth: a big loss hurts more than the extra drift of the good years helps.
    """
    f = np.asarray(f, dtype=float)
    return growth_rate(f, mu, sigma, r) - lam * J * f + lam * np.log1p(f * J)


def kelly_fraction_jump(mu, sigma, r, lam, J, f_hi=10.0):
    """Kelly fraction with crashes: bisection on the slope of growth_rate_jump (IVT). Returns (f, history).

    The slope (mu - r) - f sigma^2 - lam J + lam J / (1 + f J) is continuous, equals mu - r > 0 at
    f = 0, and goes to -infinity where one crash would wipe out the portfolio (1 + f J -> 0).
    """
    if J < 0:
        f_hi = min(f_hi, 0.999999 / -J)
    slope = lambda f: (mu - r) - f * sigma**2 - lam * J + lam * J / (1.0 + f * J)
    return bisection(slope, 0.0, f_hi, tol=1e-10)


def optimal_fraction_jump(mu, sigma, r, T, L, lam, J, alpha=0.95, f_max=2.0):
    """Best growth when crashes can happen: min(f*_crash, f_crash, f_max), at least 0."""
    if mu <= r:
        return 0.0
    f_star, _ = kelly_fraction_jump(mu, sigma, r, lam, J)
    f_var, _ = max_fraction_under_var_jump(L, mu, sigma, r, T, lam, J, alpha, f_max)
    return max(0.0, min(f_star, f_var, f_max))


# --- Expected Shortfall: how bad the bad years are on average ---------------------
# VaR says where the worst 1 - alpha of years start; ES is the average loss in them.

def expected_shortfall(f, mu, sigma, r, T, alpha=0.95):
    """Expected Shortfall: the average loss in the worst 1 - alpha of cases (never below the VaR).

    log(W_T / W_0) ~ N(m, s^2) with m = g(f) T and s = |f| sigma sqrt(T). Completing the square (as
    for E[e^X]) gives E[W_T / W_0 | worst 1 - alpha] = e^(m + s^2/2) Phi(-z_alpha - s) / (1 - alpha).
    """
    z = norm_ppf(alpha)
    f = np.asarray(f, dtype=float)
    m = growth_rate(f, mu, sigma, r) * T
    s = np.abs(f) * sigma * math.sqrt(T)
    return 1.0 - np.exp(m + 0.5 * s**2) * norm_cdf(-z - s) / (1.0 - alpha)


def expected_shortfall_jump(f, mu, sigma, r, T, lam, J, alpha=0.95, k_max=30):
    """Expected Shortfall when rare crashes can happen. With lam = 0 it equals expected_shortfall.

    Below the VaR point x_q = ln(1 - VaR_crash) (found by bisection), the bell curve of k crashes adds
    e^(m_k + s^2/2) Phi((x_q - m_k - s^2) / s), with m_k = (g(f) - lam J f) T + k ln(1 + f J).
    """
    if f == 0:
        return 1.0 - math.exp(r * T)
    x_q = math.log(1.0 - value_at_risk_jump(f, mu, sigma, r, T, lam, J, alpha, k_max))
    s = abs(f) * sigma * math.sqrt(T)
    k = np.arange(k_max + 1)
    weights = np.array([poisson_pmf(i, lam * T) for i in k])
    m_k = (growth_rate(f, mu, sigma, r) - lam * J * f) * T + k * math.log(1.0 + f * J)
    tail = np.sum(weights * np.exp(m_k + 0.5 * s**2) * norm_cdf((x_q - m_k - s**2) / s))
    return float(1.0 - tail / (1.0 - alpha))


def max_fraction_under_es(L, mu, sigma, r, T, alpha=0.95, f_max=2.0):
    """Largest f in [0, f_max] with ES(f) <= L: the same IVT + bisection as for VaR. Returns (f, history).

    ES is continuous in f, and ln(1 - ES(f)) = (r + f (mu - r)) T + ln Phi(-z - f sigma sqrt(T)) - ln(1 - alpha)
    is concave (Phi is log-concave), so the root is unique, as for VaR.
    """
    h = lambda f: float(expected_shortfall(f, mu, sigma, r, T, alpha)) - L
    if h(0.0) > 0:
        raise ValueError("risk limit L is below the ES of the risk-free portfolio")
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


def simulate_wealth_jump(f, mu, sigma, r, T, lam, J, n_steps, n_paths, rng, W0=1.0):
    """Like simulate_wealth, but each step a crash happens with probability lam * dt (a coin flip).

    As dt -> 0 the number of crashes in a year, Binomial(n_steps, lam dt), tends to Poisson(lam T).
    """
    dt = T / n_steps
    rf_return = math.exp(r * dt) - 1.0
    W = np.empty((n_paths, n_steps + 1))
    W[:, 0] = W0
    for k in range(n_steps):
        Z = rng.standard_normal(n_paths)
        crash = rng.random(n_paths) < lam * dt
        stock_return = np.exp((mu - lam * J - 0.5 * sigma**2) * dt + sigma * math.sqrt(dt) * Z) * (1.0 + J * crash) - 1.0
        W[:, k + 1] = W[:, k] * (1.0 + f * stock_return + (1.0 - f) * rf_return)
    return W


def empirical_var(final_wealth, W0=1.0, alpha=0.95):
    """VaR estimated from simulated final wealth."""
    return 1.0 - np.quantile(np.asarray(final_wealth) / W0, 1.0 - alpha)


def empirical_es(final_wealth, W0=1.0, alpha=0.95):
    """ES estimated from simulated final wealth: 1 - the average of the worst 1 - alpha share."""
    x = np.sort(np.asarray(final_wealth) / W0)
    k = max(1, int(math.floor((1.0 - alpha) * len(x) + 0.5)))
    return 1.0 - x[:k].mean()
