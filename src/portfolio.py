"""Module 5 - More than two assets: the best mix of stocks, gold, bonds and the bank.

Put weights w = (w_1, ..., w_n) in n risky assets (GBMs with drifts mu_i, volatilities
sigma_i, correlations rho_ij) and 1 - sum(w) in the bank (rate r). Rebalancing
continuously, as in Module 4, log-wealth after T years is normal with
    growth  g(w) = r + w.(mu - r) - w'Cw / 2,      spread  sd(w) = sqrt(w'Cw),
where C_ij = rho_ij sigma_i sigma_j is the covariance matrix.

Key idea: a fixed mix d of the risky assets (d_i >= 0, sum d_i = 1) behaves like ONE
fund with
    mu_d = r + d.(mu - r),      sigma_d = sqrt(d'Cd),
because w = f d gives exactly g(w) = r + f (mu_d - r) - f^2 sigma_d^2 / 2. So putting a
fraction f of the money in that fund and 1 - f in the bank is Module 4, and the
fund's best f is min(f*, f_VaR, 1), with f_VaR found by the IVT + bisection.
Trying every fund on a grid (steps of 1%) and keeping the fastest-growing one gives
the best mix without borrowing (f <= 1) and without short selling (d_i >= 0).

The problem is convex: g is concave and sd(w) is a norm, so g(w) - z sd(w) is concave and the
mixes with VaR <= L form a convex set. So there is exactly one best mix. A grid twice as fine
contains the coarser one, so its best growth can only go up, and as the grid step goes to 0 it
converges to the true maximum: Module 5's own limit (best_mix_on_grid, refine_best_mix).
"""
import math

import numpy as np

from src.bisection import bisection
from src.normal import norm_cdf, norm_ppf
from src.risk import growth_rate, kelly_fraction, var_root_formula


def covariance(sigma, rho):
    """C_ij = rho_ij sigma_i sigma_j."""
    sigma = np.asarray(sigma, dtype=float)
    return np.outer(sigma, sigma) * np.asarray(rho, dtype=float)


def portfolio_growth(w, mu, sigma, rho, r):
    """Long-run growth g(w) = r + w.(mu - r) - w'Cw / 2 and spread sd(w) = sqrt(w'Cw)."""
    w, mu = np.asarray(w, dtype=float), np.asarray(mu, dtype=float)
    v = float(w @ covariance(sigma, rho) @ w)
    return r + float(w @ (mu - r)) - 0.5 * v, math.sqrt(v)


def portfolio_var(w, mu, sigma, rho, r, T, alpha=0.95):
    """Multi-asset Value-at-Risk. With one asset it equals value_at_risk (Module 4)."""
    g, sd = portfolio_growth(w, mu, sigma, rho, r)
    return 1.0 - math.exp(g * T - norm_ppf(alpha) * sd * math.sqrt(T))


def portfolio_es(w, mu, sigma, rho, r, T, alpha=0.95):
    """Multi-asset Expected Shortfall. With one asset it equals risk.expected_shortfall (Module 4)."""
    g, sd = portfolio_growth(w, mu, sigma, rho, r)
    s = sd * math.sqrt(T)
    return 1.0 - math.exp(g * T + 0.5 * s**2) * float(norm_cdf(-norm_ppf(alpha) - s)) / (1.0 - alpha)


def fund(d, mu, sigma, rho, r):
    """The mix d of the risky assets seen as one fund: (mu_d, sigma_d)."""
    d, mu = np.asarray(d, dtype=float), np.asarray(mu, dtype=float)
    return r + float(d @ (mu - r)), math.sqrt(float(d @ covariance(sigma, rho) @ d))


def fund_fraction(d, L, mu, sigma, rho, r, T, alpha=0.95):
    """Module 4 for the fund d: f = min(f*, f_VaR, 1). Returns (f, growth, history)."""
    return _fund_fraction(d, L, mu, covariance(sigma, rho), r, T, norm_ppf(alpha))


def _fund_fraction(d, L, mu, C, r, T, z):
    # Same as risk.optimal_fraction with f_max = 1, but C and z_alpha are computed once,
    # because best_mix_under_var runs this for thousands of funds.
    mu_d = r + float(d @ (mu - r))
    s_d = math.sqrt(float(d @ C @ d))
    if mu_d <= r or s_d == 0.0:                      # never beats the bank
        return 0.0, r, np.empty((0, 4))
    h = lambda f: 1.0 - math.exp(growth_rate(f, mu_d, s_d, r) * T - z * f * s_d * math.sqrt(T)) - L
    if h(0.0) > 0:
        raise ValueError("risk limit L is below the VaR of the risk-free portfolio")
    if h(1.0) <= 0:
        f_var, history = 1.0, np.empty((0, 4))
    else:
        f_var, history = bisection(h, 0.0, 1.0)
    f = min(kelly_fraction(mu_d, s_d, r), f_var, 1.0)
    return f, growth_rate(f, mu_d, s_d, r), history


def mix_grid(n_assets, steps=100):
    """Every mix d of n_assets with d_i in {0, 1/steps, ..., 1} and sum 1, as rows of a matrix.

    Rows come in a fixed order, first entry slowest (the same order as the MATLAB and JS code).
    """
    rows, left = np.zeros((1, 0), dtype=int), np.array([steps])
    for _ in range(n_assets - 1):                   # the next asset gets 0, 1, ..., whatever is left
        new = np.concatenate([np.arange(k + 1) for k in left])
        rows = np.column_stack([np.repeat(rows, left + 1, axis=0), new])
        left = np.repeat(left, left + 1) - new
    return np.column_stack([rows, left]) / steps


def best_mix_under_var(L, mu, sigma, rho, r, T, alpha=0.95, steps=100):
    """Fastest-growing mix with VaR <= L, no borrowing, no short selling.

    Returns (w, f, d): the weights in the risky assets (the bank gets 1 - sum(w)),
    the fraction f in the winning fund and the fund d itself (w = f d).
    """
    mu, C, z = np.asarray(mu, dtype=float), covariance(sigma, rho), norm_ppf(alpha)
    best_g, best_f, best_d = r, 0.0, None
    for d in mix_grid(len(mu), steps):
        f, g, _ = _fund_fraction(d, L, mu, C, r, T, z)
        if g > best_g + 1e-15:
            best_g, best_f, best_d = g, f, d
    if best_d is None:
        return np.zeros(len(mu)), 0.0, None
    return best_f * best_d, best_f, best_d


# --- Module 5's limit: a finer and finer grid ------------------------------------
# For the bell curve, VaR(f) = L is a quadratic equation (risk.var_root_formula), so f_VaR of many
# funds can be computed at once. That makes grids with millions of mixes fast; the bisection
# above gives the same answer on the 1% grid (tested).

def fund_fractions(D, L, mu, sigma, rho, r, T, alpha=0.95):
    """Module 4 for many funds at once (one per row of D): f = min(f*, f_VaR, 1). Returns (f, growth)."""
    D, mu = np.asarray(D, dtype=float), np.asarray(mu, dtype=float)
    a = D @ (mu - r)                                                        # mu_d - r
    s = np.sqrt(np.einsum("ij,jk,ik->i", D, covariance(sigma, rho), D))   # sigma_d
    f = np.zeros(len(D))
    ok = (a > 0) & (s > 0)                                                  # the others never beat the bank
    if ok.any():
        f_var = var_root_formula(L, r + a[ok], s[ok], r, T, alpha)
        f[ok] = np.minimum(np.minimum(a[ok] / s[ok]**2, f_var), 1.0)
    return f, r + f * a - 0.5 * f**2 * s**2


def best_mix_on_grid(L, mu, sigma, rho, r, T, alpha=0.95, steps=100):
    """The fastest-growing fund on the grid with step 1/steps. Returns (w, f, d, growth).

    With steps = 100 it is the mix of best_mix_under_var. If one grid's steps divide the other's, the
    finer grid contains the coarser one, so its best growth is the same or higher.
    """
    D = mix_grid(len(mu), steps)
    f, g = fund_fractions(D, L, mu, sigma, rho, r, T, alpha)
    k = int(np.argmax(g))
    if g[k] <= r + 1e-15:                                                   # nothing beats the bank
        return np.zeros(len(mu)), 0.0, None, r
    return f[k] * D[k], f[k], D[k], float(g[k])


def refine_best_mix(L, mu, sigma, rho, r, T, alpha=0.95, steps=100, levels=40, points=100):
    """Zoom in on the best mix: let the grid step go to 0 around the best fund. Returns (w, f, d, growth).

    Start from the best fund on the 1/steps grid. Each level tries funds on a small grid around the
    current best (+-10 coarse steps at first, `points` intervals per direction) and then halves that
    window. The current best is always tried again, so the growth never goes down.
    """
    mu = np.asarray(mu, dtype=float)
    w, f, d, g = best_mix_on_grid(L, mu, sigma, rho, r, T, alpha, steps)
    if d is None:
        return w, f, d, g
    offsets = np.linspace(-10.0 / steps, 10.0 / steps, points + 1)        # contains 0 (points is even)
    for _ in range(levels):
        shift = np.stack([x.ravel() for x in np.meshgrid(*[offsets] * (len(mu) - 1), indexing="ij")], axis=1)
        D = np.column_stack([d[:-1] + shift, d[-1] - shift.sum(axis=1)])   # moves that keep sum(d) = 1
        D = np.clip(D[(D > -1e-12).all(axis=1)], 0.0, None)
        fs, gs = fund_fractions(D, L, mu, sigma, rho, r, T, alpha)
        k = int(np.argmax(gs))
        if gs[k] > g:
            f, d, g = fs[k], D[k], float(gs[k])
        offsets = offsets / 2
    return f * d, f, d, g


def near_optimal_ranges(L, mu, sigma, rho, r, T, alpha=0.95, steps=100, tol=1e-4):
    """Lowest and highest weight of each asset over the grid mixes growing within tol of the best.

    The top of the growth surface is flat, so many mixes are almost as good. Returns (low, high, count).
    """
    D = mix_grid(len(mu), steps)
    f, g = fund_fractions(D, L, mu, sigma, rho, r, T, alpha)
    W = f[:, None] * D
    near = g >= g.max() - tol
    return W[near].min(axis=0), W[near].max(axis=0), int(near.sum())


def simulate_wealth_mix(w, mu, sigma, rho, r, T, n_steps, n_paths, rng, W0=1.0):
    """Monte Carlo wealth paths, rebalanced to the weights w every dt = T / n_steps.

    Each step the risky assets get correlated GBM returns (Cholesky factor of C) and
    the rest earns e^(r dt). As dt -> 0 the VaR approaches portfolio_var.
    Returns an (n_paths, n_steps + 1) array.
    """
    w, mu, sigma = np.asarray(w, dtype=float), np.asarray(mu, dtype=float), np.asarray(sigma, dtype=float)
    chol, dt = np.linalg.cholesky(covariance(sigma, rho)), T / n_steps
    W = np.empty((n_paths, n_steps + 1))
    W[:, 0] = W0
    for k in range(n_steps):
        z = rng.standard_normal((n_paths, len(mu))) @ chol.T
        gross = np.exp((mu - 0.5 * sigma**2) * dt + z * math.sqrt(dt))
        W[:, k + 1] = W[:, k] * (gross @ w + (1.0 - w.sum()) * math.exp(r * dt))
    return W
