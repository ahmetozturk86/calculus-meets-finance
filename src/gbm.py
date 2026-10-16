"""Module 3 - A coin-flip random walk converges to Geometric Brownian Motion.

Discrete model: every dt years the log-price moves up or down by sigma sqrt(dt),
    log S_{k+1} = log S_k +/- sigma sqrt(dt),     P(up) = p(dt),
with p chosen so the average drift matches GBM. After n = T/dt steps,
log(S_T / S_0) is a scaled binomial random variable.

As dt -> 0 (n -> inf) its distribution converges to the GBM limit
    log(S_T / S_0) ~ N((mu - sigma^2/2) T, sigma^2 T),
and the Kolmogorov-Smirnov distance between the two CDFs goes to 0 like 1/sqrt(n).
"""
import math

import numpy as np

from src.normal import norm_cdf


def up_probability(mu, sigma, dt):
    """P(up) that makes the walk's log-drift equal (mu - sigma^2/2) per year."""
    return 0.5 + 0.5 * (mu - 0.5 * sigma**2) / sigma * math.sqrt(dt)


def simulate_walk_paths(S0, mu, sigma, T, n_steps, n_paths, rng):
    """Simulate price paths of the random walk. Returns an (n_paths, n_steps+1) array."""
    dt = T / n_steps
    p = up_probability(mu, sigma, dt)
    steps = np.where(rng.random((n_paths, n_steps)) < p, 1.0, -1.0)
    log_paths = np.cumsum(sigma * math.sqrt(dt) * steps, axis=1)
    log_paths = np.hstack([np.zeros((n_paths, 1)), log_paths])
    return S0 * np.exp(log_paths)


def walk_log_return_distribution(mu, sigma, T, n_steps):
    """Exact distribution of log(S_T/S_0) for the walk: support points x and probabilities pmf."""
    n = n_steps
    dt = T / n
    p = up_probability(mu, sigma, dt)
    k = np.arange(n + 1)                      # number of up moves
    # log C(n, k) built from C(n, k) = C(n, k-1) * (n - k + 1) / k
    log_comb = np.concatenate(([0.0], np.cumsum(np.log((n - k[1:] + 1) / k[1:]))))
    pmf = np.exp(log_comb + k * math.log(p) + (n - k) * math.log1p(-p))
    x = sigma * math.sqrt(dt) * (2 * k - n)
    return x, pmf


def gbm_log_return_cdf(x, mu, sigma, T):
    """CDF of the GBM limit  N((mu - sigma^2/2) T, sigma^2 T)."""
    m = (mu - 0.5 * sigma**2) * T
    s = sigma * math.sqrt(T)
    return norm_cdf((np.asarray(x) - m) / s)


def ks_distance_to_gbm(mu, sigma, T, n_steps):
    """Exact Kolmogorov-Smirnov distance  sup_x |F_walk(x) - F_GBM(x)|.

    F_walk is a step function, so the sup is reached just before or at a jump.
    """
    x, pmf = walk_log_return_distribution(mu, sigma, T, n_steps)
    cdf_at = np.cumsum(pmf)                   # F_walk(x_k)
    cdf_before = cdf_at - pmf                 # F_walk(x_k^-)
    F = gbm_log_return_cdf(x, mu, sigma, T)
    return float(np.max(np.maximum(np.abs(cdf_at - F), np.abs(cdf_before - F))))
