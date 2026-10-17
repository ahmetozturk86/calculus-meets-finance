"""Checks every module against known values and the limits it is supposed to show."""
import math

import numpy as np
import pytest

from src.bisection import bisection
from src.compounding import compound, compounding_error, continuous, error_estimate
from src.gbm import ks_distance_to_gbm, walk_log_return_distribution
from src.normal import norm_cdf, norm_ppf
from src.option_pricing import binomial_price, black_scholes
from src.risk import (binomial_pmf, empirical_var, growth_rate, kelly_fraction, max_fraction_under_var,
                      max_fraction_under_var_jump, numerical_kelly_fraction, optimal_fraction, poisson_pmf,
                      simulate_wealth, simulate_wealth_jump, value_at_risk, value_at_risk_jump)

MU, SIGMA, R, T, ALPHA, L = 0.10, 0.20, 0.03, 1.0, 0.95, 0.15
LAM, J = 0.10, -0.25


# --- bisection / normal -------------------------------------------------------

def test_bisection_finds_sqrt2():
    root, history = bisection(lambda x: x**2 - 2, 0.0, 2.0)
    assert root == pytest.approx(math.sqrt(2), abs=1e-11)
    widths = history[:, 1] - history[:, 0]
    assert np.all(np.diff(widths) < 0)          # interval shrinks every step


def test_bisection_needs_sign_change():
    with pytest.raises(ValueError):
        bisection(lambda x: x**2 + 1, -1.0, 1.0)


def test_normal_helpers():
    assert float(norm_cdf(0.0)) == pytest.approx(0.5)
    assert norm_ppf(0.95) == pytest.approx(1.6448536269514722, abs=1e-10)


# --- module 1 -----------------------------------------------------------------

def test_continuous_compounding_value():
    assert continuous(1.0, 0.1, 1.0) == pytest.approx(1.105170918, abs=1e-9)


def test_compounding_converges_to_limit():
    ns = np.array([1, 12, 365, 10**6])
    values = compound(1.0, 0.1, 1.0, ns)
    assert np.all(np.diff(values) > 0)                      # more often -> more money
    assert compounding_error(1.0, 0.1, 1.0, 10**6) < 1e-7  # ... but bounded by the limit


def test_error_estimate_is_first_order():
    n = 10**5
    ratio = compounding_error(1.0, 0.1, 1.0, n) / error_estimate(1.0, 0.1, 1.0, n)
    assert ratio == pytest.approx(1.0, abs=1e-4)


# --- module 2 -----------------------------------------------------------------

def test_black_scholes_reference_and_parity():
    call = black_scholes(100, 100, 0.05, 0.2, 1.0, "call")
    put = black_scholes(100, 100, 0.05, 0.2, 1.0, "put")
    assert call == pytest.approx(10.4506, abs=1e-4)
    assert call - put == pytest.approx(100 - 100 * math.exp(-0.05), abs=1e-10)


def test_binomial_converges_to_black_scholes():
    bs = black_scholes(100, 100, 0.05, 0.2, 1.0)
    assert abs(binomial_price(100, 100, 0.05, 0.2, 1.0, 1000) - bs) < 0.01
    bs_put = black_scholes(100, 100, 0.05, 0.2, 1.0, "put")
    assert abs(binomial_price(100, 100, 0.05, 0.2, 1.0, 1000, "put") - bs_put) < 0.01


# --- module 3 -----------------------------------------------------------------

def test_walk_distribution_matches_gbm_moments():
    x, pmf = walk_log_return_distribution(MU, SIGMA, T, 500)
    assert pmf.sum() == pytest.approx(1.0, abs=1e-12)
    assert np.dot(x, pmf) == pytest.approx((MU - 0.5 * SIGMA**2) * T, abs=1e-12)


def test_ks_distance_goes_to_zero():
    ks = [ks_distance_to_gbm(MU, SIGMA, T, n) for n in [1, 100, 10_000]]
    assert ks[0] > ks[1] > ks[2]
    assert ks[2] < 0.01


# --- module 4 -----------------------------------------------------------------

def test_kelly_fraction_analytic_and_numerical():
    assert kelly_fraction(MU, SIGMA, R) == pytest.approx(1.75)
    f_num, _ = numerical_kelly_fraction(MU, SIGMA, R)
    assert f_num == pytest.approx(1.75, abs=1e-6)
    f = np.linspace(0, 3.5, 1001)
    assert f[np.argmax(growth_rate(f, MU, SIGMA, R))] == pytest.approx(1.75, abs=0.01)


def test_var_of_risk_free_portfolio_is_a_gain():
    assert float(value_at_risk(0.0, MU, SIGMA, R, T, ALPHA)) == pytest.approx(1 - math.exp(R * T))


def test_var_root_hits_the_limit():
    f_var, _ = max_fraction_under_var(L, MU, SIGMA, R, T, ALPHA)
    assert abs(float(value_at_risk(f_var, MU, SIGMA, R, T, ALPHA)) - L) < 1e-8
    assert optimal_fraction(MU, SIGMA, R, T, L, ALPHA) == pytest.approx(min(1.75, f_var))


def test_var_limit_not_binding_returns_f_max():
    f_var, _ = max_fraction_under_var(0.99, MU, SIGMA, R, T, ALPHA, f_max=2.0)
    assert f_var == 2.0


def test_monte_carlo_var_matches_formula():
    f = 0.7
    rng = np.random.default_rng(0)
    W = simulate_wealth(f, MU, SIGMA, R, T, 52, 50_000, rng)
    assert empirical_var(W[:, -1], 1.0, ALPHA) == pytest.approx(
        float(value_at_risk(f, MU, SIGMA, R, T, ALPHA)), abs=0.01)


# --- module 4b: fat tails (rare crashes) -----------------------------------------

def test_no_crashes_gives_the_normal_var():
    for f in [0.0, 0.3, 0.705, 1.5]:
        for a in [0.9, 0.95, 0.99]:
            assert value_at_risk_jump(f, MU, SIGMA, R, T, 0.0, J, a) == pytest.approx(
                float(value_at_risk(f, MU, SIGMA, R, T, a)), abs=1e-10)


def test_binomial_crash_count_converges_to_poisson():
    errors = [max(abs(binomial_pmf(k, n, LAM / n) - poisson_pmf(k, LAM)) for k in range(31))
              for n in [10, 1_000, 100_000]]
    assert errors[0] > errors[1] > errors[2] and errors[2] < 1e-6


def test_crash_root_hits_the_limit_and_is_smaller():
    f_var, _ = max_fraction_under_var(L, MU, SIGMA, R, T, ALPHA)
    f_crash, _ = max_fraction_under_var_jump(L, MU, SIGMA, R, T, LAM, J, ALPHA)
    assert abs(value_at_risk_jump(f_crash, MU, SIGMA, R, T, LAM, J, ALPHA) - L) < 1e-9
    assert f_crash < f_var
    assert value_at_risk_jump(f_var, MU, SIGMA, R, T, LAM, J, ALPHA) > L


def test_monte_carlo_var_with_crashes_matches_formula():
    f = 0.7
    rng = np.random.default_rng(0)
    W = simulate_wealth_jump(f, MU, SIGMA, R, T, LAM, J, 252, 50_000, rng)
    assert empirical_var(W[:, -1], 1.0, ALPHA) == pytest.approx(
        value_at_risk_jump(f, MU, SIGMA, R, T, LAM, J, ALPHA), abs=0.01)
