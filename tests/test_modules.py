"""Checks every module against known values and the limits it is supposed to show."""
import math

import numpy as np
import pytest

from src.bisection import bisection
from src.compounding import compound, compounding_error, continuous, error_estimate
from src.gbm import ks_distance_to_gbm, up_probability, walk_log_return_distribution
from src.normal import norm_cdf, norm_ppf
from src.option_pricing import binomial_price, black_scholes
from src.portfolio import (best_mix_on_grid, best_mix_under_var, fund_fraction, fund_fractions, mix_grid,
                           near_optimal_ranges, portfolio_es, portfolio_growth, portfolio_var, refine_best_mix,
                           simulate_wealth_mix)
from src.risk import (binomial_pmf, empirical_es, empirical_var, expected_shortfall, expected_shortfall_jump,
                      growth_rate, growth_rate_jump, jump_log_wealth_cdf, kelly_fraction, kelly_fraction_jump,
                      max_fraction_under_es, max_fraction_under_var, max_fraction_under_var_jump,
                      numerical_kelly_fraction, optimal_fraction, optimal_fraction_jump, poisson_pmf,
                      simulate_wealth, simulate_wealth_jump, value_at_risk, value_at_risk_jump, var_root_formula)

MU, SIGMA, R, T, ALPHA, L = 0.10, 0.20, 0.03, 1.0, 0.95, 0.15
LAM, J = 0.10, -0.25
MU3 = np.array([MU, math.log(1.05), math.log(1.045)])
SIG3 = np.array([SIGMA, 0.15, 0.06])
RHO3 = np.array([[1.0, 0.0, -0.2], [0.0, 1.0, 0.3], [-0.2, 0.3, 1.0]])


# --- bisection / normal -------------------------------------------------------

def test_bisection_finds_sqrt2():
    root, history = bisection(lambda x: x**2 - 2, 0.0, 2.0)
    assert root == pytest.approx(math.sqrt(2), abs=1e-11)
    widths = history[:, 1] - history[:, 0]
    assert np.all(np.diff(widths) < 0)          # interval shrinks every step


def test_bisection_needs_sign_change():
    with pytest.raises(ValueError):
        bisection(lambda x: x**2 + 1, -1.0, 1.0)


def test_bisection_sees_the_signs_of_tiny_values():
    # fun(a) * fun(mid) would round to 0 here and hide the sign change
    root, _ = bisection(lambda x: 1e-200 * (x - 0.3), 0.0, 1.0)
    assert root == pytest.approx(0.3, abs=1e-11)


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


def test_walk_step_must_give_a_probability():
    with pytest.raises(ValueError):
        up_probability(0.30, 0.10, 1.0)         # one step a year is too long here: P(up) would be above 1


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


def test_var_root_is_unique_for_any_horizon():
    for horizon in [1.0, 10.0, 30.0]:           # past 22 years VaR(f) first falls, but it still crosses L once
        f_var, _ = max_fraction_under_var(L, MU, SIGMA, R, horizon, ALPHA, f_max=5.0)
        assert f_var == pytest.approx(var_root_formula(L, MU, SIGMA, R, horizon, ALPHA), abs=1e-9)
        f = np.linspace(0, 5, 2001)
        ok = value_at_risk(f, MU, SIGMA, R, horizon, ALPHA) <= L
        assert np.all(ok[f < f_var - 1e-6]) and not np.any(ok[f > f_var + 1e-6])


def test_expected_shortfall():
    f = np.array([0.0, 0.3, 0.705, 1.5])
    es, var = expected_shortfall(f, MU, SIGMA, R, T, ALPHA), value_at_risk(f, MU, SIGMA, R, T, ALPHA)
    assert es[0] == pytest.approx(var[0]) and np.all(es[1:] > var[1:])    # bad years lose more than the VaR on average
    for x in [0.3, 0.705, 1.5]:
        assert expected_shortfall_jump(x, MU, SIGMA, R, T, 0.0, J, ALPHA) == pytest.approx(
            float(expected_shortfall(x, MU, SIGMA, R, T, ALPHA)), abs=1e-9)
    assert expected_shortfall_jump(0.705, MU, SIGMA, R, T, LAM, J, ALPHA) > float(expected_shortfall(0.705, MU, SIGMA, R, T, ALPHA))
    f_es, _ = max_fraction_under_es(L, MU, SIGMA, R, T, ALPHA)
    assert abs(float(expected_shortfall(f_es, MU, SIGMA, R, T, ALPHA)) - L) < 1e-9
    assert f_es < max_fraction_under_var(L, MU, SIGMA, R, T, ALPHA)[0]


def test_monte_carlo_es_matches_formula():
    rng = np.random.default_rng(0)
    W = simulate_wealth(0.7, MU, SIGMA, R, T, 52, 50_000, rng)
    assert empirical_es(W[:, -1], 1.0, ALPHA) == pytest.approx(float(expected_shortfall(0.7, MU, SIGMA, R, T, ALPHA)), abs=0.01)
    W = simulate_wealth_jump(0.7, MU, SIGMA, R, T, LAM, J, 52, 50_000, rng)
    assert empirical_es(W[:, -1], 1.0, ALPHA) == pytest.approx(
        expected_shortfall_jump(0.7, MU, SIGMA, R, T, LAM, J, ALPHA), abs=0.01)


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


def test_jump_cdf_without_noise():
    # f = 0: everything in the bank, so log-wealth is exactly r T (no division by zero)
    assert jump_log_wealth_cdf(R * T - 1e-9, 0.0, MU, SIGMA, R, T, LAM, J) == 0.0
    assert jump_log_wealth_cdf(R * T + 1e-9, 0.0, MU, SIGMA, R, T, LAM, J) == pytest.approx(1.0)


def test_crashes_lower_growth_and_the_kelly_fraction():
    f = np.linspace(0.1, 2.0, 20)
    assert np.allclose(growth_rate_jump(f, MU, SIGMA, R, 0.0, J), growth_rate(f, MU, SIGMA, R))
    assert np.all(growth_rate_jump(f, MU, SIGMA, R, LAM, J) < growth_rate(f, MU, SIGMA, R))   # same average, lower growth
    assert kelly_fraction_jump(MU, SIGMA, R, 0.0, J)[0] == pytest.approx(1.75, abs=1e-8)
    f_star, _ = kelly_fraction_jump(MU, SIGMA, R, LAM, J)
    slope = (growth_rate_jump(f_star + 1e-6, MU, SIGMA, R, LAM, J) - growth_rate_jump(f_star - 1e-6, MU, SIGMA, R, LAM, J)) / 2e-6
    assert abs(slope) < 1e-8 and f_star == pytest.approx(1.4098, abs=1e-4)
    f_crash, _ = max_fraction_under_var_jump(L, MU, SIGMA, R, T, LAM, J, ALPHA)
    assert optimal_fraction_jump(MU, SIGMA, R, T, L, LAM, J, ALPHA) == pytest.approx(f_crash)


# --- module 5: more assets ------------------------------------------------------

def test_one_asset_mix_is_module4():
    for f in [0.0, 0.3, 0.705, 1.0]:
        assert portfolio_var([f, 0, 0], MU3, SIG3, RHO3, R, T, ALPHA) == pytest.approx(
            float(value_at_risk(f, MU, SIGMA, R, T, ALPHA)), abs=1e-12)
    f, _, _ = fund_fraction([1, 0, 0], L, MU3, SIG3, RHO3, R, T, ALPHA)
    assert f == pytest.approx(max_fraction_under_var(L, MU, SIGMA, R, T, ALPHA, 1.0)[0], abs=1e-12)


def test_mix_grid_covers_the_triangle():
    grid = mix_grid(3, 100)
    assert len(grid) == 5151 and np.allclose(grid.sum(axis=1), 1) and grid.min() >= 0
    # first entry slowest: the same order as matlab/mix_grid.m and MIX_GRID in the web pages
    assert np.array_equal(mix_grid(3, 2), [[0, 0, 1], [0, 0.5, 0.5], [0, 1, 0], [0.5, 0, 0.5], [0.5, 0.5, 0], [1, 0, 0]])


def test_formula_grid_matches_the_bisection_grid():
    w, f, d = best_mix_under_var(L, MU3, SIG3, RHO3, R, T, ALPHA)
    w2, f2, d2, _ = best_mix_on_grid(L, MU3, SIG3, RHO3, R, T, ALPHA, 100)
    assert np.allclose(w, w2, atol=1e-10) and np.array_equal(d, d2)
    D = np.random.default_rng(2).dirichlet(np.ones(3), 200)
    fs, gs = fund_fractions(D, L, MU3, SIG3, RHO3, R, T, ALPHA)
    for d_, f_, g_ in zip(D, fs, gs):
        f1, g1, _ = fund_fraction(d_, L, MU3, SIG3, RHO3, R, T, ALPHA)
        assert f_ == pytest.approx(f1, abs=1e-10) and g_ == pytest.approx(g1, abs=1e-12)


def test_finer_grids_converge_to_the_best_mix():
    gs = [best_mix_on_grid(L, MU3, SIG3, RHO3, R, T, ALPHA, s)[3] for s in [25, 50, 100, 200, 400]]
    assert np.all(np.diff(gs) >= 0)                                 # each grid contains the last: never worse
    w, f, d, g = refine_best_mix(L, MU3, SIG3, RHO3, R, T, ALPHA)
    assert g >= gs[-1] and g == pytest.approx(0.0746887, abs=2e-7)
    assert portfolio_var(w, MU3, SIG3, RHO3, R, T, ALPHA) <= L + 1e-9 and w.min() >= 0 and w.sum() <= 1 + 1e-9
    for x in np.random.default_rng(3).dirichlet(np.ones(4), 3000)[:, :3]:      # no random safe mix beats it
        if portfolio_var(x, MU3, SIG3, RHO3, R, T, ALPHA) <= L:
            assert portfolio_growth(x, MU3, SIG3, RHO3, R)[0] <= g + 1e-12


def test_near_optimal_mixes_form_a_flat_band():
    lo, hi, n = near_optimal_ranges(L, MU3, SIG3, RHO3, R, T, ALPHA, 100, 1e-4)
    w, _, _ = best_mix_under_var(L, MU3, SIG3, RHO3, R, T, ALPHA)
    assert n > 1 and np.all(lo <= w + 1e-12) and np.all(w <= hi + 1e-12)
    assert hi[1] - lo[1] >= 0.05 and hi[0] - lo[0] <= 0.01        # stocks pinned near 73%, gold anywhere in a band


def test_one_asset_es_is_module4():
    for f in [0.0, 0.3, 0.705, 1.0]:
        assert portfolio_es([f, 0, 0], MU3, SIG3, RHO3, R, T, ALPHA) == pytest.approx(
            float(expected_shortfall(f, MU, SIGMA, R, T, ALPHA)), abs=1e-12)


def test_best_mix_respects_the_limit_and_beats_stocks_only():
    w, f, d = best_mix_under_var(L, MU3, SIG3, RHO3, R, T, ALPHA)
    assert np.allclose(w, f * d) and w.min() >= 0 and w.sum() <= 1 + 1e-12
    assert portfolio_var(w, MU3, SIG3, RHO3, R, T, ALPHA) <= L + 1e-9
    g_best = portfolio_growth(w, MU3, SIG3, RHO3, R)[0]
    for i in range(3):                                  # better than any single asset + bank
        assert g_best >= fund_fraction(np.eye(3)[i], L, MU3, SIG3, RHO3, R, T, ALPHA)[1]
    rng = np.random.default_rng(1)                      # and no random mix within the limit does better
    for x in rng.dirichlet(np.ones(4), 2000)[:, :3]:
        if portfolio_var(x, MU3, SIG3, RHO3, R, T, ALPHA) <= L:
            assert portfolio_growth(x, MU3, SIG3, RHO3, R)[0] <= g_best + 2e-4


def test_best_mix_is_stocks_only_when_others_do_not_pay():
    mu = np.array([MU, R, R])                           # gold and bonds earn only the bank rate
    w, f, _ = best_mix_under_var(L, mu, SIG3, np.eye(3), R, T, ALPHA)
    assert np.allclose(w, [f, 0, 0]) and f == pytest.approx(0.705014, abs=1e-6)


def test_monte_carlo_mix_var():
    w = np.array([0.5, 0.1, 0.3])
    W = simulate_wealth_mix(w, MU3, SIG3, RHO3, R, T, 52, 50_000, np.random.default_rng(0))
    assert empirical_var(W[:, -1], 1.0, ALPHA) == pytest.approx(portfolio_var(w, MU3, SIG3, RHO3, R, T, ALPHA), abs=0.01)
