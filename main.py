"""Limits in Finance - runs all four modules, prints the results and saves the figures.

    python main.py
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.compounding import compound, compounding_error, continuous, error_estimate
from src.gbm import (gbm_log_return_cdf, ks_distance_to_gbm, simulate_walk_paths,
                     walk_log_return_distribution)
from src.option_pricing import binomial_price, black_scholes
from src.risk import (difference_quotient, empirical_var, growth_rate, kelly_fraction,
                      max_fraction_under_var, numerical_kelly_fraction, optimal_fraction,
                      simulate_wealth, value_at_risk)

FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
SEED = 42

# Module 1: savings account
P, R_SAVE, T_SAVE = 10_000.0, 0.10, 10.0
# Module 2: European call option
S0, K, R_OPT, SIGMA_OPT, T_OPT = 100.0, 100.0, 0.05, 0.20, 1.0
# Modules 3-4: stock market and risk-free asset (annual numbers)
MU, SIGMA, R, T = 0.10, 0.20, 0.03, 1.0
ALPHA, LOSS_LIMIT, F_MAX = 0.95, 0.15, 2.0


def header(title):
    print("\n" + "=" * 70 + "\n" + title + "\n" + "=" * 70)


def save(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    path = os.path.join(FIG_DIR, name)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  saved {os.path.relpath(path)}")


def module1_compounding():
    header("Module 1 - Continuous compounding:  lim (1 + r/n)^(nt) = e^(rt)")
    A_inf = continuous(P, R_SAVE, T_SAVE)
    print(f"  P = {P:,.0f}, r = {R_SAVE:.0%}, t = {T_SAVE:.0f} years  ->  P e^(rt) = {A_inf:,.4f}")
    print(f"  {'periods/year n':>16} {'A_n':>14} {'A - A_n':>12} {'estimate':>12}")
    for label, n in [("yearly", 1), ("quarterly", 4), ("monthly", 12), ("weekly", 52),
                     ("daily", 365), ("hourly", 8760), ("every minute", 525600)]:
        print(f"  {label:>16} {compound(P, R_SAVE, T_SAVE, n):14,.4f} "
              f"{compounding_error(P, R_SAVE, T_SAVE, n):12.4f} "
              f"{error_estimate(P, R_SAVE, T_SAVE, n):12.4f}")

    ns = np.unique(np.logspace(0, 7, 60).astype(int))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.semilogx(ns, compound(P, R_SAVE, T_SAVE, ns), "o-", ms=3, label=r"$A_n = P(1+r/n)^{nt}$")
    ax1.axhline(A_inf, color="C3", ls="--", label=r"limit $Pe^{rt}$")
    ax1.set(xlabel="compounding periods per year $n$", ylabel="final value",
            title="Discrete compounding converges to $Pe^{rt}$")
    ax1.legend()
    ax2.loglog(ns, compounding_error(P, R_SAVE, T_SAVE, ns), "o", ms=3, label="actual gap $A - A_n$")
    ax2.loglog(ns, error_estimate(P, R_SAVE, T_SAVE, ns), "-", label=r"estimate $Pe^{rt}\,r^2t/(2n)$")
    ax2.set(xlabel="$n$", ylabel="gap", title="The gap shrinks like $1/n$")
    ax2.legend()
    save(fig, "fig1_compounding.png")


def module2_options():
    header("Module 2 - Binomial tree -> Black-Scholes:  lim_{N->inf} C_N = C_BS")
    bs = black_scholes(S0, K, R_OPT, SIGMA_OPT, T_OPT)
    print(f"  S = K = {S0:.0f}, r = {R_OPT:.0%}, sigma = {SIGMA_OPT:.0%}, T = {T_OPT:.0f}  ->  C_BS = {bs:.6f}")
    print(f"  {'N':>6} {'C_N':>12} {'|C_N - C_BS|':>14}")
    for N in [1, 2, 5, 10, 50, 100, 500, 1000]:
        c = binomial_price(S0, K, R_OPT, SIGMA_OPT, T_OPT, N)
        print(f"  {N:>6} {c:12.6f} {abs(c - bs):14.6f}")

    Ns = np.arange(1, 201)
    prices = np.array([binomial_price(S0, K, R_OPT, SIGMA_OPT, T_OPT, N) for N in Ns])
    Ns_log = np.unique(np.logspace(0, 3.3, 40).astype(int))
    errors = np.array([abs(binomial_price(S0, K, R_OPT, SIGMA_OPT, T_OPT, N) - bs) for N in Ns_log])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.plot(Ns, prices, "-", lw=1, label="binomial price $C_N$")
    ax1.axhline(bs, color="C3", ls="--", label=f"Black-Scholes {bs:.4f}")
    ax1.set(xlabel="tree steps $N$", ylabel="call price", title="Binomial price converges to Black-Scholes")
    ax1.legend()
    ax2.loglog(Ns_log, errors, "o", ms=4, label="$|C_N - C_{BS}|$")
    ax2.loglog(Ns_log, errors[0] / Ns_log, "--", label="reference $\\propto 1/N$")
    ax2.set(xlabel="$N$", ylabel="error", title="Error goes to 0 like $1/N$")
    ax2.legend()
    save(fig, "fig2_binomial_to_black_scholes.png")


def module3_random_walk():
    header("Module 3 - Random walk -> Geometric Brownian Motion as dt -> 0")
    print(f"  mu = {MU:.0%}, sigma = {SIGMA:.0%}, T = {T:.0f} year")
    print(f"  {'steps n':>8} {'dt (years)':>12} {'KS distance':>12}")
    for n in [1, 4, 12, 52, 252, 2520, 25200]:
        print(f"  {n:>8} {T / n:12.6f} {ks_distance_to_gbm(MU, SIGMA, T, n):12.6f}")

    rng = np.random.default_rng(SEED)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, n in zip(axes[0], [12, 252]):
        paths = simulate_walk_paths(S0, MU, SIGMA, T, n, 30, rng)
        ax.plot(np.linspace(0, T, n + 1), paths.T, lw=0.8, alpha=0.8)
        ax.set(xlabel="time (years)", ylabel="price", title=f"Random-walk paths, n = {n} steps")

    ax = axes[1, 0]
    x_grid = np.linspace(-0.7, 0.8, 400)
    for n in [2, 12, 252]:
        x, pmf = walk_log_return_distribution(MU, SIGMA, T, n)
        # pad with F = 0 on the left and F = 1 on the right so the step CDF spans the axis
        ax.step(np.concatenate(([-1.0], x, [1.0])), np.concatenate(([0.0], np.cumsum(pmf), [1.0])),
                where="post", lw=1.2, label=f"walk, n = {n}")
    ax.plot(x_grid, gbm_log_return_cdf(x_grid, MU, SIGMA, T), "k--", lw=1.5, label="GBM limit")
    ax.set(xlim=(-0.7, 0.8), xlabel=r"log-return $\log(S_T/S_0)$", ylabel="CDF",
           title="Distribution converges to the GBM normal law")
    ax.legend()

    ax = axes[1, 1]
    ns = np.unique(np.logspace(0, 4.5, 30).astype(int))
    ks = np.array([ks_distance_to_gbm(MU, SIGMA, T, n) for n in ns])
    ax.loglog(ns, ks, "o", ms=4, label="KS distance")
    ax.loglog(ns, ks[-1] * np.sqrt(ns[-1] / ns), "--", label=r"reference $\propto 1/\sqrt{n}$")
    ax.set(xlabel="steps $n = T/dt$", ylabel=r"$\sup_x |F_n(x) - F(x)|$", title="KS distance goes to 0")
    ax.legend()
    save(fig, "fig3_random_walk_to_gbm.png")


def module4_risk():
    header("Module 4 - Risk optimization: best growth under a Value-at-Risk limit")
    print(f"  stock: mu = {MU:.0%}, sigma = {SIGMA:.0%};  risk-free r = {R:.0%};  horizon T = {T:.0f} year")
    print(f"  rule: with {ALPHA:.0%} confidence, lose at most {LOSS_LIMIT:.0%} in one year;  leverage <= {F_MAX:.0f}x")

    g = lambda f: growth_rate(f, MU, SIGMA, R)
    f_test = 1.0
    print(f"\n  Derivative as a limit: g'({f_test}) = (mu - r) - f sigma^2 = {(MU - R) - f_test * SIGMA**2:.6f}")
    print(f"  {'h':>10} {'(g(f+h) - g(f))/h':>20}")
    for h in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
        print(f"  {h:>10.0e} {difference_quotient(g, f_test, h):20.10f}")

    f_star = kelly_fraction(MU, SIGMA, R)
    f_star_num, _ = numerical_kelly_fraction(MU, SIGMA, R)
    print(f"\n  Max growth (Kelly):  f* = (mu - r)/sigma^2 = {f_star:.6f}   numerical (bisection on g') = {f_star_num:.6f}")
    print(f"  growth at f*: g(f*) = {g(f_star):.4%} per year,  but VaR(f*) = {float(value_at_risk(f_star, MU, SIGMA, R, T, ALPHA)):.2%}")

    f_var, hist = max_fraction_under_var(LOSS_LIMIT, MU, SIGMA, R, T, ALPHA, F_MAX)
    print(f"\n  IVT + bisection on h(f) = VaR(f) - L over [0, {F_MAX}]:")
    print(f"  {'k':>4} {'a':>10} {'b':>10} {'mid':>10} {'h(mid)':>12}")
    for k, (a, b, mid, hm) in enumerate(hist):
        if k < 8 or k == len(hist) - 1:
            print(f"  {k + 1:>4} {a:10.6f} {b:10.6f} {mid:10.6f} {hm:12.2e}")
        elif k == 8:
            print(f"  {'...':>4}")
    f_opt = optimal_fraction(MU, SIGMA, R, T, LOSS_LIMIT, ALPHA, F_MAX)
    print(f"\n  f_VaR = {f_var:.6f}  ->  f_opt = min(f*, f_VaR, f_max) = {f_opt:.6f}")
    print(f"  => invest {f_opt:.1%} in the stock and {1 - f_opt:.1%} risk-free;"
          f" growth {g(f_opt):.2%}/yr, VaR {float(value_at_risk(f_opt, MU, SIGMA, R, T, ALPHA)):.2%}")

    rng = np.random.default_rng(SEED)
    print(f"\n  Monte Carlo check (20,000 paths): VaR with rebalancing every dt -> analytic {LOSS_LIMIT:.2%}")
    for label, n in [("yearly", 1), ("quarterly", 4), ("monthly", 12), ("weekly", 52), ("daily", 252)]:
        W = simulate_wealth(f_opt, MU, SIGMA, R, T, n, 20_000, rng)
        print(f"  {label:>10}: VaR = {empirical_var(W[:, -1], 1.0, ALPHA):.2%}")

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    f_grid = np.linspace(0, 2 * f_star, 400)

    ax = axes[0, 0]
    ax.plot(f_grid, g(f_grid), lw=2)
    ax.axvline(f_star, color="C3", ls="--", label=f"Kelly $f^*$ = {f_star:.2f}")
    ax.axvline(f_opt, color="C2", ls="-", label=f"chosen $f_{{opt}}$ = {f_opt:.2f}")
    ax.set(xlabel="fraction in stock $f$", ylabel="growth rate $g(f)$", title="Growth: continuous, concave, max at $f^*$")
    ax.legend()

    ax = axes[0, 1]
    var_grid = value_at_risk(f_grid, MU, SIGMA, R, T, ALPHA)
    ax.plot(f_grid, var_grid, lw=2, label=f"VaR$_{{{ALPHA:.0%}}}(f)$".replace("%", r"\%"))
    ax.axhline(LOSS_LIMIT, color="k", ls=":", label=f"limit L = {LOSS_LIMIT:.0%}")
    ax.axvline(f_var, color="C2", label=f"IVT root $f_{{VaR}}$ = {f_var:.3f}")
    ax.fill_between(f_grid, var_grid.min(), var_grid.max(), where=f_grid <= f_var, color="C2", alpha=0.1)
    ax.set(xlabel="fraction in stock $f$", ylabel="1-year VaR", title="Risk: VaR(f) - L changes sign -> root exists")
    ax.legend()

    ax = axes[1, 0]
    t_grid = np.linspace(0, T, 253)
    for f, color, name in [(f_star, "C3", "Kelly $f^*$"), (f_opt, "C2", "$f_{opt}$")]:
        W = simulate_wealth(f, MU, SIGMA, R, T, 252, 5_000, rng)
        lo, mid, hi = np.quantile(W, [1 - ALPHA, 0.5, ALPHA], axis=0)
        ax.fill_between(t_grid, lo, hi, color=color, alpha=0.2)
        ax.plot(t_grid, mid, color=color, label=f"{name}: median and 5-95% band")
    ax.set(xlabel="time (years)", ylabel="wealth (start = 1)", title="Simulated wealth, daily rebalancing")
    ax.legend()

    ax = axes[1, 1]
    W = simulate_wealth(f_opt, MU, SIGMA, R, T, 252, 20_000, rng)
    loss = 1.0 - W[:, -1]
    ax.hist(loss, bins=80, density=True, color="C0", alpha=0.7)
    ax.axvline(LOSS_LIMIT, color="k", ls=":", label=f"analytic VaR = {LOSS_LIMIT:.1%}")
    ax.axvline(empirical_var(W[:, -1], 1.0, ALPHA), color="C3", ls="--",
               label=f"simulated VaR = {empirical_var(W[:, -1], 1.0, ALPHA):.1%}")
    ax.set(xlabel="1-year loss (negative = gain)", ylabel="density", title=f"Loss distribution at $f_{{opt}}$ = {f_opt:.2f}")
    ax.legend()
    save(fig, "fig4_risk_optimization.png")


if __name__ == "__main__":
    module1_compounding()
    module2_options()
    module3_random_walk()
    module4_risk()
