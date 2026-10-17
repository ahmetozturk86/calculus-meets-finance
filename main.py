"""Limits in Finance - runs all modules, prints the results and saves the figures.

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
from src.portfolio import (best_mix_under_var, fund_fraction, mix_grid, portfolio_growth, portfolio_var,
                           simulate_wealth_mix)
from src.risk import (binomial_pmf, difference_quotient, empirical_var, growth_rate, kelly_fraction,
                      jump_log_wealth_cdf, max_fraction_under_var, max_fraction_under_var_jump, numerical_kelly_fraction,
                      optimal_fraction, poisson_pmf, simulate_wealth, simulate_wealth_jump, value_at_risk,
                      value_at_risk_jump)

FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
SEED = 42

# Module 1: savings account
P, R_SAVE, T_SAVE = 10_000.0, 0.10, 10.0
# Module 2: European call option
S0, K, R_OPT, SIGMA_OPT, T_OPT = 100.0, 100.0, 0.05, 0.20, 1.0
# Modules 3-4: stock market and risk-free asset (annual numbers)
MU, SIGMA, R, T = 0.10, 0.20, 0.03, 1.0
ALPHA, LOSS_LIMIT, F_MAX = 0.95, 0.15, 2.0
# Module 4b: rare crashes (fat tails): the stock drops 25%, on average once every 10 years
CRASH_RATE, CRASH_SIZE = 0.10, -0.25
# Module 5: more assets. Stocks as above; gold and bonds are illustrative long-run numbers (5% and 4.5% a year as
# simple rates, so ln(1.05) and ln(1.045) as continuous rates), not data. Same numbers as the "Normal" market on the website.
ASSET_NAMES = ["stocks", "gold", "bonds"]
MU_ASSETS = np.array([MU, np.log(1.05), np.log(1.045)])
SIGMA_ASSETS = np.array([SIGMA, 0.15, 0.06])
RHO_ASSETS = np.array([[1.0, 0.0, -0.2], [0.0, 1.0, 0.3], [-0.2, 0.3, 1.0]])


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


def module4b_crashes():
    header("Module 4b - Fat tails: the bell curve plus rare crashes")
    lam, J = CRASH_RATE, CRASH_SIZE
    print(f"  a crash makes the stock fall {-J:.0%}; it happens on average {lam:g} times a year"
          f" (once every {1 / lam:.0f} years); the stock's average return is still {MU:.0%}")

    print(f"\n  Limit: crash coin-flips every dt = 1/n year -> Poisson law of the number of crashes")
    print(f"  {'flips/year n':>14} {'max |Binomial(n, lam/n) - Poisson(lam)|':>40}")
    for n in [1, 12, 52, 252, 10_000]:
        err = max(abs(binomial_pmf(k, n, lam / n) - poisson_pmf(k, lam)) for k in range(31))
        print(f"  {n:>14,} {err:40.2e}")

    f_var, _ = max_fraction_under_var(LOSS_LIMIT, MU, SIGMA, R, T, ALPHA, F_MAX)
    print(f"\n  Value-at-Risk at the normal answer f_VaR = {f_var:.6f}:")
    print(f"  {'confidence':>12} {'bell curve only':>16} {'with crashes':>14}")
    for a in [0.95, 0.99]:
        print(f"  {a:>12.0%} {float(value_at_risk(f_var, MU, SIGMA, R, T, a)):16.2%}"
              f" {value_at_risk_jump(f_var, MU, SIGMA, R, T, lam, J, a):14.2%}")

    f_crash, hist = max_fraction_under_var_jump(LOSS_LIMIT, MU, SIGMA, R, T, lam, J, ALPHA, F_MAX)
    print(f"\n  IVT + bisection on h(f) = VaR_crash(f) - L (each VaR_crash is itself a bisection on the CDF):")
    print(f"  {'k':>4} {'a':>10} {'b':>10} {'mid':>10} {'h(mid)':>12}")
    for k, (a, b, mid, hm) in enumerate(hist):
        if k < 5 or k == len(hist) - 1:
            print(f"  {k + 1:>4} {a:10.6f} {b:10.6f} {mid:10.6f} {hm:12.2e}")
        elif k == 5:
            print(f"  {'...':>4}")
    print(f"\n  f_crash = {f_crash:.6f}  ->  with crashes, invest {f_crash:.1%} in the stock instead of {f_var:.1%}")

    rng = np.random.default_rng(SEED)
    print(f"\n  Monte Carlo check (40,000 paths, crash coin-flips each step): VaR at f_crash -> {LOSS_LIMIT:.2%}")
    for label, n in [("yearly", 1), ("monthly", 12), ("weekly", 52), ("daily", 252)]:
        W = simulate_wealth_jump(f_crash, MU, SIGMA, R, T, lam, J, n, 40_000, rng)
        print(f"  {label:>10}: VaR = {empirical_var(W[:, -1], 1.0, ALPHA):.2%}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    losses = np.array([0.10, 0.20, 0.30])           # chance of losing more than this in one year, at f_VaR
    p_bell = [jump_log_wealth_cdf(np.log(1 - l), f_var, MU, SIGMA, R, T, 0.0, J) for l in losses]
    p_crash = [jump_log_wealth_cdf(np.log(1 - l), f_var, MU, SIGMA, R, T, lam, J) for l in losses]
    pos = np.arange(len(losses))
    for shift, p, color, name in [(-0.2, p_bell, "C0", "bell curve only"), (0.2, p_crash, "C3", "with rare crashes")]:
        bars = ax.bar(pos + shift, np.array(p) * 100, width=0.4, color=color, label=name)
        ax.bar_label(bars, labels=[f"{v:.1%}" for v in p], fontsize=9)
    ax.set_xticks(pos, [f"lose more than {l:.0%}" for l in losses])
    ax.set(ylabel="chance in one year (%)", title=f"Same average return, more big losses (f = {f_var:.2f})")
    ax.legend()

    ax = axes[1]
    f_grid = np.linspace(0, 1.5, 61)
    ax.plot(f_grid, value_at_risk(f_grid, MU, SIGMA, R, T, ALPHA), color="C0", lw=2, label="bell curve only")
    ax.plot(f_grid, [value_at_risk_jump(v, MU, SIGMA, R, T, lam, J, ALPHA) for v in f_grid], color="C3", lw=2,
            label="with rare crashes")
    ax.axhline(LOSS_LIMIT, color="k", ls=":", label=f"limit L = {LOSS_LIMIT:.0%}")
    ax.axvline(f_var, color="C0", ls="--", label=f"$f_{{VaR}}$ = {f_var:.3f}")
    ax.axvline(f_crash, color="C3", ls="--", label=f"$f_{{crash}}$ = {f_crash:.3f}")
    ax.set(xlabel="fraction in stock $f$", ylabel="1-year VaR", title="Crashes move the IVT root to the left")
    ax.legend()
    save(fig, "fig5_fat_tails.png")


def module5_portfolio():
    header("Module 5 - More assets: the best mix of stocks, gold, bonds and the bank")
    for name, m, sd in zip(ASSET_NAMES, MU_ASSETS, SIGMA_ASSETS):
        print(f"  {name:>7}: mu = {m:.2%}, sigma = {sd:.0%}")
    print(f"  bank r = {R:.0%};  rule: with {ALPHA:.0%} confidence lose at most {LOSS_LIMIT:.0%} in one year; no borrowing")
    print("  correlations (stocks, gold, bonds):", RHO_ASSETS.tolist())

    print("\n  One risky asset + bank (Module 4 for each asset, f <= 1):")
    print(f"  {'asset':>7} {'f':>10} {'growth':>10}")
    for i, name in enumerate(ASSET_NAMES):
        f, g, _ = fund_fraction(np.eye(3)[i], LOSS_LIMIT, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA)
        print(f"  {name:>7} {f:10.6f} {g:10.4%}")

    grid = mix_grid(3)
    w, f, d = best_mix_under_var(LOSS_LIMIT, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA)
    g, sd = portfolio_growth(w, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R)
    print(f"\n  Every mix d of the risky assets is one fund; Module 4 (IVT + bisection) gives each fund's f.")
    print(f"  Tried {len(grid):,} funds (1% steps). Fastest-growing one under the limit:")
    print("  " + ", ".join(f"{name} {x:.0%}" for name, x in zip(ASSET_NAMES, w)) + f", bank {1 - w.sum():.0%}")
    print(f"  growth {g:.4%}/yr, VaR {portfolio_var(w, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA):.4%}")
    f4, g4, _ = fund_fraction(np.eye(3)[0], LOSS_LIMIT, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA)
    print(f"  vs stocks + bank (Module 4): growth {g4:.4%}/yr  ->  +{g - g4:.2%} a year for the same risk")

    rng = np.random.default_rng(SEED)
    print("\n  Monte Carlo check (40,000 one-year paths, rebalanced daily): VaR of the best mix")
    W = simulate_wealth_mix(w, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, 252, 40_000, rng)
    print(f"       daily: VaR = {empirical_var(W[:, -1], 1.0, ALPHA):.2%}  (analytic {portfolio_var(w, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA):.2%})")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    colors = ["C0", "#c08a00", "C2", "0.6"]
    plans = [("stocks + bank\n(Module 4)", np.array([f4, 0, 0])), ("best mix\n(Module 5)", w)]
    for k, (label, wk) in enumerate(plans):
        left = 0.0
        for i, part in enumerate(list(wk) + [1 - wk.sum()]):
            ax.barh(k, part * 100, left=left, color=colors[i], label=(ASSET_NAMES + ["bank"])[i] if k == 1 else None)
            if part >= 0.06:
                ax.text(left + part * 50, k, f"{part:.0%}", ha="center", va="center", color="white", fontsize=10)
            left += part * 100
        gk = portfolio_growth(wk, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R)[0]
        ax.text(101, k, f"growth {gk:.2%}", va="center", fontsize=10)
    ax.set_yticks([0, 1], [p[0] for p in plans])
    ax.set(xlim=(0, 130), ylim=(-0.6, 1.9), xlabel="share of the money (%)", title=f"Same risk (VaR = {LOSS_LIMIT:.0%}), more growth")
    ax.legend(loc="upper center", ncol=4, fontsize=9)

    ax = axes[1]
    sample = rng.dirichlet(np.ones(4), 3000)[:, :3]                  # random mixes of the 3 assets + bank
    pts = np.array([[portfolio_var(x, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA),
                     portfolio_growth(x, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R)[0]] for x in sample])
    ok = pts[:, 0] <= LOSS_LIMIT
    ax.scatter(pts[~ok, 0], pts[~ok, 1], s=6, color="0.75", label="breaks the limit")
    ax.scatter(pts[ok, 0], pts[ok, 1], s=6, color="C2", alpha=0.6, label="within the limit")
    ax.axvline(LOSS_LIMIT, color="k", ls=":", label=f"limit L = {LOSS_LIMIT:.0%}")
    ax.scatter([LOSS_LIMIT], [g4], marker="^", s=90, color="C0", zorder=3, label="stocks + bank")
    ax.scatter([portfolio_var(w, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA)], [g], marker="*", s=220,
               color="C3", zorder=3, label="best mix")
    ax.set(xlabel="bad-year loss (1-year VaR)", ylabel="growth rate g", title="Random mixes: the best is the highest green one")
    ax.legend(fontsize=9)
    save(fig, "fig6_best_mix.png")


if __name__ == "__main__":
    module1_compounding()
    module2_options()
    module3_random_walk()
    module4_risk()
    module4b_crashes()
    module5_portfolio()
