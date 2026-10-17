# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

School project + GitHub portfolio piece: finance modules, each showing a discrete model converging to a continuous one (limits), ending in a VaR-constrained portfolio allocation found via the Intermediate Value Theorem + bisection (Module 4: stocks + bank; Module 5: the best mix of stocks, gold, bonds and bank). The school requires **MATLAB**; Python is the reference implementation.

## Commands

```bash
pip install -r requirements.txt           # numpy, matplotlib, pytest
python main.py                            # prints all result tables, writes figures/*.png
python -m pytest                          # all tests (pytest.ini sets pythonpath=. and testpaths=tests)
python -m pytest tests/test_modules.py::test_var_root_hits_the_limit   # single test
```

MATLAB (from the `matlab/` folder): `run_tests` then `main` (figures go to `matlab/figures/`). MATLAB is not installed on the dev machine, so MATLAB changes cannot be run locally — review them carefully.

## Architecture

- `src/bisection.py` is the shared core: `bisection(fun, a, b)` returns `(root, history)` and raises `ValueError` when the IVT sign condition fails. It is used by `normal.norm_ppf` (inverse normal CDF), `risk.numerical_kelly_fraction` (root of g'), `risk.max_fraction_under_var` (root of VaR(f) − L), `risk.value_at_risk_jump` (quantile of the crash mixture) `risk.max_fraction_under_var_jump` (root of VaR_jump(f) − L, so a bisection whose function is itself a bisection) and `portfolio.fund_fraction` (one bisection per fund in Module 5). The history rows `[a, b, mid, fun(mid)]` are printed in `main.py` as the IVT demonstration.
- `src/normal.py` implements the normal CDF/PDF/PPF with `math.erf` on purpose — no scipy — so code ports to base MATLAB without the Statistics Toolbox.
- Module chain: `compounding.py` (1) → `option_pricing.py` (2) → `gbm.py` (3) → `risk.py` (4). Module 4's closed-form wealth/VaR formulas assume the continuous limits established in modules 1 and 3. Module 4b (fat tails, `main.module4b_crashes`) adds Merton-style crashes: size J, rate λ, drift compensated by λJ so the average return stays μ. With λ = 0 every crash function must equal its Module 4 twin (tested). Module 5 (`portfolio.py`, `main.module5_portfolio`): a fixed mix d of the risky assets is one fund (μ_d = r + d·(μ − r), σ_d = √(dᵀCd)), so `fund_fraction` is Module 4 for that fund (f = min(f*, f_VaR, 1)); `best_mix_under_var` tries every mix on the 1% grid `mix_grid(3)` (5,151 funds, fixed order) and keeps the fastest-growing one. No borrowing, no short selling. A stocks-only mix must equal Module 4 (tested). Gold/bond numbers in `main.py` match the website's "Normal" market.
- `gbm.ks_distance_to_gbm` computes the KS distance **exactly** from binomial probabilities (no Monte Carlo), so its convergence plot is noise-free.
- `main.py` holds all model parameters as constants at the top (savings, option, market, risk-limit) and one function per module that prints tables and saves one figure.

## Python ↔ MATLAB parity

`matlab/` mirrors `src/` one function per `.m` file with the **same names and argument order**; `matlab/main.m` mirrors `main.py` and `matlab/run_tests.m` mirrors `tests/test_modules.py`. When changing a Python function, make the same change in its `.m` file (and vice versa). Differences to keep in mind:
- Python random functions take an `rng` argument; MATLAB versions drop it and use the global stream (`rng(seed)`), so Monte Carlo numbers differ between languages while deterministic numbers must match.
- `quantile_linear.m` reimplements `numpy.quantile` (linear interpolation) to avoid the Statistics Toolbox.
- Keep Python code MATLAB-portable: plain functions, numpy arrays only, no classes/pandas/scipy.

## Web pages

`web/index.html` (explainer, uses Chart.js from a CDN) and `web/investor.html` (real-life tool, no external libraries, SVG charts) re-implement the deterministic formulas in JavaScript: `bisection`, `normPpf`, `growthRate`, `kellyFraction`, `valueAtRisk`, `poissonPmf`, `jumpLogWealthCdf`, `valueAtRiskJump`, `maxFractionUnderVar`, `quantileLinear`, `fundOf`, `fundFraction`, `bestMix` and `portfolioVar` mirror `src/bisection.py`, `src/normal.py`, `src/risk.py` and `src/portfolio.py`. Keep them in sync when the Python changes. The JS numbers must match Python, for example f_VaR = 0.705014 with the default parameters and f_max = 2, f_crash = 0.617447 with λ = 0.1, J = −0.25, and the default plan 73% stocks / 4% gold / 23% bonds with g = 0.0746824 (JS `erf` is an approximation, so agreement is to about 6 digits). Uploaded CSVs are parsed client-side only, so never add network calls. There is no JS test runner; verify by opening the page in headless Edge with an isolated `--user-data-dir` and an injected test script.

`investor.html` is one minimalist page, no tabs. At first only the inputs card shows (no heading): money, yearly saving, interest rate, inflation, market buttons, risk slider, **your guess** (stocks/gold/bonds sliders, savings = the rest; `guessW()`), everything else in a closed "More settings" `<details>`, and a **Calculate my plan** button. Inputs only run `onInputChange()` → `renderInputs()` (live, never reads `cur`); the button runs `calculate()`, which builds everything once into `cur` (`compute()`, including snapshots of the uploaded `data` and the surprise switches) and reveals `#results`. Changing an input afterwards dims the results (`.stale`) until the button is pressed again. Results: (2) the plan (`bestMix` for today's market, donut, tiles, "Your guess vs the math" via `mixStats`, crash-safe callout); (3) the plan in real life: main strategies are your plan, your guess and the savings account (`main: true` in `lifeStrategies`), extras (crash-safe, stocks + savings only, all stocks, all gold) only in the closed "Compare more strategies" table; `renderLifeMany` runs the 500 seeds in chunks (returns a promise, `life.many`) and shows who ended ahead, plan or guess; (4) "How we found your plan", closed `<details>` with the math (mix bars, bell curve, IVT + bisection, Module 4b crash card, history check). Resizing redraws charts only when the `innerWidth < 560` switch in `drawChart` flips, without recomputing.
- **`mixStats(w, m, z)`:** bell-curve growth and VaR of any split, built from `fundOf` + `growthRate` like `fundFraction`, so `mixStats(plan.w)` equals `plan.g` and `plan.var` exactly (tested).
- **`marketStates(p, lam, J)`:** builds the four states. The typed interest rate and inflation belong to the picked market (`p.mood`); other states shift them by the `MOODS` differences. The stock μ, σ are for a normal market and are shifted the same way. Picking a market button fills in typical rate/inflation (`pickMood`). Crashes use the rate and size from More settings (`crashSettings`), and the stock drift is compensated so the average return stays μ. With markets not changing, no crashes and no events, the years are the user's GBM plus savings account exactly (tested). `mathState(s)` gives the bell-curve numbers (full drift `mucAvg`) used by `bestMix`.
- **Crash-safe plan:** `bestMixCrash` (JS only) uses `riskEstimate` (one possible crash year) on a 2% grid; growth includes the crash. Shown only when crashes are on.
- **`set.events`:** with `events` off, the life events are drawn but ignored, so seeds stay comparable. Strikes (years over the limit) are counted; nobody stops.
- **Moods:** the `MOODS` values are templates (shifts for the stock, rate and inflation, absolute values for gold and bonds). They are illustrative, simulated simple yearly rates for universal assets (stocks, gold, bonds, savings account), roughly like long-run developed-market history, not data or any one country. They are converted to continuous rates (`rc`, `muc` = ln(1 + x)) before use.
- **`portfolioVar`:** the multi-asset generalisation of `value_at_risk`. With only stocks it must equal `valueAtRisk`.
- **Markets change over time:** each year's state (`regime`) follows the illustrative Markov chain `TRANSITION` over `REGIME_KEYS`. It is drawn in `generateYears` from the seeded rng.
- **Fair comparisons:** every strategy goes through `stepYear` on the same pre-generated, seeded `years`. Keep it that way.
- **Animations:** they must respect `prefers-reduced-motion`. Tests set `window.INVESTOR_NO_ANIM = true`.

The numbers quoted in `README.md` and `report/report.md` (e.g. f* = 1.75, f_VaR = 0.705014, C_BS = 10.4506) come from `python main.py` output; update both docs if parameters change.
