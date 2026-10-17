# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

School project + GitHub portfolio piece: four finance modules, each showing a discrete model converging to a continuous one (limits), ending in a VaR-constrained portfolio allocation found via the Intermediate Value Theorem + bisection. The school requires **MATLAB**; Python is the reference implementation.

## Commands

```bash
pip install -r requirements.txt           # numpy, matplotlib, pytest
python main.py                            # prints all result tables, writes figures/*.png
python -m pytest                          # all tests (pytest.ini sets pythonpath=. and testpaths=tests)
python -m pytest tests/test_modules.py::test_var_root_hits_the_limit   # single test
```

MATLAB (from the `matlab/` folder): `run_tests` then `main` (figures go to `matlab/figures/`). MATLAB is not installed on the dev machine, so MATLAB changes cannot be run locally — review them carefully.

## Architecture

- `src/bisection.py` is the shared core: `bisection(fun, a, b)` returns `(root, history)` and raises `ValueError` when the IVT sign condition fails. It is used by `normal.norm_ppf` (inverse normal CDF), `risk.numerical_kelly_fraction` (root of g'), `risk.max_fraction_under_var` (root of VaR(f) − L), `risk.value_at_risk_jump` (quantile of the crash mixture) and `risk.max_fraction_under_var_jump` (root of VaR_jump(f) − L, so a bisection whose function is itself a bisection). The history rows `[a, b, mid, fun(mid)]` are printed in `main.py` as the IVT demonstration.
- `src/normal.py` implements the normal CDF/PDF/PPF with `math.erf` on purpose — no scipy — so code ports to base MATLAB without the Statistics Toolbox.
- Module chain: `compounding.py` (1) → `option_pricing.py` (2) → `gbm.py` (3) → `risk.py` (4). Module 4's closed-form wealth/VaR formulas assume the continuous limits established in modules 1 and 3. Module 4b (fat tails, `main.module4b_crashes`) adds Merton-style crashes: size J, rate λ, drift compensated by λJ so the average return stays μ. With λ = 0 every crash function must equal its Module 4 twin (tested).
- `gbm.ks_distance_to_gbm` computes the KS distance **exactly** from binomial probabilities (no Monte Carlo), so its convergence plot is noise-free.
- `main.py` holds all model parameters as constants at the top (savings, option, market, risk-limit) and one function per module that prints tables and saves one figure.

## Python ↔ MATLAB parity

`matlab/` mirrors `src/` one function per `.m` file with the **same names and argument order**; `matlab/main.m` mirrors `main.py` and `matlab/run_tests.m` mirrors `tests/test_modules.py`. When changing a Python function, make the same change in its `.m` file (and vice versa). Differences to keep in mind:
- Python random functions take an `rng` argument; MATLAB versions drop it and use the global stream (`rng(seed)`), so Monte Carlo numbers differ between languages while deterministic numbers must match.
- `quantile_linear.m` reimplements `numpy.quantile` (linear interpolation) to avoid the Statistics Toolbox.
- Keep Python code MATLAB-portable: plain functions, numpy arrays only, no classes/pandas/scipy.

## Web pages

`web/index.html` (explainer, uses Chart.js from a CDN) and `web/investor.html` (real-life tool + game, no external libraries, SVG charts) re-implement the deterministic formulas in JavaScript: `bisection`, `normPpf`, `growthRate`, `kellyFraction`, `valueAtRisk`, `poissonPmf`, `jumpLogWealthCdf`, `valueAtRiskJump`, `maxFractionUnderVar` and `quantileLinear` mirror `src/bisection.py`, `src/normal.py` and `src/risk.py`. Keep them in sync when the Python changes. The JS numbers must match Python, for example f_VaR = 0.705014 with the default parameters and f_max = 2, and f_crash = 0.617447 with λ = 0.1, J = −0.25 (JS `erf` is an approximation, so agreement is to about 6 digits). Uploaded CSVs are parsed client-side only, so never add network calls. There is no JS test runner; verify by opening the page in headless Edge with an isolated `--user-data-dir` and an injected test script.

Tab 1 of `investor.html` holds every input. Tab 2 plays the answer through simulated years (`lifeScenario`, `drawLife`, `renderLifeMany` over 500 seeds), and tab 4 is the game. Both use the same engine and settings (`lifeSettings`):
- **`marketStates(p, lam, J)`:** builds the four states from tab 1. "Normal" is exactly the user's μ, σ, r. Calm, crisis and high inflation are shifted by the `MOODS` differences. Crashes use the tab 1 rate and size (`crashSettings`), and the stock drift is compensated so the average return stays μ. With markets not changing, no crashes and no events, the years are the user's GBM plus bank exactly (tested).
- **`set.fire` / `set.events`:** tab 2 sets `fire: false`, so strikes are counted but nobody stops. With `events` off, the life events are drawn but ignored, so seeds stay comparable.
- **Moods:** the `MOODS` values are templates (shifts for the stock and bank, absolute values for gold, bonds and inflation). They are illustrative, simulated simple yearly rates for universal assets (stocks, gold, bonds, savings account), roughly like long-run developed-market history, not data or any one country. They are converted to continuous rates (`rc`, `muc` = ln(1 + x)) before use.
- **`portfolioVar`:** the multi-asset generalisation of `value_at_risk`. With only stocks it must equal `valueAtRisk`.
- **Markets change over time:** each year's state (`regime`) follows the illustrative Markov chain `TRANSITION` over `REGIME_KEYS`. It is drawn in `generateYears` from the seeded rng.
- **`riskEstimate(w, mix, alpha)`:** takes either one mood or a list of [mood, probability]. It mixes the possible states, each with or without a market crash (`mood.crash`), and finds the quantile of the mixture by bisection.
- **`forecast(g)`:** what the player knows: last year's state plus its transition row.
- **`adaptiveSafe`:** the project's rule, re-applied each year for last year's state.
- **Fair comparisons:** every strategy goes through `stepYear` on the same pre-generated, seeded `years`, so comparisons and "replay" are fair. Keep it that way.
- **Animations:** they must respect `prefers-reduced-motion`. Tests set `window.INVESTOR_NO_ANIM = true`.

The numbers quoted in `README.md` and `report/report.md` (e.g. f* = 1.75, f_VaR = 0.705014, C_BS = 10.4506) come from `python main.py` output; update both docs if parameters change.
