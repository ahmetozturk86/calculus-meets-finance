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

- `src/bisection.py` is the shared core: `bisection(fun, a, b)` returns `(root, history)` and raises `ValueError` when the IVT sign condition fails. It is used by `normal.norm_ppf` (inverse normal CDF), `risk.numerical_kelly_fraction` (root of g'), and `risk.max_fraction_under_var` (root of VaR(f) − L). The history rows `[a, b, mid, fun(mid)]` are printed in `main.py` as the IVT demonstration.
- `src/normal.py` implements the normal CDF/PDF/PPF with `math.erf` on purpose — no scipy — so code ports to base MATLAB without the Statistics Toolbox.
- Module chain: `compounding.py` (1) → `option_pricing.py` (2) → `gbm.py` (3) → `risk.py` (4). Module 4's closed-form wealth/VaR formulas assume the continuous limits established in modules 1 and 3.
- `gbm.ks_distance_to_gbm` computes the KS distance **exactly** from binomial probabilities (no Monte Carlo), so its convergence plot is noise-free.
- `main.py` holds all model parameters as constants at the top (savings, option, market, risk-limit) and one function per module that prints tables and saves one figure.

## Python ↔ MATLAB parity

`matlab/` mirrors `src/` one function per `.m` file with the **same names and argument order**; `matlab/main.m` mirrors `main.py` and `matlab/run_tests.m` mirrors `tests/test_modules.py`. When changing a Python function, make the same change in its `.m` file (and vice versa). Differences to keep in mind:
- Python random functions take an `rng` argument; MATLAB versions drop it and use the global stream (`rng(seed)`), so Monte Carlo numbers differ between languages while deterministic numbers must match.
- `quantile_linear.m` reimplements `numpy.quantile` (linear interpolation) to avoid the Statistics Toolbox.
- Keep Python code MATLAB-portable: plain functions, numpy arrays only, no classes/pandas/scipy.

The numbers quoted in `README.md` and `report/report.md` (e.g. f* = 1.75, f_VaR = 0.705014, C_BS = 10.4506) come from `python main.py` output; update both docs if parameters change.
