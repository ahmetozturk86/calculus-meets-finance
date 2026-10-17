# Limits in Finance — from discrete models to continuous risk optimization

**How much of my savings should I put in the stock market?**
This project answers that question with calculus. Every module takes a discrete
financial model, lets a parameter go to a **limit** (n → ∞, Δt → 0), and shows
the result converges to a continuous model. The last module uses the
**continuity** of the risk function (Intermediate Value Theorem + bisection) to
pick the allocation with the highest growth that still respects a risk limit.

Implemented in **Python** (with tests) and ported 1:1 to **MATLAB**.

| # | Module | Limit / continuity idea | Real-life question |
|---|--------|------------------------|--------------------|
| 1 | Continuous compounding | lim (1 + r/n)ⁿᵗ = eʳᵗ | Is daily compounding worth more than monthly? |
| 2 | Binomial tree → Black–Scholes | lim_{N→∞} C_N = C_BS | What is a fair price for a stock option? |
| 3 | Random walk → Geometric Brownian Motion | distribution converges as Δt → 0 | Why do we model stock prices with GBM? |
| 4 | **Risk optimization** | g'(f) as a limit, IVT on VaR(f) − L | How much should I invest in stocks? |

## Results

### 1. Continuous compounding
10,000 at 10% for 10 years. More frequent compounding pays more, but the gain is
bounded: the gap to the limit `P e^{rt} = 27,182.82` shrinks like `1/n`.

| compounding | value | gap to limit |
|---|---|---|
| yearly | 25,937.42 | 1,245.39 |
| monthly | 27,070.41 | 112.40 |
| daily | 27,179.10 | 3.72 |
| every minute | 27,182.82 | 0.003 |

![compounding](figures/fig1_compounding.png)

### 2. Binomial option price → Black–Scholes
A call option with S = K = 100, r = 5%, σ = 20%, T = 1. The N-step binomial price
oscillates around the Black–Scholes price 10.4506 and the error goes to 0 like `1/N`.

![binomial](figures/fig2_binomial_to_black_scholes.png)

### 3. Random walk → GBM
A coin-flip walk with steps ±σ√Δt. As Δt → 0 the distribution of log-returns
converges to the normal law of GBM. The Kolmogorov–Smirnov distance (computed
exactly, not simulated) falls like `1/√n`.

![gbm](figures/fig3_random_walk_to_gbm.png)

### 4. Risk optimization
A stock with μ = 10%, σ = 20%, and a risk-free asset at r = 3%. Investing a
fraction f in the stock gives the growth rate

    g(f) = r + f(μ − r) − ½ f² σ²

* **Max growth (Kelly):** g is continuous and concave, so solving g'(f) = 0 gives
  f* = (μ − r)/σ² = **1.75** (175%, borrowed money). The 95% one-year VaR is then **38%**.
* **Risk rule:** "with 95% confidence, lose at most 15% in a year." The function
  h(f) = VaR(f) − 0.15 is continuous with h(0) < 0 < h(2), so by the
  **Intermediate Value Theorem** it has a root, which bisection finds in 41 steps:
  f_VaR = **0.705**.
* **Decision:** f_opt = min(f*, f_VaR) = **0.705**, so invest **70.5% in the stock and
  29.5% risk-free**. Expected growth is 6.94%/year with exactly 15% VaR. A Monte
  Carlo simulation with daily rebalancing confirms it (simulated VaR ≈ 14.9–15.1%).

![risk](figures/fig4_risk_optimization.png)

## Run it

### Python
```bash
pip install -r requirements.txt
python main.py          # prints all tables, writes figures/*.png
python -m pytest        # 15 tests
```

### MATLAB
```matlab
cd matlab
main          % same tables, figures saved to matlab/figures/
run_tests     % same checks as the Python tests
```
Only base MATLAB is needed (R2018b or newer): the normal CDF and quantiles are
written by hand, so no Statistics Toolbox is required. Deterministic numbers
match Python exactly. Monte Carlo numbers differ slightly because the random
generators are different.

### Web pages (no install, work offline)
Open these files in a browser:
- `web/index.html`: interactive explainer. Move a slider and watch each discrete model turn into its continuous limit.
- `web/investor.html`: **The Investor**, the real-life version of Module 4. You enter your savings, loss limit and
  market numbers by hand, or upload a CSV of prices (`web/sample_prices.csv` is a simulated example).
  It then:
  - estimates μ and σ from the prices;
  - computes how much to put in stocks with the same IVT + bisection method;
  - checks the answer against historical one-year returns (historical VaR and a backtest);
  - includes a game where you manage the money year by year and are compared with the safe and greedy strategies.

  The uploaded file is read only inside the page and is never sent anywhere. `investor.html` uses no external libraries.

## Project layout
```
src/                Python modules
  bisection.py      IVT-based root finding (used everywhere)
  normal.py         normal CDF / inverse CDF without scipy
  compounding.py    Module 1
  option_pricing.py Module 2
  gbm.py            Module 3
  risk.py           Module 4
main.py             runs everything, saves figures
tests/              pytest checks against known values
matlab/             one .m file per Python function + main.m, run_tests.m
report/report.md    full write-up with derivations
web/                interactive explainer (index.html) and The Investor tool and game (investor.html)
```

## Limitations
The model assumes constant μ and σ, normally distributed log-returns, free
continuous rebalancing, and borrowing at the risk-free rate. Real markets have
fat tails, so a normal VaR underestimates extreme losses. The point of the
project is the mathematics (limits and continuity), not investment advice.
