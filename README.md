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
| 4b | Fat tails: rare crashes | Binomial(n, λ/n) → Poisson(λ), IVT twice | What if the market sometimes crashes? |
| 5 | More assets: stocks, gold, bonds, bank | one IVT + bisection per mix | What is the best mix of several assets? |

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

### 4b. Fat tails: rare crashes
The bell curve makes big losses almost impossible, but real markets sometimes crash.
Here the stock also falls **25% about once every 10 years** (a Merton jump-diffusion).
The average return stays 10%, so only the left tail gets fatter.

* **New limit:** crash coin-flips with probability λ/n in each of n steps become the
  Poisson law as n → ∞. The largest table difference falls from 9.5·10⁻³ (n = 1) to
  8.6·10⁻⁷ (n = 10,000).
* **IVT twice:** the loss distribution is now a mixture of bell curves with no inverse
  formula. Bisection finds its 5% quantile (the VaR), and a second bisection finds
  where VaR(f) hits the limit.

| at f = 0.705 | bell curve only | with rare crashes |
|---|---|---|
| 95% VaR | 15.00% | 17.47% |
| 99% VaR | 22.79% | 27.29% |
| years losing more than 15% | 1 in 20 | 1 in 14 |

With crashes the same rule gives **f_crash = 0.617, i.e. 61.7% in stocks** instead of
70.5%. A Monte Carlo check with daily crash coin-flips gives VaR = 15.0% at f_crash.

![fat tails](figures/fig5_fat_tails.png)

### 5. More assets: stocks, gold and bonds
Now gold (4.9% a year, σ = 15%) and bonds (4.4%, σ = 6%) are available too; the
numbers are illustrative. **Key idea:** any fixed mix of stocks, gold and bonds
behaves like one fund, so Module 4 (IVT + bisection) gives the best share for
that fund. Trying every mix in 1% steps (5,151 bisections) and keeping the
fastest-growing one gives the best mix, with no borrowing and no short selling.

| choice | allocation | growth per year |
|---|---|---|
| stocks + bank (Module 4) | 70.5% stocks | 6.94% |
| **best mix (Module 5)** | **73% stocks, 4% gold, 23% bonds** | **7.47%** |

Same 15% bad-year risk, **+0.53% growth a year**. Bonds replace the bank because
they pay a bit more and tend to move against stocks.

![best mix](figures/fig6_best_mix.png)

## Run it

### Python
```bash
pip install -r requirements.txt
python main.py          # prints all tables, writes figures/*.png
python -m pytest        # 24 tests
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
- `web/investor.html`: **The Investor**, the real-life version of Modules 4 and 5, on one simple page.
  1. **Your situation:** the money you have, the money you add each year, the savings account interest rate,
     inflation, how the market is right now (calm, normal, crisis or high inflation) and how much you can lose
     in a bad year. Everything else is in a closed **More settings** section: years, goal, the stock numbers,
     surprises, and an optional CSV of real prices (`web/sample_prices.csv` is a simulated example).
  2. **Your plan:** the best split between stocks, gold, bonds and the savings account (Module 5), with the money
     in each, the typical growth, the bad-year loss and the typical result after the chosen years. If rare
     crashes are switched on, it also shows a safer crash-aware mix.
  3. **Your plan in real life:** the plan is played year by year in a simulated future and compared, on exactly
     the same years, with stocks + savings only (Module 4), savings only, all stocks, an equal 25% mix and all gold.
     A table over 500 possible futures shows whether the plan keeps its promise ("a year over your limit at most
     1 in 20"). With no surprises the promise holds; with crises and crashes you see what the bell-curve math leaves out.
  4. **How we found your plan:** closed sections with the math: each asset alone vs the mix, the bell curve,
     the IVT + bisection, rare crashes (Module 4b) and the check against real history.

  **The simulated market.** The interest rate and inflation you type belong to the market you pick. The stock
  numbers are for a normal market, and calm, crisis and high-inflation markets are better or worse versions of
  them. With "markets change over time" on, they follow each other through a Markov chain, so crises tend to last
  a while. Gold, bonds and inflation use illustrative values, roughly like long-run history in developed markets.
  None of this is a forecast. If you uploaded prices, the stock years come from your real data.

  The uploaded file is read only inside the page and is never sent anywhere. `investor.html` uses no external libraries.

## Project layout
```
src/                Python modules
  bisection.py      IVT-based root finding (used everywhere)
  normal.py         normal CDF / inverse CDF without scipy
  compounding.py    Module 1
  option_pricing.py Module 2
  gbm.py            Module 3
  risk.py           Module 4 and 4b
  portfolio.py      Module 5 (more assets)
main.py             runs everything, saves figures
tests/              pytest checks against known values
matlab/             one .m file per Python function + main.m, run_tests.m
report/report.md    full write-up with derivations
web/                interactive explainer (index.html) and The Investor tool (investor.html)
```

## Limitations
The model assumes constant μ and σ, free continuous rebalancing, and borrowing at
the risk-free rate. Module 4 assumes normally distributed log-returns. Real
markets have fat tails, so a normal VaR underestimates extreme losses. Module 4b
adds rare crashes to fix part of this, but its crash size and frequency are
assumptions, not estimates. The point of the
project is the mathematics (limits and continuity), not investment advice.
The gold and bond numbers of Module 5 are assumptions too, and real average
returns are hard to estimate, so the best mix can change a lot with small input
changes. The real-life test in `investor.html` breaks the assumptions on purpose
(changing markets, crash years), so the plan can break its promise there: it
shows where the simple model stops being enough.
