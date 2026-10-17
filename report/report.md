# Limits and Continuity in Finance: Choosing a Safe Investment Allocation

## 1. Problem

An investor has savings and two options: a risk-free account paying r = 3% per
year and a stock index with expected return μ = 10% and volatility σ = 20%.
**What fraction f of the savings should go into the stock?**

Too little and the money grows slowly. Too much and a bad year can wipe out a
large part of the savings. We want the **fastest long-run growth** subject to a
**risk rule**: *with 95% confidence, do not lose more than 15% in one year.*

Real investing happens in discrete steps (interest paid monthly, prices quoted
every second, portfolios rebalanced daily). Discrete models are hard to analyse,
so we take **limits** to get continuous models with clean formulas, and then
use **continuity** to guarantee and compute the answer. The project builds this
up in five modules, plus an extension with market crashes (Module 4b).

All rates are continuously compounded, as in Module 1: μ = 10% means an average
gain of $e^{0.1}-1 = 10.5\%$ a year, and r = 3% means 3.05% a year.

---

## 2. Module 1 — Continuous compounding

**Model.** P invested at annual rate r, compounded n times per year for t years:

$$A_n = P\left(1+\frac{r}{n}\right)^{nt}.$$

**Limit.** Write $A_n = P\exp\big(nt\ln(1+r/n)\big)$. Using the Taylor series
$\ln(1+x) = x - x^2/2 + O(x^3)$ with $x = r/n$:

$$nt\ln\left(1+\frac{r}{n}\right) = rt - \frac{r^2t}{2n} + O\!\left(\frac{1}{n^2}\right) \xrightarrow[n\to\infty]{} rt.$$

Because exp is **continuous**, the limit can move inside it:

$$\lim_{n\to\infty} A_n = P\exp\Big(\lim_{n\to\infty} nt\ln(1+r/n)\Big) = Pe^{rt}.$$

**Rate of convergence.** From the same expansion,
$A - A_n \approx Pe^{rt}\,\dfrac{r^2t}{2n}$, so the gap shrinks like $1/n$.

**Result** (P = 10,000, r = 10%, t = 10): yearly 25,937.42, monthly 27,070.41,
daily 27,179.10, limit 27,182.82. Going from monthly to daily compounding earns
about 109 more; going from daily to continuous earns only 3.72. The estimate
$Pe^{rt}r^2t/(2n)$ matches the real gap almost perfectly (Figure 1).

---

## 3. Module 2 — Binomial model → Black–Scholes

**Model.** Split [0, T] into N steps of length Δt = T/N. In each step the stock
moves up by $u = e^{\sigma\sqrt{\Delta t}}$ or down by $d = 1/u$. With the
risk-neutral probability $p = (e^{r\Delta t}-d)/(u-d)$, an option is priced by
working back through the tree:
$V = e^{-r\Delta t}\,[\,pV_{up} + (1-p)V_{down}\,]$.

**Limit.** As N → ∞ the number of up moves is binomial and, by the central
limit theorem, $\ln S_T$ becomes normally distributed. The binomial price
converges to the Black–Scholes formula

$$C_{BS} = S\,\Phi(d_1) - Ke^{-rT}\Phi(d_2),\qquad d_{1,2} = \frac{\ln(S/K) + (r\pm\sigma^2/2)T}{\sigma\sqrt{T}}.$$

**Result** (S = K = 100, r = 5%, σ = 20%, T = 1): $C_{BS} = 10.4506$.
N = 10 gives 10.2534, N = 100 gives 10.4306, N = 1000 gives 10.4486. The error
is about 2/N, so it goes to 0 like $1/N$ (Figure 2).

---

## 4. Module 3 — Random walk → Geometric Brownian Motion

**Model.** In each step of length Δt the log-price moves by $\pm\sigma\sqrt{\Delta t}$,
going up with probability
$p = \tfrac12 + \tfrac12\,\frac{\mu-\sigma^2/2}{\sigma}\sqrt{\Delta t}$.
This choice gives the walk the right mean drift.

**Limit.** After n = T/Δt steps, $\ln(S_T/S_0) = \sigma\sqrt{\Delta t}\,(2K-n)$
with $K\sim\text{Binomial}(n,p)$. Its mean is exactly $(\mu-\sigma^2/2)T$ and its
variance tends to $\sigma^2T$. As Δt → 0 the distribution converges to

$$\ln(S_T/S_0) \sim N\big((\mu-\sigma^2/2)T,\ \sigma^2T\big),$$

which is Geometric Brownian Motion, the continuous model used in Module 4.

**Measuring the convergence.** The Kolmogorov–Smirnov distance
$D_n = \sup_x |F_n(x) - F(x)|$ is computed **exactly** from the binomial
probabilities (no simulation noise): $D_1 = 0.426$, $D_{12} = 0.119$,
$D_{252} = 0.025$, $D_{25200} = 0.0025$. On a log-log plot it falls on a line
with slope −1/2, i.e. $D_n \propto 1/\sqrt{n}$ (Figure 3).

---

## 5. Module 4 — Risk optimization (main result)

### 5.1 Continuous-time portfolio

Keep a fraction f in the stock (GBM) and 1 − f in the risk-free asset, and
rebalance continuously (the Δt → 0 limit). Then

$$W_T = W_0 \exp\big(g(f)\,T + f\sigma\sqrt{T}\,Z\big),\qquad Z\sim N(0,1),$$
$$g(f) = r + f(\mu - r) - \tfrac12 f^2\sigma^2 .$$

### 5.2 Maximum growth: the derivative as a limit

g is a polynomial, so it is continuous and differentiable. Its derivative is the limit

$$g'(f) = \lim_{h\to 0}\frac{g(f+h)-g(f)}{h} = (\mu - r) - f\sigma^2 .$$

Numerically, at f = 1 the difference quotient gives 0.0280 (h = 0.1), 0.0298
(h = 0.01), 0.02998 (h = 0.001), … → **0.03** = exact value. The quotient's error
is $-h\sigma^2/2$, which goes to 0 with h.

Since $g''(f) = -\sigma^2 < 0$, g is concave and the maximum is where g'(f) = 0:

$$f^* = \frac{\mu - r}{\sigma^2} = \frac{0.07}{0.04} = 1.75 .$$

This is the **Kelly / Merton fraction**. Bisection on g' gives the same value
numerically. But f* = 1.75 means borrowing 75% extra, and its 95% one-year VaR
is **38.4%**, far above the 15% rule.

### 5.3 The risk constraint: continuity and the Intermediate Value Theorem

The 95% Value-at-Risk is the loss exceeded only 5% of the time:

$$\text{VaR}(f) = 1 - \exp\big(g(f)\,T - z_{0.95}\,f\sigma\sqrt{T}\big),\qquad z_{0.95} = 1.6449 .$$

Define $h(f) = \text{VaR}(f) - L$ with L = 0.15.

* h is **continuous** on [0, 2], because it is built from polynomials and exp.
* $h(0) = (1 - e^{0.03}) - 0.15 = -0.180 < 0$ (the risk-free portfolio never loses).
* $h(2) = 0.283 > 0$.

By the **Intermediate Value Theorem** there is an $f_{VaR}\in(0,2)$ with
$h(f_{VaR}) = 0$, i.e. VaR exactly equal to the limit.

**The root is unique, for any horizon.** $\text{VaR}(f) = L$ means
$\ln(1-L) = g(f)T - zf\sigma\sqrt{T}$, which is a quadratic equation:

$$\tfrac12\sigma^2T\,f^2 + \big(z\sigma\sqrt{T} - (\mu - r)T\big)f + \big(\ln(1-L) - rT\big) = 0 .$$

Its last coefficient is negative exactly when $h(0) < 0$, so the product of its two
roots, $(\ln(1-L) - rT)/(\tfrac12\sigma^2T)$, is negative: one root is negative and
one is positive. So $f_{VaR}$ is the only crossing with $f \ge 0$, and every
$f \le f_{VaR}$ satisfies the rule. (For T = 1, VaR is even increasing, because
$(\mu - r)\sqrt{T} < z\sigma$, here 0.07 < 0.329. For horizons beyond 22 years it first
falls, but the root is still unique.) The quadratic formula gives 0.705014 as well,
an independent check on the bisection.

**Bisection** turns the IVT into an algorithm. Each step halves the interval
and keeps the half where h changes sign. After k steps the root lies in an
interval of length $2/2^k \to 0$:

| k | a | b | mid | h(mid) |
|---|---|---|---|---|
| 1 | 0 | 2 | 1 | +0.0704 |
| 2 | 0 | 1 | 0.5 | −0.0508 |
| 3 | 0.5 | 1 | 0.75 | +0.0109 |
| 4 | 0.5 | 0.75 | 0.625 | −0.0196 |
| … | | | | |
| 41 | 0.705014 | 0.705014 | **0.705014** | 1e−13 |

The normal quantile $z_{0.95}$ itself is found the same way: Φ is continuous
and increasing, so Φ(z) − 0.95 has a root by the IVT.

### 5.4 Decision

Because g is increasing on [0, f*], the best feasible choice is the largest
allowed f:

$$f_{opt} = \min(f^*, f_{VaR}, f_{max}) = \min(1.75,\ 0.705,\ 2) = 0.705 .$$

**Invest 70.5% in the stock and 29.5% in the risk-free asset.** The typical
(median) growth is 6.94% per year (compared with 9.13% at Kelly), and the VaR is
exactly 15%. The *average* wealth grows faster, at $r + f(\mu - r) = 7.94\%$, but g is
what a typical investor gets in the long run: $g = E[\ln(W_T/W_0)]/T$.

### 5.5 Verification by simulation

Simulating 20,000 one-year paths with discrete rebalancing gives VaR ≈ 14.6%
(yearly), 15.1% (quarterly), 14.8% (monthly), 14.9% (weekly), 14.9% (daily).
These are all within Monte Carlo noise of the analytic 15%. Figure 4 shows the
growth curve, the VaR curve crossing the limit, the wealth bands (much narrower
for $f_{opt}$ than for Kelly), and the simulated loss distribution.

### 5.6 How bad are the bad years? Expected Shortfall

VaR says where the worst 5% of years *start*. It does not say how bad they are. The
**Expected Shortfall** is the average loss in those years. With
$\ln(W_T/W_0)\sim N(m, s^2)$, $m = g(f)T$, $s = f\sigma\sqrt T$, completing the square
(as for $E[e^X]$) gives

$$\text{ES}(f) = 1 - e^{m + s^2/2}\,\frac{\Phi(-z_\alpha - s)}{1-\alpha} \;\ge\; \text{VaR}(f).$$

At $f_{VaR} = 0.705$ the worst 5% of years lose **19.76% on average**, although VaR
only says "more than 15%". A simulation with daily rebalancing gives 19.91%.

ES is continuous in f, so the same IVT + bisection works for the stricter rule
"lose at most 15% on average in the worst 5%": $f_{ES} = 0.547$. This root is unique
too: $\ln(1 - \text{ES}(f)) = (r + f(\mu - r))T + \ln\Phi(-z - f\sigma\sqrt T) - \ln(1-\alpha)$
is concave because $\Phi$ is log-concave.

### 5.7 Fat tails: adding rare crashes

The bell curve makes big losses almost impossible, but real markets sometimes
crash. Following Merton's jump-diffusion model, the stock now also falls by
J = −25% on rare days. These crashes arrive at an average rate of λ = 0.1 per
year, i.e. about once every 10 years. The drift is lowered by λJ, so the
stock's **average return is still μ**: only the shape of the distribution
changes, with a fatter left tail.

**A new limit: crash coin-flips → Poisson.** Split the year into n steps and
let a crash happen in each step with probability λ/n. The number of crashes is
Binomial(n, λ/n). As n → ∞ it converges to the Poisson law,
$P(N=k) = e^{-\lambda}\lambda^k/k!$. The largest difference between the two
probability tables is $9.5\cdot10^{-3}$ for n = 1, $3.4\cdot10^{-5}$ for n = 252
and $8.6\cdot10^{-7}$ for n = 10,000.

**Wealth.** With continuous rebalancing, a crash multiplies the portfolio by
$(1+fJ)$, so with $N \sim \text{Poisson}(\lambda T)$ crashes

$$\ln\frac{W_T}{W_0} = \big(g(f) - \lambda J f\big)T + N\ln(1+fJ) + f\sigma\sqrt{T}\,Z .$$

Its CDF is a **mixture of bell curves**, one for each number of crashes k:

$$F(x) = \sum_{k\ge 0} P(N=k)\,\Phi\!\left(\frac{x - (g(f)-\lambda J f)T - k\ln(1+fJ)}{f\sigma\sqrt{T}}\right).$$

**Two uses of the IVT.**

1. The mixture has no inverse formula, but F is continuous and increasing from
   0 to 1. So F(x) = 0.05 has exactly one solution, and bisection finds it.
   This gives $\text{VaR}_{crash}(f) = 1 - e^{x}$.
2. $h(f) = \text{VaR}_{crash}(f) - L$ is continuous, with $h(0) = -0.180 < 0$
   and $h(2) = 0.361 > 0$. So again a root exists, and bisection finds it in
   41 steps:

| k | a | b | mid | h(mid) |
|---|---|---|---|---|
| 1 | 0 | 2 | 1 | +0.1060 |
| 2 | 0 | 1 | 0.5 | −0.0335 |
| 3 | 0.5 | 1 | 0.75 | +0.0373 |
| 4 | 0.5 | 0.75 | 0.625 | +0.0021 |
| … | | | | |
| 41 | 0.617447 | 0.617447 | **0.617447** | 6e−14 |

With λ = 0 every formula reduces exactly to Section 5.3.

**Same average return, lower growth.** The long-run growth with crashes is

$$g_{crash}(f) = \frac{E[\ln(W_T/W_0)]}{T} = g(f) - \lambda Jf + \lambda\ln(1 + fJ) \;<\; g(f),$$

because $\ln(1 + x) < x$: a big loss hurts more than the extra drift of the good
years helps. Its slope $(\mu - r) - f\sigma^2 - \lambda J + \lambda J/(1 + fJ)$ is
continuous, equals $\mu - r > 0$ at f = 0 and goes to $-\infty$ where one crash would
wipe out the portfolio, so by the IVT bisection finds the **crash-aware Kelly
fraction $f^*_{crash} = 1.41$** (instead of 1.75). The decision is again
$f_{opt} = \min(f^*_{crash}, f_{crash}, f_{max})$.

**Result.** Crashes matter most for the biggest losses:

| at f = 0.705 | bell curve only | with rare crashes |
|---|---|---|
| 95% VaR | 15.00% | 17.47% |
| 99% VaR | 22.79% | 27.29% |
| average loss in the worst 5% (ES) | 19.76% | 23.48% |
| years losing more than 15% | 1 in 20 | 1 in 14 |
| P(lose more than 30%) | 0.1% | 0.6% |
| typical growth per year | 6.94% | 6.76% |

The 70.5% answer therefore breaks the rule in a world with crashes. The
crash-aware answer is **$f_{opt} = \min(1.41,\ 0.617,\ 2) = f_{crash}$ = 61.7% in stocks**. A Monte Carlo check
with crash coin-flips at every step gives a VaR of 14.5% (yearly steps), 15.0%
(monthly), 15.0% (weekly) and 15.0% (daily) at $f_{crash}$, converging to 15%.
This is the binomial → Poisson limit at work. See Figure 5.

### 5.8 More assets: stocks, gold, bonds and the bank

So far the only choice was stocks or the bank. Now the investor can also hold
gold and bonds. The numbers are illustrative long-run values, not data:

| asset | μ (continuous) | σ |
|---|---|---|
| stocks | 10% | 20% |
| gold | ln 1.05 = 4.88% | 15% |
| bonds | ln 1.045 = 4.40% | 6% |

The correlations are 0 (stocks–gold), −0.2 (stocks–bonds) and 0.3 (gold–bonds).
With weights $w$ in the three risky assets, $1 - \sum w_i$ in the bank,
no borrowing ($\sum w_i \le 1$) and no short selling ($w_i \ge 0$), continuous
rebalancing gives a normal log-wealth with

$$g(w) = r + w\cdot(\mu - r) - \tfrac12 w^\top C w, \qquad \mathrm{sd}(w) = \sqrt{w^\top C w},$$

where $C_{ij} = \rho_{ij}\sigma_i\sigma_j$. The VaR has the same form as in
Section 5.3, with $f\sigma$ replaced by $\mathrm{sd}(w)$.

**One mix = one fund.** Fix a mix $d$ of the risky assets ($d_i \ge 0$,
$\sum d_i = 1$) and put a share $f$ of the money in it, so $w = f d$. Then

$$g(fd) = r + f(\mu_d - r) - \tfrac12 f^2\sigma_d^2, \qquad \mu_d = r + d\cdot(\mu - r), \quad \sigma_d = \sqrt{d^\top C d}.$$

This is exactly Module 4 with one "stock" $(\mu_d, \sigma_d)$. So for every
fund, VaR is continuous in $f$, the IVT guarantees the root $f_{VaR}(d)$, and
bisection finds it; the fund's best share is $\min(f^*_d, f_{VaR}(d), 1)$.
Trying every mix on a 1% grid (5,151 funds, so 5,151 bisections) and keeping
the fastest-growing one gives the best mix.

**Result.** With the 15% / 95% rule:

| choice | allocation | growth g |
|---|---|---|
| bonds + bank | 100% bonds | 4.22% |
| gold + bank | 81.2% gold | 3.78% |
| stocks + bank (Module 4) | 70.5% stocks | 6.94% |
| **best mix** | **73% stocks, 4% gold, 23% bonds** | **7.47%** |

The best mix uses all the money (nothing in the bank) and has a VaR of 15.0%,
so it gets **0.53% more growth per year for the same risk**. Bonds replace the
bank because they pay a bit more and tend to move against stocks. A Monte Carlo
check with daily rebalancing gives a VaR of 14.95%. See Figure 6. With gold and
bonds paying only the bank rate, the method returns exactly the Module 4 answer
(70.5% stocks), which is tested.

**Exactly one best mix.** The problem is convex: $g$ is concave, and
$\mathrm{sd}(w) = \sqrt{w^\top Cw}$ is a norm, so $g(w)T - z\,\mathrm{sd}(w)\sqrt T$ is concave
and the mixes with VaR ≤ L form a convex set (a superlevel set of a concave function).
A strictly concave function has at most one maximum on a convex set, so the best mix is
unique. Its top is flat, though: every mix on the 1% grid within 0.01% a year of the best
has 73% stocks, but anywhere from 1% to 7% gold (bonds 20–26%). **About 73% stocks is
the message**; the split between gold and bonds barely matters.

**Module 5's own limit: grid step → 0.** Let $G_s$ be the best growth on the grid with
step 1/s. If s divides s′, the finer grid contains the coarser one, so
$G_s \le G_{s'}$: the sequence $G_{25}, G_{50}, G_{100}, \dots$ is increasing and
bounded by the true maximum, so it converges (monotone convergence). The limit *is* the
true maximum, because each fund's growth depends continuously on the mix and the grids
fill the triangle. For the bell curve, $f_{VaR}$ of every fund is also given by the
quadratic formula of Section 5.3, so even grids with a million mixes take a second:

| grid step | mixes | best mix (stocks / gold / bonds) | growth g |
|---|---|---|---|
| 4% | 351 | 72 / 12 / 16 | 7.45461% |
| 2% | 1,326 | 72 / 12 / 16 | 7.45461% |
| 1% | 5,151 | 73 / 4 / 23 | 7.46824% |
| 0.5% | 20,301 | 72.98 / 4.50 / 22.49 | 7.46828% |
| 0.25% | 80,601 | 73 / 4.25 / 22.75 | 7.46875% |
| 0.125% | 321,201 | 73 / 4.25 / 22.75 | 7.46875% |
| 0.0625% | 1,282,401 | 72.94 / 5 / 22.06 | 7.46882% |
| → 0 (zooming in) | | 72.9 / 4.9 / 22.2 | **7.46887%** |

The last row zooms in on the best fund with a small grid whose step halves 40 times. The
1% grid is only 0.0006% a year short of the exact best, but its gold share (4%) is not
precise: the exact best has 4.9%.

---

## 6. Conclusion

| Concept | Where it is used |
|---|---|
| Limit of a sequence | $(1+r/n)^{nt} \to e^{rt}$; binomial price → Black–Scholes |
| Limit of distributions | random walk → GBM as Δt → 0, measured by the KS distance |
| Derivative as a limit | g'(f) from difference quotients → Kelly fraction |
| Continuity + IVT | existence of $f_{VaR}$, $f_{ES}$ and $z_{0.95}$; bisection computes them |
| Concavity | uniqueness of $f_{VaR}$ and $f_{ES}$ for any horizon; exactly one best mix |
| Poisson limit + IVT twice | crash coin-flips → Poisson; quantile of a mixture, $f_{crash}$ and $f^*_{crash}$ by bisection |
| IVT for every fund | best mix of stocks, gold and bonds: one bisection per mix |
| Monotone convergence | finer and finer grids of mixes → the exact best mix |

Limits turn messy discrete models into continuous formulas, and continuity
guarantees the answer exists, so it can be found reliably. The result is a
concrete recommendation: **70.5% in stocks** under a "max 15% loss at 95%
confidence" rule. If rare crashes are taken into account, the same method gives
**61.7%**. With gold and bonds available, the best mix is about **73% stocks**, with
the rest in bonds and a little gold (4% on the 1% grid, 4.9% exactly), which grows
0.53% a year faster for the same risk.

**Limitations.** Constant μ and σ, no transaction costs, and borrowing at r. Fat
tails are handled only by the simple crash model of Section 5.7, with an
assumed crash size and frequency. VaR ignores how bad the worst years are; the
Expected Shortfall of Section 5.6 shows it. Extensions: estimate μ and σ (and the crash
rate) from real price data, or use historical or Student-t VaR. The gold and
bond numbers of Section 5.7 are assumptions, and real average returns are very
hard to estimate, so the best mix can change a lot with small changes in the
inputs.

The real-life test in `web/investor.html` breaks these assumptions on purpose:
markets change from year to year and crash years happen. There the plan can
break its "1 year in 20" promise more often (about 1 year in 13 with the default
surprises) and can lose to other strategies in some futures. This does not
contradict the result above; it shows what the assumptions leave out. The
simple rule is a starting point, and the test shows where it stops being enough.

## References
1. Hull, J. C. *Options, Futures, and Other Derivatives.*
2. Cox, J., Ross, S., Rubinstein, M. (1979). Option pricing: a simplified approach. *J. Financial Economics* 7, 229–263.
3. Merton, R. C. (1969). Lifetime portfolio selection under uncertainty. *Rev. Economics and Statistics* 51, 247–257.
4. Kelly, J. L. (1956). A new interpretation of information rate. *Bell System Technical Journal* 35, 917–926.
5. Stewart, J. *Calculus* — Intermediate Value Theorem and limits.
6. Merton, R. C. (1976). Option pricing when underlying stock returns are discontinuous. *J. Financial Economics* 3, 125–144.
