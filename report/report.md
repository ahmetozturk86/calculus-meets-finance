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
up in four steps.

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
$h(f_{VaR}) = 0$, i.e. VaR exactly equal to the limit. Moreover $h'(f) > 0$
whenever $(\mu - r)\sqrt{T} < z\sigma$ (here 0.07 < 0.329), so VaR is increasing,
the root is unique, and every $f \le f_{VaR}$ satisfies the rule.

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

**Invest 70.5% in the stock and 29.5% in the risk-free asset.** Expected growth
is 6.94% per year (compared with 9.13% at Kelly), and the VaR is exactly 15%.

### 5.5 Verification by simulation

Simulating 20,000 one-year paths with discrete rebalancing gives VaR ≈ 14.6%
(yearly), 15.1% (quarterly), 14.8% (monthly), 14.9% (weekly), 14.9% (daily).
These are all within Monte Carlo noise of the analytic 15%. Figure 4 shows the
growth curve, the VaR curve crossing the limit, the wealth bands (much narrower
for $f_{opt}$ than for Kelly), and the simulated loss distribution.

### 5.6 Fat tails: adding rare crashes

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

**Result.** Crashes matter most for the biggest losses:

| at f = 0.705 | bell curve only | with rare crashes |
|---|---|---|
| 95% VaR | 15.00% | 17.47% |
| 99% VaR | 22.79% | 27.29% |
| years losing more than 15% | 1 in 20 | 1 in 14 |
| P(lose more than 30%) | 0.1% | 0.6% |

The 70.5% answer therefore breaks the rule in a world with crashes. The
crash-aware answer is **$f_{crash}$ = 61.7% in stocks**. A Monte Carlo check
with crash coin-flips at every step gives a VaR of 14.5% (yearly steps), 15.0%
(monthly), 15.0% (weekly) and 15.0% (daily) at $f_{crash}$, converging to 15%.
This is the binomial → Poisson limit at work. See Figure 5.

---

## 6. Conclusion

| Concept | Where it is used |
|---|---|
| Limit of a sequence | $(1+r/n)^{nt} \to e^{rt}$; binomial price → Black–Scholes |
| Limit of distributions | random walk → GBM as Δt → 0, measured by the KS distance |
| Derivative as a limit | g'(f) from difference quotients → Kelly fraction |
| Continuity + IVT | existence of $f_{VaR}$ and $z_{0.95}$; bisection computes them |
| Poisson limit + IVT twice | crash coin-flips → Poisson; quantile of a mixture and $f_{crash}$ by bisection |

Limits turn messy discrete models into continuous formulas, and continuity
guarantees the answer exists, so it can be found reliably. The result is a
concrete recommendation: **70.5% in stocks** under a "max 15% loss at 95%
confidence" rule. If rare crashes are taken into account, the same method gives
**61.7%**.

**Limitations.** Constant μ and σ, no transaction costs, and borrowing at r. Fat
tails are handled only by the simple crash model of Section 5.6, with an
assumed crash size and frequency. Extensions: estimate μ and σ (and the crash
rate) from real price data, use historical or Student-t VaR, or add several
stocks (an efficient frontier).

The game in `web/investor.html` breaks these assumptions on purpose: markets
change from year to year, crash years happen, and gold and bonds are available.
There the 70.5% rule can lose to other strategies. This does not contradict the
result above; it shows what the assumptions leave out. The simple rule is a
starting point, and the game shows where it stops being enough.

## References
1. Hull, J. C. *Options, Futures, and Other Derivatives.*
2. Cox, J., Ross, S., Rubinstein, M. (1979). Option pricing: a simplified approach. *J. Financial Economics* 7, 229–263.
3. Merton, R. C. (1969). Lifetime portfolio selection under uncertainty. *Rev. Economics and Statistics* 51, 247–257.
4. Kelly, J. L. (1956). A new interpretation of information rate. *Bell System Technical Journal* 35, 917–926.
5. Stewart, J. *Calculus* — Intermediate Value Theorem and limits.
6. Merton, R. C. (1976). Option pricing when underlying stock returns are discontinuous. *J. Financial Economics* 3, 125–144.
