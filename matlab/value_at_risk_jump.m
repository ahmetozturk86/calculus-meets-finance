function v = value_at_risk_jump(f, mu, sigma, r, T, lam, J, alpha, k_max)
%VALUE_AT_RISK_JUMP VaR at level alpha when rare crashes can happen.
%   With lam = 0 it equals value_at_risk. The mixture CDF has no inverse
%   formula, but it is continuous and increasing, so by the IVT it crosses
%   1 - alpha exactly once, and bisection finds that quantile.
if nargin < 8, alpha = 0.95; end
if nargin < 9, k_max = 30; end
if f == 0
    v = 1 - exp(r * T);
    return
end
if 1 + f * J <= 0
    error('risk:wipeout', 'a crash would wipe out the portfolio: need 1 + f J > 0');
end
m = (growth_rate(f, mu, sigma, r) - lam * J * f) * T;
sd = abs(f) * sigma * sqrt(T);
lo = m + k_max * min(0, log(1 + f * J)) - 12 * sd;
hi = m + max(0, k_max * log(1 + f * J)) + 12 * sd;
x = bisection(@(x) jump_log_wealth_cdf(x, f, mu, sigma, r, T, lam, J, k_max) - (1 - alpha), ...
    lo, hi, 1e-12);
v = 1 - exp(x);
end
