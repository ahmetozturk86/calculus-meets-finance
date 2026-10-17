function es = expected_shortfall_jump(f, mu, sigma, r, T, lam, J, alpha, k_max)
%EXPECTED_SHORTFALL_JUMP Expected Shortfall when rare crashes can happen (lam = 0: EXPECTED_SHORTFALL).
%   Below the VaR point x_q = ln(1 - VaR_crash) (found by bisection), the bell curve of k crashes adds
%   e^(m_k + s^2/2) Phi((x_q - m_k - s^2) / s), with m_k = (g(f) - lam J f) T + k ln(1 + f J).
if nargin < 8, alpha = 0.95; end
if nargin < 9, k_max = 30; end
if f == 0
    es = 1 - exp(r * T);
    return
end
x_q = log(1 - value_at_risk_jump(f, mu, sigma, r, T, lam, J, alpha, k_max));
s = abs(f) * sigma * sqrt(T);
k = 0:k_max;
weights = arrayfun(@(i) poisson_pmf(i, lam * T), k);
m_k = (growth_rate(f, mu, sigma, r) - lam * J * f) * T + k * log(1 + f * J);
tail = sum(weights .* exp(m_k + 0.5 * s^2) .* norm_cdf((x_q - m_k - s^2) / s));
es = 1 - tail / (1 - alpha);
end
