function p = jump_log_wealth_cdf(x, f, mu, sigma, r, T, lam, J, k_max)
%JUMP_LOG_WEALTH_CDF P(log(W_T / W_0) <= x) when rare crashes can happen.
%   A crash multiplies the stock by (1 + J) and happens on average lam times a
%   year; the drift is compensated (mu - lam J) so the average return is still mu.
%   With N ~ Poisson(lam T) crashes,
%       log(W_T / W_0) = (g(f) - lam J f) T + N ln(1 + f J) + f sigma sqrt(T) Z,
%   so the CDF is a mixture of bell curves, one for each number of crashes k.
if nargin < 9, k_max = 30; end
m = (growth_rate(f, mu, sigma, r) - lam * J * f) * T;
sd = abs(f) * sigma * sqrt(T);
jump = log(1 + f * J);
k = 0:k_max;
weights = arrayfun(@(i) poisson_pmf(i, lam * T), k);
p = sum(weights .* norm_cdf((x - m - k * jump) / sd));
end
