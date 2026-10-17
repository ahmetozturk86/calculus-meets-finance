function [f, g, history] = fund_fraction(d, L, mu, sigma, rho, r, T, alpha, z)
%FUND_FRACTION Module 4 for the fund d: f = min(f*, f_VaR, 1), with f_VaR found by bisection.
%   Same as optimal_fraction with f_max = 1. The optional z = norm_ppf(alpha) is
%   passed in by best_mix_under_var, which runs this for thousands of funds.
if nargin < 8, alpha = 0.95; end
if nargin < 9, z = norm_ppf(alpha); end
[mu_d, s_d] = fund(d, mu, sigma, rho, r);
history = zeros(0, 4);
if mu_d <= r || s_d == 0          % never beats the bank
    f = 0; g = r;
    return
end
h = @(f) 1 - exp(growth_rate(f, mu_d, s_d, r) * T - z * f * s_d * sqrt(T)) - L;
if h(0) > 0
    error('risk:limit', 'risk limit L is below the VaR of the risk-free portfolio');
end
if h(1) <= 0
    f_var = 1;
else
    [f_var, history] = bisection(h, 0, 1);
end
f = min([kelly_fraction(mu_d, s_d, r), f_var, 1]);
g = growth_rate(f, mu_d, s_d, r);
end
