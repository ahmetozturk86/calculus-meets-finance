function f_opt = optimal_fraction(mu, sigma, r, T, L, alpha, f_max)
%OPTIMAL_FRACTION Best growth subject to the VaR limit: min(f*, f_VaR, f_max), at least 0.
if nargin < 6, alpha = 0.95; end
if nargin < 7, f_max = 2; end
f_star = kelly_fraction(mu, sigma, r);
f_var = max_fraction_under_var(L, mu, sigma, r, T, alpha, f_max);
f_opt = max(0, min([f_star, f_var, f_max]));
end
