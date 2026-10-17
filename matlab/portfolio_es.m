function v = portfolio_es(w, mu, sigma, rho, r, T, alpha)
%PORTFOLIO_ES Multi-asset Expected Shortfall. With one asset it equals EXPECTED_SHORTFALL (Module 4).
if nargin < 7, alpha = 0.95; end
[g, sd] = portfolio_growth(w, mu, sigma, rho, r);
s = sd * sqrt(T);
v = 1 - exp(g * T + 0.5 * s^2) * norm_cdf(-norm_ppf(alpha) - s) / (1 - alpha);
end
