function v = portfolio_var(w, mu, sigma, rho, r, T, alpha)
%PORTFOLIO_VAR Multi-asset Value-at-Risk. With one asset it equals value_at_risk (Module 4).
if nargin < 7, alpha = 0.95; end
[g, sd] = portfolio_growth(w, mu, sigma, rho, r);
v = 1 - exp(g * T - norm_ppf(alpha) * sd * sqrt(T));
end
