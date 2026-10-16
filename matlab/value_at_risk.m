function v = value_at_risk(f, mu, sigma, r, T, alpha)
%VALUE_AT_RISK Loss (fraction of starting wealth) exceeded with probability 1 - alpha.
%   The (1 - alpha) quantile of W_T / W_0 is exp(g(f) T - z_alpha |f| sigma sqrt(T)).
%   A negative VaR means even the bad case is a gain.
if nargin < 6, alpha = 0.95; end
z = norm_ppf(alpha);
v = 1 - exp(growth_rate(f, mu, sigma, r) * T - z * abs(f) * sigma * sqrt(T));
end
