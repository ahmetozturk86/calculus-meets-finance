function es = expected_shortfall(f, mu, sigma, r, T, alpha)
%EXPECTED_SHORTFALL The average loss in the worst 1 - alpha of cases (never below the VaR).
%   log(W_T/W_0) ~ N(m, s^2) with m = g(f) T and s = |f| sigma sqrt(T). Completing the square (as
%   for E[e^X]) gives E[W_T/W_0 | worst 1 - alpha] = e^(m + s^2/2) Phi(-z_alpha - s) / (1 - alpha).
if nargin < 6, alpha = 0.95; end
z = norm_ppf(alpha);
m = growth_rate(f, mu, sigma, r) * T;
s = abs(f) * sigma * sqrt(T);
es = 1 - exp(m + 0.5 * s.^2) .* norm_cdf(-z - s) / (1 - alpha);
end
