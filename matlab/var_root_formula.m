function f = var_root_formula(L, mu, sigma, r, T, alpha)
%VAR_ROOT_FORMULA The root of VaR(f) = L from the quadratic formula, a check on the bisection.
%   VaR(f) = L  <=>  (sigma^2 T/2) f^2 + (z sigma sqrt(T) - (mu - r) T) f + (ln(1 - L) - r T) = 0.
%   The last coefficient is negative exactly when h(0) < 0, so the product of the two roots is
%   negative: one root is negative and one positive, for every horizon T. The positive one is
%   the only f >= 0 with VaR(f) = L. Works on arrays of mu and sigma (one fund per entry).
if nargin < 6, alpha = 0.95; end
if 1 - exp(r * T) - L > 0
    error('risk:limit', 'risk limit L is below the VaR of the risk-free portfolio');
end
a = 0.5 * sigma.^2 * T;
b = norm_ppf(alpha) * sigma * sqrt(T) - (mu - r) * T;
if L >= 1                                  % no limit: a loss of 100% never has 5% chance
    f = inf(size(a + b));
    return
end
c = log(1 - L) - r * T;
disc = sqrt(b.^2 - 4 * a * c);
f = (disc - b) ./ (2 * a);                 % two equal forms; each avoids cancellation
pos = b >= 0;
f(pos) = -2 * c ./ (b(pos) + disc(pos));
end
