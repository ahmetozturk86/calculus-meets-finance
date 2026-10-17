function [f, g] = fund_fractions(D, L, mu, sigma, rho, r, T, alpha)
%FUND_FRACTIONS Module 4 for many funds at once (one per row of D): f = min(f*, f_VaR, 1) and the growth.
%   Same answer as FUND_FRACTION, but f_VaR comes from the quadratic formula (VAR_ROOT_FORMULA)
%   instead of bisection, so grids with millions of mixes are fast.
if nargin < 8, alpha = 0.95; end
mu = mu(:);
a = D * (mu - r);                                       % mu_d - r
s = sqrt(sum((D * covariance(sigma, rho)) .* D, 2));    % sigma_d
f = zeros(size(D, 1), 1);
ok = a > 0 & s > 0;                                     % the others never beat the bank
if any(ok)
    f(ok) = min(min(a(ok) ./ s(ok).^2, var_root_formula(L, r + a(ok), s(ok), r, T, alpha)), 1);
end
g = r + f .* a - 0.5 * f.^2 .* s.^2;
end
