function [w, f, d] = best_mix_under_var(L, mu, sigma, rho, r, T, alpha, steps)
%BEST_MIX_UNDER_VAR Fastest-growing mix with VaR <= L, no borrowing, no short selling.
%   Every mix d of the risky assets is one fund; fund_fraction (Module 4: IVT +
%   bisection) gives each fund's best fraction f. Returns the weights w = f d in the
%   risky assets (the bank gets 1 - sum(w)), f and the winning fund d.
if nargin < 7, alpha = 0.95; end
if nargin < 8, steps = 100; end
z = norm_ppf(alpha);
grid = mix_grid(numel(mu), steps);
best_g = r; f = 0; d = [];
for k = 1:size(grid, 1)
    [fk, gk] = fund_fraction(grid(k, :), L, mu, sigma, rho, r, T, alpha, z);
    if gk > best_g + 1e-15
        best_g = gk; f = fk; d = grid(k, :);
    end
end
if isempty(d)
    w = zeros(1, numel(mu));
else
    w = f * d;
end
end
