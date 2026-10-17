function [lo, hi, count] = near_optimal_ranges(L, mu, sigma, rho, r, T, alpha, steps, tol)
%NEAR_OPTIMAL_RANGES Lowest and highest weight of each asset over the grid mixes growing within tol of the best.
%   The top of the growth surface is flat, so many mixes are almost as good.
if nargin < 7, alpha = 0.95; end
if nargin < 8, steps = 100; end
if nargin < 9, tol = 1e-4; end
D = mix_grid(numel(mu), steps);
[f, g] = fund_fractions(D, L, mu, sigma, rho, r, T, alpha);
W = f .* D;
near = g >= max(g) - tol;
lo = min(W(near, :), [], 1);
hi = max(W(near, :), [], 1);
count = sum(near);
end
