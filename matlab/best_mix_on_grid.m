function [w, f, d, g] = best_mix_on_grid(L, mu, sigma, rho, r, T, alpha, steps)
%BEST_MIX_ON_GRID The fastest-growing fund on the grid with step 1/steps (FUND_FRACTIONS for all of them).
%   With steps = 100 it is the mix of BEST_MIX_UNDER_VAR. If one grid's steps divide the other's, the
%   finer grid contains the coarser one, so its best growth is the same or higher.
if nargin < 7, alpha = 0.95; end
if nargin < 8, steps = 100; end
D = mix_grid(numel(mu), steps);
[fs, gs] = fund_fractions(D, L, mu, sigma, rho, r, T, alpha);
[g, k] = max(gs);
if g <= r + 1e-15                                       % nothing beats the bank
    w = zeros(1, numel(mu)); f = 0; d = []; g = r;
    return
end
f = fs(k);
d = D(k, :);
w = f * d;
end
