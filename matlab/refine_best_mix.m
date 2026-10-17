function [w, f, d, g] = refine_best_mix(L, mu, sigma, rho, r, T, alpha, steps, levels, points)
%REFINE_BEST_MIX Zoom in on the best mix: let the grid step go to 0 around the best fund.
%   Start from the best fund on the 1/steps grid. Each level tries funds on a small grid around the
%   current best (+-10 coarse steps at first, POINTS intervals per direction) and then halves that
%   window. The current best is always tried again, so the growth never goes down.
if nargin < 7, alpha = 0.95; end
if nargin < 8, steps = 100; end
if nargin < 9, levels = 40; end
if nargin < 10, points = 100; end
[w, f, d, g] = best_mix_on_grid(L, mu, sigma, rho, r, T, alpha, steps);
if isempty(d), return; end
n = numel(mu);
offsets = linspace(-10 / steps, 10 / steps, points + 1);   % contains 0 (points is even)
for level = 1:levels
    grids = cell(1, n - 1);
    [grids{end:-1:1}] = ndgrid(offsets);                   % first move slowest, as in Python
    shift = cell2mat(cellfun(@(x) x(:), grids, 'UniformOutput', false));
    D = [d(1:end-1) + shift, d(end) - sum(shift, 2)];      % moves that keep sum(d) = 1
    D = max(D(all(D > -1e-12, 2), :), 0);
    [fs, gs] = fund_fractions(D, L, mu, sigma, rho, r, T, alpha);
    [gk, k] = max(gs);
    if gk > g
        f = fs(k); d = D(k, :); g = gk;
    end
    offsets = offsets / 2;
end
w = f * d;
end
