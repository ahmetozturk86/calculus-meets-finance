function grid = mix_grid(n_assets, steps)
%MIX_GRID Every mix of n_assets with shares in {0, 1/steps, ..., 1} that add up to 1.
%   One mix per row, in the same order as src/portfolio.py (first share slowest).
if nargin < 2, steps = 100; end
grid = counts(n_assets, steps) / steps;
end

function c = counts(n, total)
if n == 1
    c = total;
    return
end
c = zeros(0, n);
for i = 0:total
    rest = counts(n - 1, total - i);
    c = [c; repmat(i, size(rest, 1), 1), rest]; %#ok<AGROW>
end
end
