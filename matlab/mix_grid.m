function grid = mix_grid(n_assets, steps)
%MIX_GRID Every mix of n_assets with shares in {0, 1/steps, ..., 1} that add up to 1.
%   One mix per row, in the same order as src/portfolio.py (first share slowest).
if nargin < 2, steps = 100; end
rows = zeros(1, 0);
left = steps;
for a = 1:n_assets - 1                     % the next asset gets 0, 1, ..., whatever is left
    parts = arrayfun(@(k) (0:k)', left, 'UniformOutput', false);
    new = vertcat(parts{:});
    rows = [repelem(rows, left + 1, 1), new];
    left = repelem(left, left + 1, 1) - new;
end
grid = [rows, left] / steps;
end
