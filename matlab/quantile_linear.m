function v = quantile_linear(x, q)
%QUANTILE_LINEAR Sample quantiles with linear interpolation (same as numpy.quantile).
%   Avoids the Statistics Toolbox. Returns a column vector, one value per q.
x = sort(x(:));
q = q(:);
n = numel(x);
pos = (n - 1) * q + 1;
lo = floor(pos);
hi = min(lo + 1, n);
v = x(lo) + (pos - lo) .* (x(hi) - x(lo));
end
