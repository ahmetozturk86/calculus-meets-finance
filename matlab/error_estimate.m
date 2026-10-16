function e = error_estimate(P, r, t, n)
%ERROR_ESTIMATE First-order estimate of the gap:  A r^2 t / (2n).
%   From n log(1 + r/n) = r - r^2 / (2n) + O(1/n^2).
e = continuous(P, r, t) .* r.^2 .* t ./ (2 .* n);
end
