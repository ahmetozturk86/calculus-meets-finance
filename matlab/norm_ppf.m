function z = norm_ppf(q)
%NORM_PPF Inverse standard normal CDF: the z with Phi(z) = q.
%   Phi is continuous and increasing, so Phi(z) - q changes sign on [-10, 10]
%   and the IVT guarantees a solution, which bisection finds.
if q <= 0 || q >= 1
    error('norm_ppf:range', 'q must be in (0, 1)');
end
z = bisection(@(x) norm_cdf(x) - q, -10, 10, 1e-13);
end
