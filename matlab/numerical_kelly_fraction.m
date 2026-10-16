function [f, history] = numerical_kelly_fraction(mu, sigma, r, f_lo, f_hi)
%NUMERICAL_KELLY_FRACTION Maximize g without the formula: solve g'(f) = 0 by bisection.
if nargin < 4, f_lo = 0; end
if nargin < 5, f_hi = 10; end
g = @(f) growth_rate(f, mu, sigma, r);
[f, history] = bisection(@(f) derivative(g, f), f_lo, f_hi, 1e-10);
end
