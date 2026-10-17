function [f, history] = kelly_fraction_jump(mu, sigma, r, lam, J, f_hi)
%KELLY_FRACTION_JUMP Kelly fraction with crashes: bisection on the slope of GROWTH_RATE_JUMP (IVT).
%   The slope (mu - r) - f sigma^2 - lam J + lam J / (1 + f J) is continuous, equals mu - r > 0 at
%   f = 0, and goes to -infinity where one crash would wipe out the portfolio (1 + f J -> 0).
if nargin < 6, f_hi = 10; end
if J < 0, f_hi = min(f_hi, 0.999999 / -J); end
slope = @(f) (mu - r) - f * sigma^2 - lam * J + lam * J / (1 + f * J);
[f, history] = bisection(slope, 0, f_hi, 1e-10);
end
