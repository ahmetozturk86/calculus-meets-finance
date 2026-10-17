function [f_es, history] = max_fraction_under_es(L, mu, sigma, r, T, alpha, f_max)
%MAX_FRACTION_UNDER_ES Largest f in [0, f_max] with ES(f) <= L: the same IVT + bisection as for VaR.
%   ES is continuous in f, and ln(1 - ES(f)) = (r + f (mu - r)) T + ln Phi(-z - f sigma sqrt(T)) - ln(1 - alpha)
%   is concave (Phi is log-concave), so the root is unique, as for VaR.
if nargin < 6, alpha = 0.95; end
if nargin < 7, f_max = 2; end
h = @(f) expected_shortfall(f, mu, sigma, r, T, alpha) - L;
if h(0) > 0
    error('risk:limit', 'risk limit L is below the ES of the risk-free portfolio');
end
if h(f_max) <= 0
    f_es = f_max;
    history = zeros(0, 4);
    return
end
[f_es, history] = bisection(h, 0, f_max);
end
