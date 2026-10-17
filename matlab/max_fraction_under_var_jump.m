function [f_var, history] = max_fraction_under_var_jump(L, mu, sigma, r, T, lam, J, alpha, f_max)
%MAX_FRACTION_UNDER_VAR_JUMP Largest f in [0, f_max] with VaR_jump(f) <= L.
%   h(f) = VaR_jump(f) - L is continuous with h(0) < 0 < h(f_max), so by the
%   IVT a root exists; bisection finds it. If even f_max respects the limit
%   the constraint is not binding and f_max is returned.
if nargin < 8, alpha = 0.95; end
if nargin < 9, f_max = 2; end
h = @(f) value_at_risk_jump(f, mu, sigma, r, T, lam, J, alpha) - L;
if h(0) > 0
    error('risk:limit', 'risk limit L is below the VaR of the risk-free portfolio');
end
if h(f_max) <= 0
    f_var = f_max;
    history = zeros(0, 4);
    return
end
[f_var, history] = bisection(h, 0, f_max);
end
