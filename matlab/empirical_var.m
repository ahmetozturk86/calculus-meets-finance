function v = empirical_var(final_wealth, W0, alpha)
%EMPIRICAL_VAR VaR estimated from simulated final wealth.
if nargin < 2, W0 = 1; end
if nargin < 3, alpha = 0.95; end
v = 1 - quantile_linear(final_wealth / W0, 1 - alpha);
end
