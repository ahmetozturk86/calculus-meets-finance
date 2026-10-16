function W = simulate_wealth(f, mu, sigma, r, T, n_steps, n_paths, W0)
%SIMULATE_WEALTH Monte Carlo wealth paths, rebalanced to fraction f every dt = T / n_steps.
%   Returns an n_paths x (n_steps+1) matrix. As dt -> 0 this converges to the
%   continuous-rebalancing formula. Uses MATLAB's global random stream.
if nargin < 8, W0 = 1; end
dt = T / n_steps;
rf_return = exp(r * dt) - 1;
W = zeros(n_paths, n_steps + 1);
W(:, 1) = W0;
for k = 1:n_steps
    Z = randn(n_paths, 1);
    stock_return = exp((mu - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z) - 1;
    W(:, k + 1) = W(:, k) .* (1 + f * stock_return + (1 - f) * rf_return);
end
end
