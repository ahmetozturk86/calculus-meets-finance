function W = simulate_wealth_jump(f, mu, sigma, r, T, lam, J, n_steps, n_paths, W0)
%SIMULATE_WEALTH_JUMP Like simulate_wealth, but each step a crash happens with probability lam * dt.
%   As dt -> 0 the number of crashes in a year, Binomial(n_steps, lam dt),
%   tends to Poisson(lam T). Uses MATLAB's global random stream.
if nargin < 10, W0 = 1; end
dt = T / n_steps;
rf_return = exp(r * dt) - 1;
W = zeros(n_paths, n_steps + 1);
W(:, 1) = W0;
for k = 1:n_steps
    Z = randn(n_paths, 1);
    crash = rand(n_paths, 1) < lam * dt;
    stock_return = exp((mu - lam * J - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z) .* (1 + J * crash) - 1;
    W(:, k + 1) = W(:, k) .* (1 + f * stock_return + (1 - f) * rf_return);
end
end
