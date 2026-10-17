function W = simulate_wealth_mix(w, mu, sigma, rho, r, T, n_steps, n_paths, W0)
%SIMULATE_WEALTH_MIX Monte Carlo wealth paths, rebalanced to the weights w every dt = T / n_steps.
%   The risky assets get correlated GBM returns (Cholesky factor of C) and the rest
%   earns e^(r dt). Returns an n_paths x (n_steps+1) matrix. Uses MATLAB's global random stream.
if nargin < 9, W0 = 1; end
w = w(:); mu = mu(:)'; sigma = sigma(:)';
chol_C = chol(covariance(sigma, rho), 'lower');
dt = T / n_steps;
W = zeros(n_paths, n_steps + 1);
W(:, 1) = W0;
for k = 1:n_steps
    Z = randn(n_paths, numel(mu)) * chol_C';
    gross = exp((mu - 0.5 * sigma.^2) * dt + Z * sqrt(dt));
    W(:, k + 1) = W(:, k) .* (gross * w + (1 - sum(w)) * exp(r * dt));
end
end
