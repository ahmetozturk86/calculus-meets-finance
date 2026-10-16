function paths = simulate_walk_paths(S0, mu, sigma, T, n_steps, n_paths)
%SIMULATE_WALK_PATHS Price paths of the coin-flip random walk, size n_paths x (n_steps+1).
%   Uses MATLAB's global random stream (set it with rng(seed)).
dt = T / n_steps;
p = up_probability(mu, sigma, dt);
steps = 2 * (rand(n_paths, n_steps) < p) - 1;
log_paths = cumsum(sigma * sqrt(dt) * steps, 2);
log_paths = [zeros(n_paths, 1), log_paths];
paths = S0 * exp(log_paths);
end
