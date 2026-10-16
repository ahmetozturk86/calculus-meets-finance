function d = ks_distance_to_gbm(mu, sigma, T, n_steps)
%KS_DISTANCE_TO_GBM Exact Kolmogorov-Smirnov distance sup_x |F_walk(x) - F_GBM(x)|.
%   F_walk is a step function, so the sup is reached just before or at a jump.
[x, pmf] = walk_log_return_distribution(mu, sigma, T, n_steps);
cdf_at = cumsum(pmf);                      % F_walk(x_k)
cdf_before = cdf_at - pmf;                 % F_walk(x_k^-)
F = gbm_log_return_cdf(x, mu, sigma, T);
d = max(max(abs(cdf_at - F), abs(cdf_before - F)));
end
