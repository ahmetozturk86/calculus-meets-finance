function [x, pmf] = walk_log_return_distribution(mu, sigma, T, n_steps)
%WALK_LOG_RETURN_DISTRIBUTION Exact distribution of log(S_T/S_0) for the walk.
%   x are the possible values, pmf their probabilities (column vectors).
n = n_steps;
dt = T / n;
p = up_probability(mu, sigma, dt);
k = (0:n)';                                % number of up moves
% log C(n, k) built from C(n, k) = C(n, k-1) * (n - k + 1) / k
log_comb = [0; cumsum(log((n - k(2:end) + 1) ./ k(2:end)))];
pmf = exp(log_comb + k * log(p) + (n - k) * log1p(-p));
x = sigma * sqrt(dt) * (2 * k - n);
end
