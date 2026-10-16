function F = gbm_log_return_cdf(x, mu, sigma, T)
%GBM_LOG_RETURN_CDF CDF of the GBM limit  N((mu - sigma^2/2) T, sigma^2 T).
m = (mu - 0.5 * sigma^2) * T;
s = sigma * sqrt(T);
F = norm_cdf((x - m) / s);
end
