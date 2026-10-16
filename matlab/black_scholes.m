function price = black_scholes(S, K, r, sigma, T, option)
%BLACK_SCHOLES Closed-form price of a European option (the continuous limit).
if nargin < 6, option = 'call'; end
d1 = (log(S / K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T));
d2 = d1 - sigma * sqrt(T);
if strcmp(option, 'call')
    price = S * norm_cdf(d1) - K * exp(-r * T) * norm_cdf(d2);
elseif strcmp(option, 'put')
    price = K * exp(-r * T) * norm_cdf(-d2) - S * norm_cdf(-d1);
else
    error('black_scholes:option', 'option must be call or put');
end
end
