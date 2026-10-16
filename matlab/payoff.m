function V = payoff(ST, K, option)
%PAYOFF Option value at expiry for stock prices ST ('call' or 'put').
if nargin < 3, option = 'call'; end
if strcmp(option, 'call')
    V = max(ST - K, 0);
elseif strcmp(option, 'put')
    V = max(K - ST, 0);
else
    error('payoff:option', 'option must be call or put');
end
end
