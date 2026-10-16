function price = binomial_price(S, K, r, sigma, T, N, option)
%BINOMIAL_PRICE European option price on an N-step Cox-Ross-Rubinstein tree.
%   As N -> inf the price converges to BLACK_SCHOLES.
if nargin < 7, option = 'call'; end
dt = T / N;
u = exp(sigma * sqrt(dt));
d = 1 / u;
p = (exp(r * dt) - d) / (u - d);           % risk-neutral up probability
disc = exp(-r * dt);

j = (0:N)';                                % number of up moves
ST = S * exp(sigma * sqrt(dt) * (2 * j - N));
V = payoff(ST, K, option);
for k = 1:N                                % step back one level at a time
    V = disc * (p * V(2:end) + (1 - p) * V(1:end-1));
end
price = V(1);
end
