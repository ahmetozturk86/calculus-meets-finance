function p = up_probability(mu, sigma, dt)
%UP_PROBABILITY P(up) that makes the walk's log-drift equal (mu - sigma^2/2) per year.
%   It is a probability only while dt < sigma^2 / (mu - sigma^2/2)^2 (6.25 years for mu = 10%, sigma = 20%).
p = 0.5 + 0.5 * (mu - 0.5 * sigma^2) / sigma * sqrt(dt);
if ~(p > 0 && p < 1)
    error('gbm:step', 'step too long: need dt < sigma^2 / (mu - sigma^2/2)^2 so that 0 < P(up) < 1');
end
end
