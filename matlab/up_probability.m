function p = up_probability(mu, sigma, dt)
%UP_PROBABILITY P(up) that makes the walk's log-drift equal (mu - sigma^2/2) per year.
p = 0.5 + 0.5 * (mu - 0.5 * sigma^2) / sigma * sqrt(dt);
end
