function g = growth_rate(f, mu, sigma, r)
%GROWTH_RATE Long-run growth g(f) = r + f (mu - r) - f^2 sigma^2 / 2.
g = r + f .* (mu - r) - 0.5 * f.^2 * sigma^2;
end
