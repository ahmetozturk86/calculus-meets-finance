function [g, sd] = portfolio_growth(w, mu, sigma, rho, r)
%PORTFOLIO_GROWTH Long-run growth g(w) = r + w.(mu - r) - w'Cw / 2 and spread sd(w) = sqrt(w'Cw).
w = w(:)'; mu = mu(:)';
v = w * covariance(sigma, rho) * w';
g = r + sum(w .* (mu - r)) - 0.5 * v;
sd = sqrt(v);
end
