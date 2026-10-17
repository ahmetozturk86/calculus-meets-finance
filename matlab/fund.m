function [mu_d, s_d] = fund(d, mu, sigma, rho, r)
%FUND The mix d of the risky assets seen as one fund: mu_d = r + d.(mu - r), s_d = sqrt(d'Cd).
d = d(:)'; mu = mu(:)';
mu_d = r + sum(d .* (mu - r));
s_d = sqrt(d * covariance(sigma, rho) * d');
end
