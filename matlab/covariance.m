function C = covariance(sigma, rho)
%COVARIANCE C_ij = rho_ij sigma_i sigma_j.
sigma = sigma(:)';
C = (sigma' * sigma) .* rho;
end
