function f_opt = optimal_fraction_jump(mu, sigma, r, T, L, lam, J, alpha, f_max)
%OPTIMAL_FRACTION_JUMP Best growth when crashes can happen: min(f*_crash, f_crash, f_max), at least 0.
if nargin < 8, alpha = 0.95; end
if nargin < 9, f_max = 2; end
if mu <= r
    f_opt = 0;
    return
end
f_star = kelly_fraction_jump(mu, sigma, r, lam, J);
f_var = max_fraction_under_var_jump(L, mu, sigma, r, T, lam, J, alpha, f_max);
f_opt = max(0, min([f_star, f_var, f_max]));
end
