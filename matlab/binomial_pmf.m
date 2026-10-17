function p = binomial_pmf(k, n, p_success)
%BINOMIAL_PMF P(B = k) for B ~ Binomial(n, p).
%   As n -> inf with p = lam / n it tends to poisson_pmf(k, lam).
%   Written with gammaln (not nchoosek) so large n stays exact enough.
if p_success == 0
    p = double(k == 0);
    return
end
if k > n
    p = 0;
    return
end
p = exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) ...
    + k * log(p_success) + (n - k) * log1p(-p_success));
end
