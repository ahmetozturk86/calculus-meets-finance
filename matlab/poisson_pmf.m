function p = poisson_pmf(k, m)
%POISSON_PMF P(N = k) for N ~ Poisson(m).
if m == 0
    p = double(k == 0);
    return
end
p = exp(-m + k * log(m) - gammaln(k + 1));
end
