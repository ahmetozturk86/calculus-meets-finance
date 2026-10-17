function g = growth_rate_jump(f, mu, sigma, r, lam, J)
%GROWTH_RATE_JUMP Long-run growth with crashes, E[log(W_T/W_0)]/T = g(f) - lam J f + lam ln(1 + f J).
%   The average return is the same as without crashes, but ln(1 + x) < x for x ~= 0, so crashes
%   lower the growth: a big loss hurts more than the extra drift of the good years helps.
g = growth_rate(f, mu, sigma, r) - lam * J * f + lam * log1p(f * J);
end
