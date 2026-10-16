function A = compound(P, r, t, n)
%COMPOUND Value of P after t years at annual rate r, compounded n times per year.
%   Same as P (1 + r/n)^(n t); log1p keeps it accurate when r/n is tiny.
A = P .* exp(n .* t .* log1p(r ./ n));
end
