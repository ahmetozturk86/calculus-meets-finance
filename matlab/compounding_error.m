function e = compounding_error(P, r, t, n)
%COMPOUNDING_ERROR Exact gap A - A_n between continuous and n-times-a-year compounding.
e = continuous(P, r, t) - compound(P, r, t, n);
end
