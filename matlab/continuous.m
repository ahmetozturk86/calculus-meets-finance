function A = continuous(P, r, t)
%CONTINUOUS Continuous compounding P e^(r t), the n -> inf limit of COMPOUND.
A = P .* exp(r .* t);
end
