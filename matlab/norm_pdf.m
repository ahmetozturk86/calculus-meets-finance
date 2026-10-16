function y = norm_pdf(x)
%NORM_PDF Standard normal density phi(x).
y = exp(-0.5 * x.^2) / sqrt(2 * pi);
end
