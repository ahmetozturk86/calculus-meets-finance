function p = norm_cdf(x)
%NORM_CDF Standard normal CDF  Phi(x) = (1 + erf(x / sqrt(2))) / 2.
%   Written with erf so no Statistics Toolbox is needed.
p = 0.5 * (1 + erf(x ./ sqrt(2)));
end
