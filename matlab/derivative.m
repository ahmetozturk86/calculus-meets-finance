function d = derivative(fun, x, h)
%DERIVATIVE Numerical derivative (central difference quotient).
if nargin < 3, h = 1e-6; end
d = (fun(x + h) - fun(x - h)) / (2 * h);
end
