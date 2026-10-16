function q = difference_quotient(fun, x, h)
%DIFFERENCE_QUOTIENT (fun(x + h) - fun(x)) / h. Its limit as h -> 0 is fun'(x).
q = (fun(x + h) - fun(x)) / h;
end
