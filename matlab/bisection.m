function [root, history] = bisection(fun, a, b, tol, max_iter)
%BISECTION Find c in [a, b] with fun(c) = 0 for a continuous fun.
%   Intermediate Value Theorem: if fun is continuous on [a, b] and fun(a),
%   fun(b) have opposite signs, a root exists in (a, b). Halving the interval
%   and keeping the half with the sign change traps it in an interval of
%   length (b - a) / 2^k -> 0.
%   Each row of HISTORY is [a, b, mid, fun(mid)]. Signs are compared directly:
%   the product of two tiny values could round to 0 and hide a sign change.
if nargin < 4, tol = 1e-12; end
if nargin < 5, max_iter = 200; end
fa = fun(a);
fb = fun(b);
history = zeros(0, 4);
if fa == 0, root = a; return; end
if fb == 0, root = b; return; end
if (fa < 0) == (fb < 0)
    error('bisection:noSignChange', ...
        'IVT condition fails: fun(a) and fun(b) must have opposite signs');
end
mid = 0.5 * (a + b);
for k = 1:max_iter
    mid = 0.5 * (a + b);
    fm = fun(mid);
    history(end + 1, :) = [a, b, mid, fm]; %#ok<AGROW>
    if fm == 0 || 0.5 * (b - a) < tol
        break
    end
    if (fa < 0) ~= (fm < 0)
        b = mid;
    else
        a = mid;
        fa = fm;
    end
end
root = mid;
end
