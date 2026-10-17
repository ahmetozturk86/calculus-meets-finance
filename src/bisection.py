"""Bisection root finding - the Intermediate Value Theorem turned into an algorithm.

Intermediate Value Theorem (IVT): if h is continuous on [a, b] and h(a), h(b)
have opposite signs, then there is a c in (a, b) with h(c) = 0.

Bisection halves [a, b] and keeps the half where the sign still changes, so by
the IVT that half still contains a root. After k steps the root is trapped in
an interval of length (b - a) / 2^k, and lim_{k -> inf} (b - a) / 2^k = 0, so
the midpoints converge to the root.
"""
import numpy as np


def bisection(fun, a, b, tol=1e-12, max_iter=200):
    """Find c in [a, b] with fun(c) = 0 for a continuous fun.

    Returns (root, history). Each history row is [a, b, mid, fun(mid)], so the
    shrinking interval can be printed or plotted. Signs are compared directly: the
    product of two tiny values could round to 0 and hide a sign change.
    """
    fa = float(fun(a))
    fb = float(fun(b))
    if fa == 0.0:
        return a, np.empty((0, 4))
    if fb == 0.0:
        return b, np.empty((0, 4))
    if (fa < 0) == (fb < 0):
        raise ValueError("IVT condition fails: fun(a) and fun(b) must have opposite signs")

    history = []
    mid = 0.5 * (a + b)
    for _ in range(max_iter):
        mid = 0.5 * (a + b)
        fm = float(fun(mid))
        history.append([a, b, mid, fm])
        if fm == 0.0 or 0.5 * (b - a) < tol:
            break
        if (fa < 0) != (fm < 0):
            b = mid
        else:
            a, fa = mid, fm
    return mid, np.array(history)
