"""Module 1 - Continuous compounding as a limit.

Compounding n times per year turns a principal P into
    A_n = P * (1 + r/n)^(n t)
and because lim_{n -> inf} (1 + r/n)^n = e^r, this converges to
    A = P * e^(r t)          (continuous compounding).

Rate of convergence: n log(1 + r/n) = r - r^2 / (2n) + O(1/n^2), so
    A - A_n  ~  A * r^2 t / (2n),
i.e. the extra money from compounding more often shrinks like 1/n.
"""
import numpy as np


def compound(P, r, t, n):
    """Value of P after t years at annual rate r, compounded n times per year."""
    n = np.asarray(n, dtype=float)
    # Same as (1 + r/n)^(n t); log1p keeps it accurate when r/n is tiny.
    return P * np.exp(n * t * np.log1p(r / n))


def continuous(P, r, t):
    """Value of P after t years with continuous compounding (the n -> inf limit)."""
    return P * np.exp(r * t)


def compounding_error(P, r, t, n):
    """Exact gap A - A_n between continuous and n-times-a-year compounding."""
    return continuous(P, r, t) - compound(P, r, t, n)


def error_estimate(P, r, t, n):
    """First-order estimate of the gap:  A * r^2 t / (2n)."""
    n = np.asarray(n, dtype=float)
    return continuous(P, r, t) * r**2 * t / (2.0 * n)
