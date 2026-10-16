"""Module 2 - The binomial option model converges to Black-Scholes.

Cox-Ross-Rubinstein (CRR): split [0, T] into N steps of dt = T/N. Each step
the stock goes up by u = e^(sigma sqrt(dt)) or down by d = 1/u, with
risk-neutral probability p = (e^(r dt) - d) / (u - d). The option is priced by
working backwards through the tree.

As N -> inf (dt -> 0) the binomial price converges to the Black-Scholes price:
    lim_{N -> inf} C_N = C_BS.
"""
import math

import numpy as np

from src.normal import norm_cdf


def payoff(ST, K, option="call"):
    """Option value at expiry for stock prices ST."""
    if option == "call":
        return np.maximum(ST - K, 0.0)
    if option == "put":
        return np.maximum(K - ST, 0.0)
    raise ValueError("option must be 'call' or 'put'")


def black_scholes(S, K, r, sigma, T, option="call"):
    """Closed-form Black-Scholes price of a European option (the continuous limit)."""
    d1 = (math.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)
    if option == "call":
        return float(S * norm_cdf(d1) - K * math.exp(-r * T) * norm_cdf(d2))
    if option == "put":
        return float(K * math.exp(-r * T) * norm_cdf(-d2) - S * norm_cdf(-d1))
    raise ValueError("option must be 'call' or 'put'")


def binomial_price(S, K, r, sigma, T, N, option="call"):
    """European option price on an N-step CRR binomial tree."""
    dt = T / N
    u = math.exp(sigma * math.sqrt(dt))
    d = 1.0 / u
    p = (math.exp(r * dt) - d) / (u - d)
    disc = math.exp(-r * dt)

    j = np.arange(N + 1)                      # number of up moves
    ST = S * np.exp(sigma * math.sqrt(dt) * (2 * j - N))
    V = payoff(ST, K, option)
    for _ in range(N):                        # step back one level at a time
        V = disc * (p * V[1:] + (1.0 - p) * V[:-1])
    return float(V[0])
