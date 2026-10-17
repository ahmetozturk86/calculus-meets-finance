%% Tests for the MATLAB port - same checks as tests/test_modules.py
% Run from this folder:  >> run_tests
MU = 0.10; SIGMA = 0.20; R = 0.03; T = 1; ALPHA = 0.95; L = 0.15;
close_to = @(a, b, tol) abs(a - b) <= tol;

% bisection / normal
[root, history] = bisection(@(x) x^2 - 2, 0, 2);
assert(close_to(root, sqrt(2), 1e-11), 'bisection sqrt(2)');
assert(all(diff(history(:, 2) - history(:, 1)) < 0), 'bisection interval shrinks');
try
    bisection(@(x) x^2 + 1, -1, 1);
    error('tests:noError', 'bisection should fail without a sign change');
catch err
    assert(strcmp(err.identifier, 'bisection:noSignChange'), 'bisection sign check');
end
assert(close_to(norm_cdf(0), 0.5, 1e-15), 'norm_cdf(0)');
assert(close_to(norm_ppf(0.95), 1.6448536269514722, 1e-10), 'norm_ppf(0.95)');

% module 1
assert(close_to(continuous(1, 0.1, 1), 1.105170918, 1e-9), 'e^0.1');
values = compound(1, 0.1, 1, [1 12 365 1e6]);
assert(all(diff(values) > 0), 'compounding increases with n');
assert(compounding_error(1, 0.1, 1, 1e6) < 1e-7, 'compounding converges');
ratio = compounding_error(1, 0.1, 1, 1e5) / error_estimate(1, 0.1, 1, 1e5);
assert(close_to(ratio, 1, 1e-4), 'first-order error estimate');

% module 2
call = black_scholes(100, 100, 0.05, 0.2, 1, 'call');
put = black_scholes(100, 100, 0.05, 0.2, 1, 'put');
assert(close_to(call, 10.4506, 1e-4), 'Black-Scholes reference');
assert(close_to(call - put, 100 - 100 * exp(-0.05), 1e-10), 'put-call parity');
assert(abs(binomial_price(100, 100, 0.05, 0.2, 1, 1000, 'call') - call) < 0.01, 'binomial call');
assert(abs(binomial_price(100, 100, 0.05, 0.2, 1, 1000, 'put') - put) < 0.01, 'binomial put');

% module 3
[x, pmf] = walk_log_return_distribution(MU, SIGMA, T, 500);
assert(close_to(sum(pmf), 1, 1e-12), 'pmf sums to 1');
assert(close_to(x' * pmf, (MU - 0.5 * SIGMA^2) * T, 1e-12), 'walk mean matches GBM');
ks = arrayfun(@(n) ks_distance_to_gbm(MU, SIGMA, T, n), [1 100 10000]);
assert(ks(1) > ks(2) && ks(2) > ks(3) && ks(3) < 0.01, 'KS distance goes to 0');

% module 4
assert(close_to(kelly_fraction(MU, SIGMA, R), 1.75, 1e-12), 'Kelly fraction');
assert(close_to(numerical_kelly_fraction(MU, SIGMA, R), 1.75, 1e-6), 'numerical Kelly');
assert(close_to(value_at_risk(0, MU, SIGMA, R, T, ALPHA), 1 - exp(R * T), 1e-12), 'VaR of risk-free');
f_var = max_fraction_under_var(L, MU, SIGMA, R, T, ALPHA);
assert(abs(value_at_risk(f_var, MU, SIGMA, R, T, ALPHA) - L) < 1e-8, 'VaR root hits limit');
assert(close_to(optimal_fraction(MU, SIGMA, R, T, L, ALPHA), min(1.75, f_var), 1e-12), 'optimal fraction');
assert(max_fraction_under_var(0.99, MU, SIGMA, R, T, ALPHA, 2) == 2, 'non-binding limit');
rng(0);
W = simulate_wealth(0.7, MU, SIGMA, R, T, 52, 50000);
assert(close_to(empirical_var(W(:, end), 1, ALPHA), value_at_risk(0.7, MU, SIGMA, R, T, ALPHA), 0.01), 'Monte Carlo VaR');

% module 4b: fat tails (rare crashes)
LAM = 0.10; J = -0.25;
for f = [0 0.3 0.705 1.5]
    for a = [0.9 0.95 0.99]
        assert(close_to(value_at_risk_jump(f, MU, SIGMA, R, T, 0, J, a), value_at_risk(f, MU, SIGMA, R, T, a), 1e-10), ...
            'no crashes gives the normal VaR');
    end
end
errs = arrayfun(@(n) max(arrayfun(@(k) abs(binomial_pmf(k, n, LAM / n) - poisson_pmf(k, LAM)), 0:30)), [10 1000 100000]);
assert(errs(1) > errs(2) && errs(2) > errs(3) && errs(3) < 1e-6, 'binomial crash count converges to Poisson');
f_crash = max_fraction_under_var_jump(L, MU, SIGMA, R, T, LAM, J, ALPHA);
assert(abs(value_at_risk_jump(f_crash, MU, SIGMA, R, T, LAM, J, ALPHA) - L) < 1e-9, 'crash root hits limit');
assert(f_crash < f_var, 'crashes lower the stock fraction');
assert(value_at_risk_jump(f_var, MU, SIGMA, R, T, LAM, J, ALPHA) > L, 'crash VaR above the normal VaR');
rng(0);
W = simulate_wealth_jump(0.7, MU, SIGMA, R, T, LAM, J, 252, 50000);
assert(close_to(empirical_var(W(:, end), 1, ALPHA), value_at_risk_jump(0.7, MU, SIGMA, R, T, LAM, J, ALPHA), 0.01), ...
    'Monte Carlo VaR with crashes');

% module 5: more assets
MU3 = [MU, log(1.05), log(1.045)]; SIG3 = [SIGMA, 0.15, 0.06]; RHO3 = [1 0 -0.2; 0 1 0.3; -0.2 0.3 1];
for f = [0 0.3 0.705 1]
    assert(close_to(portfolio_var([f 0 0], MU3, SIG3, RHO3, R, T, ALPHA), value_at_risk(f, MU, SIGMA, R, T, ALPHA), 1e-12), ...
        'one-asset mix is Module 4');
end
assert(close_to(fund_fraction([1 0 0], L, MU3, SIG3, RHO3, R, T, ALPHA), max_fraction_under_var(L, MU, SIGMA, R, T, ALPHA, 1), 1e-12), ...
    'stocks-only fund is Module 4');
grid5 = mix_grid(3, 100);
assert(size(grid5, 1) == 5151 && all(abs(sum(grid5, 2) - 1) < 1e-12) && min(grid5(:)) >= 0, 'mix grid');
[w5, f5, d5] = best_mix_under_var(L, MU3, SIG3, RHO3, R, T, ALPHA);
assert(all(abs(w5 - f5 * d5) < 1e-12) && min(w5) >= 0 && sum(w5) <= 1 + 1e-12, 'best mix weights');
assert(portfolio_var(w5, MU3, SIG3, RHO3, R, T, ALPHA) <= L + 1e-9, 'best mix respects the limit');
g5 = portfolio_growth(w5, MU3, SIG3, RHO3, R);
E = eye(3);
for i = 1:3
    [~, gi] = fund_fraction(E(i, :), L, MU3, SIG3, RHO3, R, T, ALPHA);
    assert(g5 >= gi, 'best mix beats one asset + bank');
end
[w1, f1] = best_mix_under_var(L, [MU R R], SIG3, eye(3), R, T, ALPHA);
assert(all(abs(w1 - [f1 0 0]) < 1e-12) && close_to(f1, 0.705014, 1e-6), 'stocks only when the others do not pay');
rng(0);
W = simulate_wealth_mix([0.5 0.1 0.3], MU3, SIG3, RHO3, R, T, 52, 50000);
assert(close_to(empirical_var(W(:, end), 1, ALPHA), portfolio_var([0.5 0.1 0.3], MU3, SIG3, RHO3, R, T, ALPHA), 0.01), ...
    'Monte Carlo VaR of a mix');

fprintf('All MATLAB tests passed.\n');
