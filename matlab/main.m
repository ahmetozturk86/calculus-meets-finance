%% Limits in Finance - MATLAB version of main.py
% Runs all modules, prints the results and saves the figures to
% matlab/figures/. Run from this folder:  >> main
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
fig_dir = fullfile(here, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
rng(42);

% Module 1: savings account
P = 10000; R_SAVE = 0.10; T_SAVE = 10;
% Module 2: European call option
S0 = 100; K = 100; R_OPT = 0.05; SIGMA_OPT = 0.20; T_OPT = 1;
% Modules 3-4: stock market and risk-free asset (yearly, continuously compounded: mu = 10% is 10.5% a year)
MU = 0.10; SIGMA = 0.20; R = 0.03; T = 1;
ALPHA = 0.95; LOSS_LIMIT = 0.15; F_MAX = 2;
% Module 4b: rare crashes (fat tails): the stock drops 25%, on average once every 10 years
CRASH_RATE = 0.10; CRASH_SIZE = -0.25;
% Module 5: more assets. Stocks as above; gold and bonds are illustrative long-run numbers (5% and 4.5% a year as
% simple rates, so log(1.05) and log(1.045) as continuous rates), not data. Same numbers as the "Normal" market on the website.
ASSET_NAMES = {'stocks', 'gold', 'bonds'};
MU_ASSETS = [MU, log(1.05), log(1.045)];
SIGMA_ASSETS = [SIGMA, 0.15, 0.06];
RHO_ASSETS = [1 0 -0.2; 0 1 0.3; -0.2 0.3 1];

RED = [0.84 0.15 0.16]; GREEN = [0.17 0.63 0.17];

%% Module 1 - Continuous compounding
header('Module 1 - Continuous compounding:  lim (1 + r/n)^(nt) = e^(rt)');
A_inf = continuous(P, R_SAVE, T_SAVE);
fprintf('  P = %.0f, r = %.0f%%, t = %.0f years  ->  P e^(rt) = %.4f\n', P, 100 * R_SAVE, T_SAVE, A_inf);
fprintf('  %16s %14s %12s %12s\n', 'periods/year n', 'A_n', 'A - A_n', 'estimate');
labels = {'yearly', 'quarterly', 'monthly', 'weekly', 'daily', 'hourly', 'every minute'};
periods = [1 4 12 52 365 8760 525600];
for i = 1:numel(periods)
    n = periods(i);
    fprintf('  %16s %14.4f %12.4f %12.4f\n', labels{i}, compound(P, R_SAVE, T_SAVE, n), ...
        compounding_error(P, R_SAVE, T_SAVE, n), error_estimate(P, R_SAVE, T_SAVE, n));
end

ns = unique(floor(logspace(0, 7, 60)));
fig = figure('Position', [100 100 1100 420]);
subplot(1, 2, 1);
semilogx(ns, compound(P, R_SAVE, T_SAVE, ns), 'o-', 'MarkerSize', 3); hold on;
yline(A_inf, '--', 'Color', RED);
xlabel('compounding periods per year n'); ylabel('final value');
title('Discrete compounding converges to Pe^{rt}');
legend('A_n = P(1+r/n)^{nt}', 'limit Pe^{rt}', 'Location', 'southeast');
subplot(1, 2, 2);
loglog(ns, compounding_error(P, R_SAVE, T_SAVE, ns), 'o', 'MarkerSize', 3); hold on;
loglog(ns, error_estimate(P, R_SAVE, T_SAVE, ns), '-');
xlabel('n'); ylabel('gap'); title('The gap shrinks like 1/n');
legend('actual gap A - A_n', 'estimate Pe^{rt} r^2 t / (2n)');
save_fig(fig, fig_dir, 'fig1_compounding.png');

%% Module 2 - Binomial tree -> Black-Scholes
header('Module 2 - Binomial tree -> Black-Scholes:  lim_{N->inf} C_N = C_BS');
bs = black_scholes(S0, K, R_OPT, SIGMA_OPT, T_OPT, 'call');
fprintf('  S = K = %.0f, r = %.0f%%, sigma = %.0f%%, T = %.0f  ->  C_BS = %.6f\n', ...
    S0, 100 * R_OPT, 100 * SIGMA_OPT, T_OPT, bs);
fprintf('  %6s %12s %14s\n', 'N', 'C_N', '|C_N - C_BS|');
for N = [1 2 5 10 50 100 500 1000]
    c = binomial_price(S0, K, R_OPT, SIGMA_OPT, T_OPT, N, 'call');
    fprintf('  %6d %12.6f %14.6f\n', N, c, abs(c - bs));
end

Ns = 1:200;
prices = arrayfun(@(N) binomial_price(S0, K, R_OPT, SIGMA_OPT, T_OPT, N, 'call'), Ns);
Ns_log = unique(floor(logspace(0, 3.3, 40)));
errors = arrayfun(@(N) abs(binomial_price(S0, K, R_OPT, SIGMA_OPT, T_OPT, N, 'call') - bs), Ns_log);

fig = figure('Position', [100 100 1100 420]);
subplot(1, 2, 1);
plot(Ns, prices, '-'); hold on;
yline(bs, '--', 'Color', RED);
xlabel('tree steps N'); ylabel('call price'); title('Binomial price converges to Black-Scholes');
legend('binomial price C_N', sprintf('Black-Scholes %.4f', bs));
subplot(1, 2, 2);
loglog(Ns_log, errors, 'o', 'MarkerSize', 4); hold on;
loglog(Ns_log, errors(1) ./ Ns_log, '--');
xlabel('N'); ylabel('error'); title('Error goes to 0 like 1/N');
legend('|C_N - C_{BS}|', 'reference \propto 1/N');
save_fig(fig, fig_dir, 'fig2_binomial_to_black_scholes.png');

%% Module 3 - Random walk -> Geometric Brownian Motion
header('Module 3 - Random walk -> Geometric Brownian Motion as dt -> 0');
fprintf('  mu = %.0f%%, sigma = %.0f%%, T = %.0f year\n', 100 * MU, 100 * SIGMA, T);
fprintf('  %8s %12s %12s\n', 'steps n', 'dt (years)', 'KS distance');
for n = [1 4 12 52 252 2520 25200]
    fprintf('  %8d %12.6f %12.6f\n', n, T / n, ks_distance_to_gbm(MU, SIGMA, T, n));
end

fig = figure('Position', [100 100 1100 800]);
steps_list = [12 252];
for i = 1:2
    n = steps_list(i);
    subplot(2, 2, i);
    paths = simulate_walk_paths(S0, MU, SIGMA, T, n, 30);
    plot(linspace(0, T, n + 1), paths');
    xlabel('time (years)'); ylabel('price'); title(sprintf('Random-walk paths, n = %d steps', n));
end

subplot(2, 2, 3); hold on;
cdf_labels = {};
for n = [2 12 252]
    [x, pmf] = walk_log_return_distribution(MU, SIGMA, T, n);
    % pad with F = 0 on the left and F = 1 on the right so the step CDF spans the axis
    stairs([-1; x; 1], [0; cumsum(pmf); 1]);
    cdf_labels{end + 1} = sprintf('walk, n = %d', n); %#ok<SAGROW>
end
x_grid = linspace(-0.7, 0.8, 400);
plot(x_grid, gbm_log_return_cdf(x_grid, MU, SIGMA, T), 'k--', 'LineWidth', 1.5);
xlim([-0.7 0.8]); box on;
xlabel('log-return log(S_T/S_0)'); ylabel('CDF'); title('Distribution converges to the GBM normal law');
legend([cdf_labels, {'GBM limit'}], 'Location', 'northwest');

subplot(2, 2, 4);
ns = unique(floor(logspace(0, 4.5, 30)));
ks = arrayfun(@(n) ks_distance_to_gbm(MU, SIGMA, T, n), ns);
loglog(ns, ks, 'o', 'MarkerSize', 4); hold on;
loglog(ns, ks(end) * sqrt(ns(end) ./ ns), '--');
xlabel('steps n = T/dt'); ylabel('sup_x |F_n(x) - F(x)|'); title('KS distance goes to 0');
legend('KS distance', 'reference \propto 1/\surd n');
save_fig(fig, fig_dir, 'fig3_random_walk_to_gbm.png');

%% Module 4 - Risk optimization
header('Module 4 - Risk optimization: best growth under a Value-at-Risk limit');
fprintf('  stock: mu = %.0f%%, sigma = %.0f%%;  risk-free r = %.0f%%;  horizon T = %.0f year\n', ...
    100 * MU, 100 * SIGMA, 100 * R, T);
fprintf('  rule: with %.0f%% confidence, lose at most %.0f%% in one year;  leverage <= %.0fx\n', ...
    100 * ALPHA, 100 * LOSS_LIMIT, F_MAX);

g = @(f) growth_rate(f, MU, SIGMA, R);
f_test = 1.0;
fprintf('\n  Derivative as a limit: g''(%.1f) = (mu - r) - f sigma^2 = %.6f\n', f_test, (MU - R) - f_test * SIGMA^2);
fprintf('  %10s %20s\n', 'h', '(g(f+h) - g(f))/h');
for h = [1e-1 1e-2 1e-3 1e-4 1e-5 1e-6]
    fprintf('  %10.0e %20.10f\n', h, difference_quotient(g, f_test, h));
end

f_star = kelly_fraction(MU, SIGMA, R);
f_star_num = numerical_kelly_fraction(MU, SIGMA, R);
fprintf('\n  Max growth (Kelly):  f* = (mu - r)/sigma^2 = %.6f   numerical (bisection on g'') = %.6f\n', f_star, f_star_num);
fprintf('  growth at f*: g(f*) = %.4f%% per year,  but VaR(f*) = %.2f%%\n', ...
    100 * g(f_star), 100 * value_at_risk(f_star, MU, SIGMA, R, T, ALPHA));

[f_var, hist] = max_fraction_under_var(LOSS_LIMIT, MU, SIGMA, R, T, ALPHA, F_MAX);
fprintf('\n  IVT + bisection on h(f) = VaR(f) - L over [0, %.1f]:\n', F_MAX);
fprintf('  %4s %10s %10s %10s %12s\n', 'k', 'a', 'b', 'mid', 'h(mid)');
for k = 1:size(hist, 1)
    if k <= 8 || k == size(hist, 1)
        fprintf('  %4d %10.6f %10.6f %10.6f %12.2e\n', k, hist(k, 1), hist(k, 2), hist(k, 3), hist(k, 4));
    elseif k == 9
        fprintf('  %4s\n', '...');
    end
end
f_opt = optimal_fraction(MU, SIGMA, R, T, LOSS_LIMIT, ALPHA, F_MAX);
fprintf('  check with the quadratic formula: f = %.6f (the other root is negative, so this is the only crossing, for any horizon)\n', ...
    var_root_formula(LOSS_LIMIT, MU, SIGMA, R, T, ALPHA));
fprintf('\n  f_VaR = %.6f  ->  f_opt = min(f*, f_VaR, f_max) = %.6f\n', f_var, f_opt);
fprintf('  => invest %.1f%% in the stock and %.1f%% risk-free; typical growth %.2f%%/yr, VaR %.2f%%\n', ...
    100 * f_opt, 100 * (1 - f_opt), 100 * g(f_opt), 100 * value_at_risk(f_opt, MU, SIGMA, R, T, ALPHA));
es = expected_shortfall(f_opt, MU, SIGMA, R, T, ALPHA);
f_es = max_fraction_under_es(LOSS_LIMIT, MU, SIGMA, R, T, ALPHA, F_MAX);
fprintf(['  Expected Shortfall: the worst %.0f%% of years lose %.2f%% on average (VaR only says ''more than %.0f%%'');\n' ...
    '  the same IVT + bisection on ES(f) - L gives f_ES = %.6f for ''lose at most %.0f%% on average in the worst %.0f%%''\n'], ...
    100 * (1 - ALPHA), 100 * es, 100 * LOSS_LIMIT, f_es, 100 * LOSS_LIMIT, 100 * (1 - ALPHA));

fprintf('\n  Monte Carlo check (20,000 paths): VaR with rebalancing every dt -> analytic %.2f%%\n', 100 * LOSS_LIMIT);
freq_labels = {'yearly', 'quarterly', 'monthly', 'weekly', 'daily'};
freq_steps = [1 4 12 52 252];
for i = 1:numel(freq_steps)
    W = simulate_wealth(f_opt, MU, SIGMA, R, T, freq_steps(i), 20000);
    fprintf('  %10s: VaR = %.2f%%\n', freq_labels{i}, 100 * empirical_var(W(:, end), 1, ALPHA));
end
fprintf('       daily: ES  = %.2f%%  (formula %.2f%%)\n', 100 * empirical_es(W(:, end), 1, ALPHA), 100 * es);

fig = figure('Position', [100 100 1100 800]);
f_grid = linspace(0, 2 * f_star, 400);

subplot(2, 2, 1);
plot(f_grid, g(f_grid), 'LineWidth', 2); hold on;
xline(f_star, '--', 'Color', RED);
xline(f_opt, '-', 'Color', GREEN);
xlabel('fraction in stock f'); ylabel('growth rate g(f)'); title('Growth: continuous, concave, max at f^*');
legend('g(f)', sprintf('Kelly f^* = %.2f', f_star), sprintf('chosen f_{opt} = %.2f', f_opt), 'Location', 'south');

subplot(2, 2, 2);
var_grid = value_at_risk(f_grid, MU, SIGMA, R, T, ALPHA);
plot(f_grid, var_grid, 'LineWidth', 2); hold on;
yline(LOSS_LIMIT, 'k:');
xline(f_var, '-', 'Color', GREEN);
xlabel('fraction in stock f'); ylabel('1-year VaR'); title('Risk: VaR(f) - L changes sign -> root exists');
legend(sprintf('VaR_{%.0f%%}(f)', 100 * ALPHA), sprintf('limit L = %.0f%%', 100 * LOSS_LIMIT), ...
    sprintf('IVT root f_{VaR} = %.3f', f_var), 'Location', 'southeast');

subplot(2, 2, 3); hold on; box on;
t_grid = linspace(0, T, 253);
band_f = [f_star, f_opt];
band_color = {RED, GREEN};
band_handles = gobjects(1, 2);
for i = 1:2
    W = simulate_wealth(band_f(i), MU, SIGMA, R, T, 252, 5000);
    q = zeros(3, size(W, 2));
    for c = 1:size(W, 2)
        q(:, c) = quantile_linear(W(:, c), [1 - ALPHA, 0.5, ALPHA]);
    end
    fill([t_grid, fliplr(t_grid)], [q(1, :), fliplr(q(3, :))], band_color{i}, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    band_handles(i) = plot(t_grid, q(2, :), 'Color', band_color{i}, 'LineWidth', 1.5);
end
xlabel('time (years)'); ylabel('wealth (start = 1)'); title('Simulated wealth, daily rebalancing');
legend(band_handles, {'Kelly f^*: median and 5-95% band', 'f_{opt}: median and 5-95% band'}, 'Location', 'northwest');

subplot(2, 2, 4);
W = simulate_wealth(f_opt, MU, SIGMA, R, T, 252, 20000);
loss = 1 - W(:, end);
sim_var = empirical_var(W(:, end), 1, ALPHA);
histogram(loss, 80, 'Normalization', 'pdf'); hold on;
xline(LOSS_LIMIT, 'k:', 'LineWidth', 1.5);
xline(sim_var, '--', 'Color', RED, 'LineWidth', 1.5);
xlabel('1-year loss (negative = gain)'); ylabel('density');
title(sprintf('Loss distribution at f_{opt} = %.2f', f_opt));
legend('simulated losses', sprintf('analytic VaR = %.1f%%', 100 * LOSS_LIMIT), ...
    sprintf('simulated VaR = %.1f%%', 100 * sim_var), 'Location', 'northwest');
save_fig(fig, fig_dir, 'fig4_risk_optimization.png');

%% Module 4b - Fat tails: the bell curve plus rare crashes
header('Module 4b - Fat tails: the bell curve plus rare crashes');
lam = CRASH_RATE; J = CRASH_SIZE;
fprintf(['  a crash makes the stock fall %.0f%%; it happens on average %g times a year' ...
    ' (once every %.0f years); the stock''s average return is still %.0f%%\n'], -100 * J, lam, 1 / lam, 100 * MU);

fprintf('\n  Limit: crash coin-flips every dt = 1/n year -> Poisson law of the number of crashes\n');
fprintf('  %14s %40s\n', 'flips/year n', 'max |Binomial(n, lam/n) - Poisson(lam)|');
for n = [1 12 52 252 10000]
    err = max(arrayfun(@(k) abs(binomial_pmf(k, n, lam / n) - poisson_pmf(k, lam)), 0:30));
    fprintf('  %14d %40.2e\n', n, err);
end

fprintf('\n  Risk at the normal answer f_VaR = %.6f:\n', f_var);
fprintf('  %12s %16s %14s\n', '', 'bell curve only', 'with crashes');
for a = [0.95 0.99]
    fprintf('  %12s %15.2f%% %13.2f%%\n', sprintf('VaR %.0f%%', 100 * a), 100 * value_at_risk(f_var, MU, SIGMA, R, T, a), ...
        100 * value_at_risk_jump(f_var, MU, SIGMA, R, T, lam, J, a));
end
fprintf('  %12s %15.2f%% %13.2f%%   (average loss in the worst 5%%)\n', sprintf('ES %.0f%%', 100 * ALPHA), ...
    100 * expected_shortfall(f_var, MU, SIGMA, R, T, ALPHA), 100 * expected_shortfall_jump(f_var, MU, SIGMA, R, T, lam, J, ALPHA));
fprintf('  %12s %15.2f%% %13.2f%%   (same average return, but ln(1 + fJ) < fJ)\n', 'growth', ...
    100 * growth_rate(f_var, MU, SIGMA, R), 100 * growth_rate_jump(f_var, MU, SIGMA, R, lam, J));

[f_crash, hist] = max_fraction_under_var_jump(LOSS_LIMIT, MU, SIGMA, R, T, lam, J, ALPHA, F_MAX);
fprintf('\n  IVT + bisection on h(f) = VaR_crash(f) - L (each VaR_crash is itself a bisection on the CDF):\n');
fprintf('  %4s %10s %10s %10s %12s\n', 'k', 'a', 'b', 'mid', 'h(mid)');
for k = 1:size(hist, 1)
    if k <= 5 || k == size(hist, 1)
        fprintf('  %4d %10.6f %10.6f %10.6f %12.2e\n', k, hist(k, 1), hist(k, 2), hist(k, 3), hist(k, 4));
    elseif k == 6
        fprintf('  %4s\n', '...');
    end
end
f_star_jump = kelly_fraction_jump(MU, SIGMA, R, lam, J);
f_opt_jump = optimal_fraction_jump(MU, SIGMA, R, T, LOSS_LIMIT, lam, J, ALPHA, F_MAX);
fprintf('\n  f_crash = %.6f;  Kelly with crashes (bisection on the slope of the growth): f* = %.4f instead of %.2f\n', ...
    f_crash, f_star_jump, kelly_fraction(MU, SIGMA, R));
fprintf('  => f_opt = min(f*, f_crash, f_max) = %.6f: with crashes, invest %.1f%% in the stock instead of %.1f%%\n', ...
    f_opt_jump, 100 * f_opt_jump, 100 * f_var);

fprintf('\n  Monte Carlo check (40,000 paths, crash coin-flips each step): VaR at f_crash -> %.2f%%\n', 100 * LOSS_LIMIT);
freq_labels = {'yearly', 'monthly', 'weekly', 'daily'};
freq_steps = [1 12 52 252];
for i = 1:numel(freq_steps)
    W = simulate_wealth_jump(f_crash, MU, SIGMA, R, T, lam, J, freq_steps(i), 40000);
    fprintf('  %10s: VaR = %.2f%%\n', freq_labels{i}, 100 * empirical_var(W(:, end), 1, ALPHA));
end

fig = figure('Position', [100 100 1100 420]);
subplot(1, 2, 1);
losses = [0.10 0.20 0.30];          % chance of losing more than this in one year, at f_VaR
p_bell = arrayfun(@(l) jump_log_wealth_cdf(log(1 - l), f_var, MU, SIGMA, R, T, 0, J), losses);
p_crash = arrayfun(@(l) jump_log_wealth_cdf(log(1 - l), f_var, MU, SIGMA, R, T, lam, J), losses);
b = bar(100 * [p_bell; p_crash]');
b(1).FaceColor = [0.12 0.47 0.71]; b(2).FaceColor = RED;
set(gca, 'XTickLabel', arrayfun(@(l) sprintf('lose more than %.0f%%', 100 * l), losses, 'UniformOutput', false));
ylabel('chance in one year (%)');
title(sprintf('Same average return, more big losses (f = %.2f)', f_var));
legend('bell curve only', 'with rare crashes');

subplot(1, 2, 2);
f_grid = linspace(0, 1.5, 61);
plot(f_grid, value_at_risk(f_grid, MU, SIGMA, R, T, ALPHA), 'LineWidth', 2, 'Color', [0.12 0.47 0.71]); hold on;
plot(f_grid, arrayfun(@(v) value_at_risk_jump(v, MU, SIGMA, R, T, lam, J, ALPHA), f_grid), 'LineWidth', 2, 'Color', RED);
yline(LOSS_LIMIT, 'k:');
xline(f_var, '--', 'Color', [0.12 0.47 0.71]);
xline(f_crash, '--', 'Color', RED);
xlabel('fraction in stock f'); ylabel('1-year VaR'); title('Crashes move the IVT root to the left');
legend('bell curve only', 'with rare crashes', sprintf('limit L = %.0f%%', 100 * LOSS_LIMIT), ...
    sprintf('f_{VaR} = %.3f', f_var), sprintf('f_{crash} = %.3f', f_crash), 'Location', 'northwest');
save_fig(fig, fig_dir, 'fig5_fat_tails.png');

%% Module 5 - More assets: the best mix of stocks, gold, bonds and the bank
header('Module 5 - More assets: the best mix of stocks, gold, bonds and the bank');
for i = 1:3
    fprintf('  %7s: mu = %.2f%%, sigma = %.0f%%\n', ASSET_NAMES{i}, 100 * MU_ASSETS(i), 100 * SIGMA_ASSETS(i));
end
fprintf('  bank r = %.0f%%;  rule: with %.0f%% confidence lose at most %.0f%% in one year; no borrowing\n', ...
    100 * R, 100 * ALPHA, 100 * LOSS_LIMIT);

fprintf('\n  One risky asset + bank (Module 4 for each asset, f <= 1):\n');
fprintf('  %7s %10s %10s\n', 'asset', 'f', 'growth');
E = eye(3);
for i = 1:3
    [fi, gi] = fund_fraction(E(i, :), LOSS_LIMIT, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA);
    fprintf('  %7s %10.6f %9.4f%%\n', ASSET_NAMES{i}, fi, 100 * gi);
end

[w5, f5, d5] = best_mix_under_var(LOSS_LIMIT, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA);
g5 = portfolio_growth(w5, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R);
fprintf('\n  Every mix d of the risky assets is one fund; Module 4 (IVT + bisection) gives each fund''s f.\n');
fprintf('  Tried %d funds (1%% steps). Fastest-growing one under the limit:\n', size(mix_grid(3), 1));
fprintf('  stocks %.0f%%, gold %.0f%%, bonds %.0f%%, bank %.0f%%\n', 100 * w5, 100 * (1 - sum(w5)));
fprintf('  growth %.4f%%/yr, VaR %.4f%%\n', 100 * g5, 100 * portfolio_var(w5, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA));
[f4, g4] = fund_fraction([1 0 0], LOSS_LIMIT, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA);
fprintf('  vs stocks + bank (Module 4): growth %.4f%%/yr  ->  +%.2f%% a year for the same risk\n', 100 * g4, 100 * (g5 - g4));
[lo, hi, n_near] = near_optimal_ranges(LOSS_LIMIT, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA, 100, 1e-4);
fprintf(['  The top is flat: %d mixes grow within 0.01%%/yr of the best: stocks %.0f%%-%.0f%%, gold %.0f%%-%.0f%%,' ...
    ' bonds %.0f%%-%.0f%% -> about 73%% stocks is the message; the gold/bond split barely matters\n'], n_near, 100 * [lo; hi]);

fprintf(['\n  Module 5''s limit: grid step -> 0 (each grid is twice as fine and contains the last,' ...
    ' so the best growth never goes down)\n']);
fprintf('  %9s %10s %34s %10s\n', 'step', 'mixes', 'best mix: stocks / gold / bonds', 'growth');
steps_list = [25 50 100 200 400 800 1600];
limit_g = zeros(size(steps_list));
for i = 1:numel(steps_list)
    s = steps_list(i);
    [wk, ~, ~, limit_g(i)] = best_mix_on_grid(LOSS_LIMIT, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA, s);
    fprintf('  %8.4f%% %10d %34s %9.5f%%\n', 100 / s, (s + 1) * (s + 2) / 2, ...
        sprintf('%.2f%% / %.2f%% / %.2f%%', 100 * wk), 100 * limit_g(i));
end
% along the flat top the exact spot is only clear to about 0.1%, so the mix is shown with one decimal
[w_lim, ~, ~, g_lim] = refine_best_mix(LOSS_LIMIT, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA);
fprintf('  %9s %10s %34s %9.5f%%   <- the limit: the exact best mix\n', '-> 0', 'zoom in', ...
    sprintf('%.1f%% / %.1f%% / %.1f%%', 100 * w_lim), 100 * g_lim);

fprintf('\n  Monte Carlo check (40,000 one-year paths, rebalanced daily): VaR of the best mix\n');
W = simulate_wealth_mix(w5, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, 252, 40000);
fprintf('       daily: VaR = %.2f%%  (analytic %.2f%%)\n', 100 * empirical_var(W(:, end), 1, ALPHA), ...
    100 * portfolio_var(w5, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA));

fig = figure('Position', [100 100 1600 440]);
subplot(1, 3, 1);
plans = [f4 0 0 1 - f4; w5 1 - sum(w5)];
b = barh(100 * plans, 'stacked');
colors = [0.12 0.47 0.71; 0.75 0.54 0; GREEN; 0.6 0.6 0.6];
for i = 1:4, b(i).FaceColor = colors(i, :); end
set(gca, 'YTickLabel', {'stocks + bank (Module 4)', 'best mix (Module 5)'});
text(101, 1, sprintf('growth %.2f%%', 100 * g4)); text(101, 2, sprintf('growth %.2f%%', 100 * g5));
xlim([0 130]); xlabel('share of the money (%)');
title(sprintf('Same risk (VaR = %.0f%%), more growth', 100 * LOSS_LIMIT));
legend('stocks', 'gold', 'bonds', 'bank', 'Location', 'northoutside', 'Orientation', 'horizontal');

subplot(1, 3, 2);
x = -log(rand(3000, 4)); x = x ./ sum(x, 2); x = x(:, 1:3);   % random mixes of the 3 assets + bank (Dirichlet)
pts = zeros(3000, 2);
for k = 1:3000
    pts(k, :) = [portfolio_var(x(k, :), MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA), ...
                 portfolio_growth(x(k, :), MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R)];
end
ok = pts(:, 1) <= LOSS_LIMIT;
scatter(pts(~ok, 1), pts(~ok, 2), 6, [0.75 0.75 0.75], 'filled'); hold on;
scatter(pts(ok, 1), pts(ok, 2), 6, GREEN, 'filled');
xline(LOSS_LIMIT, 'k:');
plot(LOSS_LIMIT, g4, '^', 'MarkerSize', 9, 'MarkerFaceColor', [0.12 0.47 0.71], 'MarkerEdgeColor', 'none');
plot(portfolio_var(w5, MU_ASSETS, SIGMA_ASSETS, RHO_ASSETS, R, T, ALPHA), g5, 'p', 'MarkerSize', 16, ...
    'MarkerFaceColor', RED, 'MarkerEdgeColor', 'none');
xlabel('bad-year loss (1-year VaR)'); ylabel('growth rate g');
title('Random mixes: the best is the highest green one');
legend('breaks the limit', 'within the limit', sprintf('limit L = %.0f%%', 100 * LOSS_LIMIT), 'stocks + bank', ...
    'best mix', 'Location', 'northwest');

subplot(1, 3, 3);
plot(1:numel(steps_list), 100 * limit_g, 'o-', 'Color', [0.12 0.47 0.71]); hold on;
yline(100 * g_lim, '--', 'Color', RED);
set(gca, 'XTick', 1:numel(steps_list), ...
    'XTickLabel', arrayfun(@(s) sprintf('%g%%', 100 / s), steps_list, 'UniformOutput', false));
xlim([0.5, numel(steps_list) + 0.5]);
xlabel('grid step (finer to the right)'); ylabel('growth rate g (%)');
title('Module 5''s limit: finer grid -> exact best');
legend('best mix on the grid', sprintf('exact best mix: %.4f%%', 100 * g_lim), 'Location', 'southeast');
save_fig(fig, fig_dir, 'fig6_best_mix.png');

%% Local functions
function header(title_text)
fprintf('\n%s\n%s\n%s\n', repmat('=', 1, 70), title_text, repmat('=', 1, 70));
end

function save_fig(fig, fig_dir, name)
path = fullfile(fig_dir, name);
print(fig, path, '-dpng', '-r150');
fprintf('  saved %s\n', path);
end
