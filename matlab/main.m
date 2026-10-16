%% Limits in Finance - MATLAB version of main.py
% Runs all four modules, prints the results and saves the figures to
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
% Modules 3-4: stock market and risk-free asset (annual numbers)
MU = 0.10; SIGMA = 0.20; R = 0.03; T = 1;
ALPHA = 0.95; LOSS_LIMIT = 0.15; F_MAX = 2;

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
fprintf('\n  f_VaR = %.6f  ->  f_opt = min(f*, f_VaR, f_max) = %.6f\n', f_var, f_opt);
fprintf('  => invest %.1f%% in the stock and %.1f%% risk-free; growth %.2f%%/yr, VaR %.2f%%\n', ...
    100 * f_opt, 100 * (1 - f_opt), 100 * g(f_opt), 100 * value_at_risk(f_opt, MU, SIGMA, R, T, ALPHA));

fprintf('\n  Monte Carlo check (20,000 paths): VaR with rebalancing every dt -> analytic %.2f%%\n', 100 * LOSS_LIMIT);
freq_labels = {'yearly', 'quarterly', 'monthly', 'weekly', 'daily'};
freq_steps = [1 4 12 52 252];
for i = 1:numel(freq_steps)
    W = simulate_wealth(f_opt, MU, SIGMA, R, T, freq_steps(i), 20000);
    fprintf('  %10s: VaR = %.2f%%\n', freq_labels{i}, 100 * empirical_var(W(:, end), 1, ALPHA));
end

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

%% Local functions
function header(title_text)
fprintf('\n%s\n%s\n%s\n', repmat('=', 1, 70), title_text, repmat('=', 1, 70));
end

function save_fig(fig, fig_dir, name)
path = fullfile(fig_dir, name);
print(fig, path, '-dpng', '-r150');
fprintf('  saved %s\n', path);
end
