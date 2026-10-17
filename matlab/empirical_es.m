function v = empirical_es(final_wealth, W0, alpha)
%EMPIRICAL_ES ES estimated from simulated final wealth: 1 - the average of the worst 1 - alpha share.
if nargin < 2, W0 = 1; end
if nargin < 3, alpha = 0.95; end
x = sort(final_wealth(:) / W0);
k = max(1, floor((1 - alpha) * numel(x) + 0.5));
v = 1 - mean(x(1:k));
end
