function [a, b, stats] = fit_eyetracker_clock(t_stamp_at_drain, drain_t)
% Robust map  t_ptb ≈ a*t_stamp + b, fit to the minimum-latency envelope.
%
%   t_stamp_at_drain : t_stamp of the LAST sample in each drain
%   drain_t          : GetSecs() at that drain
%
% Latency = drain_t - (a*t_stamp + b) is >= 0 by construction, so the correct
% line rides the lower edge of the point cloud, not through its middle.

ok = isfinite(t_stamp_at_drain) & isfinite(drain_t);
x = t_stamp_at_drain(ok(:)); y = drain_t(ok(:));

p = polyfit(x, y, 1);                    % seed: ordinary LS
for it = 1:5
    r   = y - polyval(p, x);             % residual = latency + noise
    thr = prctile(r, 15);                % keep the fastest 15%
    keep = r <= thr;
    if nnz(keep) < 10, break, end
    p = polyfit(x(keep), y(keep), 1);
end
a = p(1); b = p(2);

r = y - polyval(p, x);
stats = struct( ...
    'drift_ppm',     (a - 1) * 1e6, ...      % clock rate mismatch
    'min_latency_s', min(r), ...
    'med_latency_s', median(r), ...
    'p95_latency_s', prctile(r, 95), ...
    'n_pairs',       numel(x));

fprintf('clock fit: drift %.1f ppm | latency min %.2f ms, med %.2f ms, p95 %.2f ms\n', ...
    stats.drift_ppm, stats.min_latency_s*1e3, ...
    stats.med_latency_s*1e3, stats.p95_latency_s*1e3);
if abs(stats.drift_ppm) > 500
    warning(['clock drift %.0f ppm is large; check that t_stamp is in ' ...
        'seconds and not ms or frame counts'], stats.drift_ppm);
end
end