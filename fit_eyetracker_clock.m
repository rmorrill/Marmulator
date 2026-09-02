function [a, b, stats] = fit_eyetracker_clock(t_stamp_at_drain, drain_t)
% Robust map  t_ptb ~= a*t_stamp + b, fit to the minimum-latency envelope.
%
%   t_stamp_at_drain : t_stamp of the LAST sample in each drain
%   drain_t          : GetSecs() at that drain
%
% Latency = drain_t - (a*t_stamp + b) is >= 0 by construction, so the correct
% line rides the lower edge of the point cloud, not through its middle.

% Force both to columns before combining. parse_eyetracker_raw passes a column
% (V(idx,1)) and a row (meta.drain_t(good)); combining those with & triggers
% implicit expansion into an NxN mask, which then overruns the array on index.
x = t_stamp_at_drain(:);
y = drain_t(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);

stats = struct('drift_ppm', NaN, 'min_latency_s', NaN, 'med_latency_s', NaN, ...
               'p95_latency_s', NaN, 'n_pairs', numel(x));

if numel(x) < 2
    if isempty(x)
        a = NaN; b = NaN;
    else
        a = 1; b = y - x;          % single pair: offset only, no drift term
        stats.drift_ppm = 0;
        stats.min_latency_s = 0; stats.med_latency_s = 0; stats.p95_latency_s = 0;
    end
    warning('fit_eyetracker_clock:tooFewPairs', ...
            ['only %d usable drain pair(s); clock map is %s. ' ...
             'Check that the capture file and drain index are both intact.'], ...
            numel(x), ternary(isempty(x), 'undefined (NaN)', 'offset-only'));
    return
end

% Center before fitting. t_stamp may be a large absolute value (uptime or epoch
% seconds) spanning only an hour, and an uncentered polyfit on that is badly
% conditioned -- the slope, which IS the drift estimate, is the part that goes.
x0 = x(1); y0 = y(1);
xc = x - x0; yc = y - y0;

p = polyfit(xc, yc, 1);                  % seed: ordinary LS
for it = 1:5
    r   = yc - polyval(p, xc);           % residual = latency + noise
    thr = prctile(r, 15);                % keep the fastest 15%
    keep = r <= thr;
    if nnz(keep) < 10, break, end
    p = polyfit(xc(keep), yc(keep), 1);
end
a = p(1);
b = y0 - a*x0 + p(2);                    % undo the centering

r = yc - polyval(p, xc);
stats.drift_ppm     = (a - 1) * 1e6;     % clock rate mismatch
stats.min_latency_s = min(r);
stats.med_latency_s = median(r);
stats.p95_latency_s = prctile(r, 95);

fprintf('clock fit: drift %.1f ppm | latency min %.2f ms, med %.2f ms, p95 %.2f ms\n', ...
    stats.drift_ppm, stats.min_latency_s*1e3, ...
    stats.med_latency_s*1e3, stats.p95_latency_s*1e3);
if abs(stats.drift_ppm) > 500
    warning(['clock drift %.0f ppm is large; check that t_stamp is in ' ...
        'seconds and not ms or frame counts'], stats.drift_ppm);
end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
