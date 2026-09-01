function H = parse_eyetracker_raw(meta)
% Parse the raw capture into the eyetrack.hires struct.
%   meta : struct returned by EyetrackData.stopCapture()

raw = fileread(meta.capfile);
raw(raw == char(13)) = [];                       % CRLF -> LF

V = sscanf(raw, '%f,%f,%f,%f,%f', [5 Inf])';     % whole file, one call
if isempty(V)
    warning('eyetracker capture parsed to zero samples: %s', meta.capfile);
end

% byte offset -> sample index, so each drain can be tied to a sample
nl        = find(raw == newline);
line_end  = nl(:);
drain_sam = arrayfun(@(w) sum(line_end < w), meta.drain_w);
good      = drain_sam >= 1 & drain_sam <= size(V,1);

[a, b, cstats] = fit_eyetracker_clock(V(drain_sam(good),1), meta.drain_t(good));

H = struct();
H.t_stamp      = V(:,1);
H.t_ptb        = a * V(:,1) + b;      % <-- align to calib.start_t etc.
H.x_raw        = V(:,2);
H.y_raw        = V(:,3);
H.radius       = V(:,4);
H.eye_open_rat = V(:,5);
H.clock_a      = a;
H.clock_b      = b;
H.clock_stats  = cstats;
H.capfile      = meta.capfile;
H.drain_t      = meta.drain_t;
H.drain_sample = drain_sam;
H.n_samples    = size(V,1);
if H.n_samples > 2
    H.median_dt   = median(diff(V(:,1)));
    H.sample_rate = 1 / H.median_dt;
    H.n_gaps      = sum(diff(V(:,1)) > 3 * H.median_dt);
end

fprintf('parsed %d hi-rate samples (%.1f Hz nominal, %d gaps)\n', ...
    H.n_samples, H.sample_rate, H.n_gaps);
end