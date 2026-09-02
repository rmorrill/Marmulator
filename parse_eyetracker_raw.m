function H = parse_eyetracker_raw(meta, eye_rect)
% Parse the raw capture into the eyetrack.hires struct.
%   meta     : struct returned by EyetrackData.stopCapture()
%   eye_rect : the same eye_rect engine.m uses for the control path. Eye
%              positions are divided by eye_rect(3)/eye_rect(4) so that
%              hires.eyepos_raw is on the identical 0-1 scale as
%              eyetrack.eyepos_raw. Omit it and the values stay in the
%              eyetracker's own pixel units (and a warning is issued).
%
% This runs at session end, after Screen('CloseAll') but BEFORE the save, so it
% must not throw -- an error here costs the whole session .mat. Every field is
% given a default up front and every degenerate case warns rather than errors.

if nargin < 2, eye_rect = []; end

H = struct('t_stamp', [], 't_ptb', [], 'eyepos_raw', [], 'eye_rect', eye_rect, ...
           'radius', [], 'eye_open_rat', [], ...
           'clock_a', NaN, 'clock_b', NaN, 'clock_stats', struct(), ...
           'capfile', '', 'drain_t', [], 'drain_sample', [], ...
           'n_samples', 0, 'median_dt', NaN, 'sample_rate', NaN, 'n_gaps', NaN);

if ~isstruct(meta) || ~isfield(meta, 'capfile') || isempty(meta.capfile)
    warning('parse_eyetracker_raw:noMeta', ...
            'no capture metadata (capture never started?); eyetrack.hires is empty');
    return
end
H.capfile = meta.capfile;
H.drain_t = meta.drain_t;
if ~isfile(meta.capfile)
    warning('parse_eyetracker_raw:noFile', ...
            'capture file not found: %s', meta.capfile);
    return
end

raw = fileread(meta.capfile);
raw(raw == char(13)) = [];                       % CRLF -> LF

% Drop a trailing partial line. sscanf with size [5 Inf] zero-pads an incomplete
% final column rather than complaining, which would invent a bogus sample.
if ~isempty(raw) && raw(end) ~= newline
    lastnl = find(raw == newline, 1, 'last');
    if isempty(lastnl)
        raw = '';
    else
        raw = raw(1:lastnl);
    end
end

V = sscanf(raw, '%f,%f,%f,%f,%f', [5 Inf])';     % whole file, one call
if isempty(V)
    warning('parse_eyetracker_raw:noSamples', ...
            'eyetracker capture parsed to zero samples: %s', meta.capfile);
    return
end

% byte offset -> sample index, so each drain can be tied to a sample.
% line_end and drain_w are both ascending, so merge the two sorted lists
% instead of the obvious arrayfun(@(w) sum(line_end < w), ...): that is
% O(n_drains * n_lines), which on an hour-long session (~200k drains x ~400k
% lines) is ~1e11 comparisons and hangs the end of the session for minutes.
nl        = find(raw == newline);
line_end  = nl(:);
dw        = meta.drain_w(:);
[dws, isort] = sort(dw);                       % robust to unsorted drain_w
n_le      = numel(line_end);
[~, ord]  = sort([line_end; dws - 0.5]);       % -0.5 makes ties resolve as "<"
is_line   = ord <= n_le;
c         = cumsum(is_line);
drain_sam = zeros(numel(dw), 1);
drain_sam(isort) = c(~is_line);
drain_sam = reshape(drain_sam, size(meta.drain_w));

good = drain_sam >= 1 & drain_sam <= size(V,1);

[a, b, cstats] = fit_eyetracker_clock(V(drain_sam(good),1), meta.drain_t(good));

H.t_stamp      = V(:,1);
H.t_ptb        = a * V(:,1) + b;      % absolute GetSecs; subtract
                                      % settings.t_start_sec to compare with
                                      % eyetrack.time and the calib_* times
% Same normalization the control path applies at engine.m:2712-2713, so that
% hires.eyepos_raw and eyetrack.eyepos_raw are directly comparable. Still
% uncalibrated: apply cX/cY/cProj offline to get screen coordinates.
if numel(eye_rect) >= 4
    H.eyepos_raw = [V(:,2)/eye_rect(3), V(:,3)/eye_rect(4)];
else
    warning('parse_eyetracker_raw:noEyeRect', ...
            ['eye_rect not supplied; hires.eyepos_raw is left in eyetracker ' ...
             'pixel units and will NOT match eyetrack.eyepos_raw']);
    H.eyepos_raw = [V(:,2), V(:,3)];
end
H.radius       = V(:,4);
H.eye_open_rat = V(:,5);
H.clock_a      = a;
H.clock_b      = b;
H.clock_stats  = cstats;
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
