classdef EyetrackData < handle
% Streaming client for iRecHS2.
%
% Design note: iRecHS2 pushes continuously. MATLAB's tcpclient buffers those
% bytes on a background I/O thread, so the main thread does not need to poll
% at the sample rate -- it only needs to DRAIN the buffer often enough that it
% does not grow without bound. drain() is O(1) in the number of samples waiting
% and performs no parsing, which is what keeps it safe to call from inside the
% Psychtoolbox flip loop.
%
% Byte bookkeeping: rawbuf holds a contiguous SUFFIX of the session stream, not
% the whole thing -- flushToDisk() writes the new bytes out and retains only a
% short tail so parseTail() can still find a complete line. Two counters keep
% that honest across flushes:
%   buf_base  : absolute (1-based) session byte index of rawbuf(1)
%   bytes_out : number of session bytes already written to fid
% so rawbuf(p) is session byte (buf_base + p - 1), and the bytes still owed to
% disk are the local range (bytes_out - buf_base + 2 : wptr - 1). drain_w is
% recorded in absolute session coordinates so that parse_eyetracker_raw can map
% each drain onto a sample index in the parsed file.

    properties (Access = public)
        eyetracker = [];
        Data       = [];        % most recent parsed line (control path)
        hold_last  = true;      % on a frame with no new sample, repeat the last
    end

    properties (Access = private)
        rawbuf     = uint8([]); % raw ASCII capture buffer (suffix of the stream)
        wptr       = 1;
        buf_base   = 1;         % absolute session byte index of rawbuf(1)
        drain_t    = [];        % GetSecs() at each drain
        drain_w    = [];        % absolute session byte index just past each drain
        dptr       = 0;
        capturing  = false;
        fid        = -1;
        capfile    = '';
        bytes_out  = 0;         % bytes already flushed to disk
        last       = struct('t_stamp',NaN,'x_pos',NaN,'y_pos',NaN, ...
                            'radius',NaN,'eye_open_rat',NaN);
    end

    properties (Constant, Access = private)
        BUFCHUNK = 8*2^20;      % 8 MB growth step (~7 min at 20 KB/s)
        DRAINCHUNK = 2^16;
        TAILKEEP = 256;         % bytes retained after a flush for parseTail
    end

    methods
        function obj = EyetrackData()
        end

        function connectTCPIP(obj, IPaddr, port, varargin)
            if numel(varargin) > 0
                connect_timeout = varargin{1};
            else
                connect_timeout = 10;
            end
            obj.eyetracker = tcpclient(IPaddr, port, ...
                                       'ConnectTimeout', connect_timeout);
            configureCallback(obj.eyetracker, 'off');   % polled, not callback
            flush(obj.eyetracker);
            write(obj.eyetracker, 'start');
        end

        % ---------- capture control ----------

        function startCapture(obj, capfile)
        % Begin high-rate capture. Call once, immediately before the trial
        % loop -- AFTER stimulus loading, so startup gaps do not count.
            obj.capfile   = capfile;
            obj.fid       = fopen(capfile, 'w');
            if obj.fid < 0
                error('EyetrackData:capture', ...
                      'could not open capture file %s', capfile);
            end
            obj.rawbuf    = zeros(1, obj.BUFCHUNK, 'uint8');
            obj.wptr      = 1;
            obj.buf_base  = 1;
            obj.drain_t   = zeros(1, 1e6);
            obj.drain_w   = zeros(1, 1e6);
            obj.dptr      = 0;
            obj.bytes_out = 0;
            flush(obj.eyetracker);          % discard pre-session backlog
            obj.capturing = true;
            fprintf('eyetracker capture started -> %s\n', capfile);
        end

        function n = drain(obj)
        % Pull everything buffered. O(1) in samples waiting; no parsing.
        % Safe to call inside the flip loop.
            n = 0;
            if isempty(obj.eyetracker), return, end
            n = obj.eyetracker.NumBytesAvailable;
            if n == 0, return, end

            if obj.capturing
                need = obj.wptr + n - 1;
                if need > numel(obj.rawbuf)
                    grow = max(obj.BUFCHUNK, need - numel(obj.rawbuf));
                    obj.rawbuf(numel(obj.rawbuf) + grow) = uint8(0);
                end
                obj.rawbuf(obj.wptr : obj.wptr+n-1) = read(obj.eyetracker, n, 'uint8');
                obj.wptr = obj.wptr + n;
                obj.dptr = obj.dptr + 1;
                obj.drain_t(obj.dptr) = GetSecs();
                obj.drain_w(obj.dptr) = obj.buf_base + obj.wptr - 1;  % absolute
            else
                % not capturing: keep only a small tail for the control path
                b = read(obj.eyetracker, n, 'uint8');
                keep = min(numel(b), obj.DRAINCHUNK);
                obj.rawbuf = b(end-keep+1:end);
                obj.wptr   = keep + 1;
                obj.buf_base = 1;
            end
        end

        function flushToDisk(obj)
        % Write out the bytes not yet on disk and reset. Call once per trial,
        % in the ITI. Only the range past bytes_out is written -- the tail
        % retained by the previous flush is already in the file, and writing it
        % again would splice a partial line into the middle of the record and
        % stop any downstream sscanf dead.
            if ~obj.capturing || obj.wptr <= 1, return, end

            lo = obj.bytes_out - obj.buf_base + 2;   % first local byte owed
            hi = obj.wptr - 1;                       % last local byte held
            if lo <= hi
                fwrite(obj.fid, obj.rawbuf(lo:hi), 'uint8');
                obj.bytes_out = obj.buf_base + hi - 1;
            end

            % keep a tail so getData can still find a complete line, and move
            % buf_base with it so drain_w stays in session coordinates
            keep = min(hi, obj.TAILKEEP);
            obj.buf_base = obj.buf_base + hi - keep;
            obj.rawbuf(1:keep) = obj.rawbuf(hi-keep+1 : hi);
            obj.wptr = keep + 1;
        end

        function meta = stopCapture(obj)
            if ~obj.capturing
                meta = struct(); return
            end
            obj.drain();
            obj.flushToDisk();
            fclose(obj.fid);
            obj.fid = -1;
            obj.capturing = false;
            meta = struct('capfile',   obj.capfile, ...
                          'drain_t',   obj.drain_t(1:obj.dptr), ...
                          'drain_w',   obj.drain_w(1:obj.dptr), ...
                          'n_bytes',   obj.bytes_out, ...
                          'n_drains',  obj.dptr);
            obj.rawbuf = uint8([]);
            fprintf('eyetracker capture stopped: %d bytes, %d drains\n', ...
                    obj.bytes_out, obj.dptr);
        end

        % ---------- control path (~60 Hz), API unchanged ----------

        function pupil = getData(obj)
            obj.drain();
            pupil = obj.parseTail();
            if isnan(pupil.x_pos) && obj.hold_last
                pupil = obj.last;       % no new sample this frame: hold
            else
                obj.last = pupil;
            end
            obj.Data = pupil;
        end

        function startSend(obj),  write(obj.eyetracker, 'start'); end
        function stopSend(obj),   write(obj.eyetracker, 'stop');  end

        function delete(obj)
            try
                if obj.capturing, obj.stopCapture(); end
            catch
            end
            if ~isempty(obj.eyetracker)
                try write(obj.eyetracker, 'stop'); catch, end
            end
            delete(obj.eyetracker);
            fprintf('eyetracker deleted\n');
        end
    end

    methods (Access = private)
        function pupil = parseTail(obj)
        % Parse ONLY the last complete line in the buffer. sscanf, not str2num.
            pupil = struct('t_stamp',NaN,'x_pos',NaN,'y_pos',NaN, ...
                           'radius',NaN,'eye_open_rat',NaN);
            if obj.wptr <= 2, return, end
            lo   = max(1, obj.wptr - 256);
            tail = obj.rawbuf(lo : obj.wptr-1);
            nl   = find(tail == uint8(10), 2, 'last');   % LF
            if numel(nl) < 2, return, end
            ln = tail(nl(1)+1 : nl(2)-1);
            ln(ln == uint8(13)) = [];                    % strip CR
            v = sscanf(char(ln), '%f,%f,%f,%f,%f');
            if numel(v) < 5, return, end
            pupil.t_stamp      = v(1);
            pupil.x_pos        = v(2);
            pupil.y_pos        = v(3);
            pupil.radius       = v(4);
            pupil.eye_open_rat = v(5);
        end
    end
end
