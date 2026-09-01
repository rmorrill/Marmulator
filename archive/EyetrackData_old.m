classdef EyetrackData < handle

    % to-do: when last packet of data comes in, currently it is held,
    % despite any calls to getData. ideally reading data would also make it
    % disappear 

    properties (Access = public)
        eyetracker = []; % for connection handle
        Data = [];
    end

    %properties (Access = public)

    
    methods
        % Constructor to initialize the class
        function obj = EyetrackData()
            %obj.Data = []; % Initialize data
        end

        % Method to add data
        function obj = addData(obj, newData)
            obj.Data = newData; % Overwrite data
        end

        function connectTCPIP(obj, IPaddr, port, varargin)
            %connectTCPIP(obj, IPaddr, port, connect_timeout)
            if numel(varargin)>0
                connect_timeout = varargin{1}; 
            else
                connect_timeout = 10; 
            end
            obj.eyetracker = tcpclient(IPaddr, port, 'ConnectTimeout', connect_timeout);             
            %configureCallback(obj.eyetracker, 'terminator', @obj.readAddData);
            configureCallback(obj.eyetracker, 'terminator', @(src, event) obj.readAddData(src, event)); 
            flush(obj.eyetracker)
            write(obj.eyetracker, 'start'); 

        end

        %function readAddData(src, ~, obj)
        function readAddData(obj, src, ~)
            % callback for tcp client
            data = readline(src); 
            obj.addData(data); 
        end

        function startSend(obj)
            write(obj.eyetracker, 'start')
        end

        function stopSend(obj)
            write(obj.eyetracker, 'stop')
        end

        function delete(obj)
            if ~isempty(obj.eyetracker)
                write(obj.eyetracker, 'stop'); 
            end
            delete(obj.eyetracker); 
            fprintf('eyetracker deleted\n')
        end

        % Method to retrieve data
        function pupil = getData(obj)
            data = obj.Data;
            if isempty(data)
                %fprintf('EyetrackData: data read is empty\n')
                pupil.t_stamp = NaN; 
                pupil.x_pos = NaN; 
                pupil.y_pos = NaN; 
                pupil.radius = NaN; 
                pupil.eye_open_rat = NaN; 
            else
                data_split = split(data, ',');
                data_vec = arrayfun(@str2num, data_split);
                pupil.t_stamp = data_vec(1);
                pupil.x_pos = data_vec(2);
                pupil.y_pos = data_vec(3);
                pupil.radius = data_vec(4);
                pupil.eye_open_rat = data_vec(5);
            end
            obj.Data = [];
        end

    end
end