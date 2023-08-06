function varargout = Marmulator(varargin)
% MARMULATOR MATLAB code for Marmulator.fig
%      MARMULATOR, by itself, creates a new MARMULATOR or raises the existing
%      singleton*.
%
%      H = MARMULATOR returns the handle to a new MARMULATOR or the handle to
%      the existing singleton*.
%
%      MARMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARMULATOR.M with the given input arguments.
%
%      MARMULATOR('Property','Value',...) creates a new MARMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Marmulator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Marmulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Marmulator

% Last Modified by GUIDE v2.5 05-Aug-2023 18:37:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Marmulator_OpeningFcn, ...
                   'gui_OutputFcn',  @Marmulator_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Marmulator is made visible.
function Marmulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Marmulator (see VARARGIN)

% Choose default command line output for Marmulator
handles.output = hObject;

% first check for setup params
m = mfilename('fullpath'); 
mdir = fileparts(m); 
setup_params_file = fullfile(mdir, 'setup_config.mat'); 
if ~exist(setup_params_file, 'file')
    % complain 
    warndlg(sprintf('No setup config found in %s, edit and run setup_marmulator.m', mdir), 'Setup not found'); 
    handles.sc = []; 
    handles.base_dir = mdir; 
    handles.default_calib_dir = mdir; 
else
    sc = load(setup_params_file); 
    set(handles.com_serial_edit, 'String', sc.serial_pump_comport); 
    set(handles.trig_arduino_com_edit, 'String', sc.arduino_triggers_comport); 
    set(handles.lick_arduino_com_edit, 'String', sc.arduino_lickometer_comport); 
    set(handles.com_edit, 'String', sc.arduino_pump_comport); 
    if isfield(sc,'expt_params_dir')
        set(handles.expt_params_edit, 'String', sc.expt_params_dir);
        handles.expt_dir = sc.expt_params_dir; 
    else
        set(handles.expt_params_edit, 'String', fullfile(sc.marmulator_base_dir, 'experiment_params'));  
    end
    handles.base_dir = sc.marmulator_base_dir;
    handles.default_calib_dir = sc.save_dir_local; 
    if isfield(sc, 'save_dir_local_extra')
        if ~isempty(sc.save_dir_local_extra)
            set(handles.data_save_local_dir_edit, 'Enable', 'on')
            set(handles.data_save_local_dir_edit, 'String', sc.save_dir_local_extra)
        else
            set(handles.data_save_local_dir_edit, 'Enable', 'off')
            set(handles.data_save_local_dir_edit, 'String', '')
        end
    else
         set(handles.data_save_local_dir_edit, 'Enable', 'off')
         set(handles.data_save_local_dir_edit, 'String', '')
    end
    handles.setup_config = sc; 
    handles.reward_list = [get(handles.reward_popup, 'String') sc.reward_types]; 
    set(handles.reward_popup, 'String', handles.reward_list); 
    handles.sc = sc; 
end

logo_file = fullfile(handles.base_dir, 'gui', 'marmulator_logo.png'); 
if isfile(logo_file)
    logo_img = imread(logo_file); 
    axes(handles.logo_axes); 
    imshow(logo_img); 
else
    delete(handles.logo_axes)
end

handles.trig_arduino_connected = 0; 
handles.lick_arduino_connected = 0; 
handles.pump_arduino_connected = 0; 
handles.pump_serial_connected = 0; 
handles.expt_params_dir = get(handles.expt_params_edit, 'String'); 
handles.reward_today = 0; % start reward counter at 0 mL 
handles.calibration_loaded = false; 
handles.calib_file = []; 

handles.subject = []; 
handles.subject_file = []; 
curr_date = datestr(now, 'yyyy-mm-dd');
handles.curr_date = curr_date; 

% ensure that directories exist 
ppdir = fullfile(handles.base_dir, 'gui'); 
if ~exist(ppdir) 
    mkdir(ppdir); 
end

% load the pump parameters
handles.ppfile = fullfile(ppdir, 'pump_params.mat');
if exist(handles.ppfile)
    pp = load(handles.ppfile);
    handles.reward_vol = pp.reward_vol; 
    handles.reward_rate = pp.reward_rate; 
    handles.syringe_diam = pp.syringe_diam; 
    set(handles.vol_edit, 'String', sprintf('%d', handles.reward_vol)); 
    set(handles.rate_edit, 'String', sprintf('%0.1f', handles.reward_rate)); 
    set(handles.diam_edit, 'String', sprintf('%0.2f', handles.syringe_diam)); 
else
    handles.reward_vol = str2double(get(handles.vol_edit, 'String'));
    handles.reward_rate = str2double(get(handles.rate_edit, 'String'));
    handles.syringe_diam = str2double(get(handles.diam_edit, 'String'));
end

% default reward pump duration for arduino
handles.reward_dur = 0.1; % around 10 microliters
set(handles.reward_dur_edit, 'String', sprintf('%0.2f', handles.reward_dur)); 

flist = dir([handles.expt_params_dir '\*.mat']); 
flist = flist(~[flist.isdir]); 

handles.expt_params_list = {flist.name}; 
popup_list_all = ['[select experiment params]', handles.expt_params_list]; 
set(handles.select_popup, 'String', popup_list_all); 

% set close request function
%keyboard
set(handles.figure1, 'CloseRequestFcn', @cleanUpFun); 

%set(handles.expt_params_edit, 'String', fullfile(handles.base_dir, 'experiment_params')); 
set(handles.expt_params_edit, 'String', handles.expt_dir); 
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Marmulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Marmulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in connect_push.
function connect_push_Callback(hObject, eventdata, handles)
% hObject    handle to connect_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.arduino_comport = get(handles.com_edit, 'String');

port = get(handles.com_edit, 'String');
if ~handles.pump_arduino_connected
    baudrate = 115200;
    readtimeout = 0.01; % 10ms 
    %[ahand, errmsg] = IOPort('OpenSerialPort', port, sprintf('BaudRate=%d ReceiveTimeout=0.1', baudrate));
    %configstr = ''; 
    configstr = sprintf('DataBits=8 Parity=None StopBits=1 DTR=1 RTS=1 FlowControl=Hardware HardwareBufferSizes=16,16 BaudRate=%d ReceiveTimeout=0.1', baudrate, readtimeout); 
    [ahand, errmsg] = IOPort('OpenSerialPort', port, configstr);
    if isempty(errmsg) % good, success
        fprintf('Pump arduino on port %s is connected\n', port);
        reward_arduino.ahand = ahand; 
        IOPort('Flush', ahand); 
        set(gcbo, 'String', 'Disconnect');
        set(handles.com_edit, 'enable', 'off');
        handles.pump_arduino_connected = 1;
        handles.reward_vol = 100 *handles.reward_dur; % microliters pump has 1.6ml/min flow rate
        reward_arduino.pump_pin = handles.sc.reward_pin;
        reward_arduino.pump_led_pin = handles.sc.reward_led_pin; 
        assign_pump_pins(reward_arduino); 
        pump_cmd = gen_pump_command(reward_arduino);
        handles.pump_cmd = pump_cmd; 
        handles.reward_arduino = reward_arduino;
        handles.arduino_comport = port;
    else
        fprintf('ERROR CONNECTING TO %s\n', port);
        fprintf('Check:\nis device plugged in?\nis device on %s?\nis device being used by another program?', port);
        handles.pump_arduino_connected = 0;
    end
else
    IOPort('Flush', handles.reward_arduino.ahand); 
    IOPort('Close', handles.reward_arduino.ahand);    
    set(handles.com_edit, 'enable', 'on');
    set(gcbo, 'String', 'Connect');
    handles.pump_arduino_connected = 0;
    fprintf('Disconnecting from port %s\n', port);
    handles.arduino_comport = []; 
end

guidata(hObject, handles);

function com_edit_Callback(hObject, eventdata, handles)
% hObject    handle to com_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of com_edit as text
%        str2double(get(hObject,'String')) returns contents of com_edit as a double


% --- Executes during object creation, after setting all properties.
function com_edit_CreateFcn(hObject, ~, handles)
% hObject    handle to com_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in reward_push.
function reward_push_Callback(hObject, eventdata, handles)
% hObject    handle to reward_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%keyboard
if handles.pump_serial_connected
    writeline(handles.reward_serial, 'RUN');
    reward_given = 1;
else
    disp('No connection to pump');
    reward_given = 0;
end

if reward_given
    curr_date = datestr(now, 'yyyy-mm-dd');
    if strcmp(handles.curr_date, curr_date)
        handles.reward_today = handles.reward_today + handles.reward_vol/1e3;
        set(handles.reward_today_txt, 'String', sprintf('%0.3f mL', handles.reward_today));
    else % leftover, so save amount and update
        if ~isempty(handles.subject_file)
            reward_today = getRewardTodayFromTxt(handles.reward_today_txt);
            updateSubjectLogTable(handles.subject_file, reward_today, handles.curr_date);
        end
        handles.reward_today = 0 + handles.reward_vol/1e3; % reset
        set(handles.reward_today_txt, 'String', sprintf('%0.3f mL', handles.reward_today));
        handles.curr_date = curr_date;
    end
end
guidata(hObject, handles);

% --- Executes on button press in run_push.
function run_push_Callback(hObject, eventdata, handles)
% hObject    handle to run_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'params_file') || isempty(handles.params_file)
    fprintf('No expt params file selected!\n'); 
    return
end

if handles.pump_arduino_connected 
    ra = handles.reward_arduino; 
elseif handles.pump_serial_connected
    ra = handles.reward_serial; 
else
    ra = []; 
end

curr_date = datestr(now, 'yyyy-mm-dd');

if ~isempty(handles.subject_file)
    handles.reward_today = str2double(char(regexp(get(handles.reward_today_txt, 'String'), '\d*\.\d*', 'match')));
    updateSubjectLogTable(handles.subject_file, handles.reward_today, handles.curr_date);
end

if ~strcmp(curr_date, handles.curr_date)
    set(handles.reward_today_txt, 'String', sprintf('%0.3f mL', 0)); 
    handles.curr_date = curr_date; 
end

subj_tmp = get(handles.subject_edit, 'String'); 
if isempty(subj_tmp)
    disp('Please enter a subject name'); 
    return
else
    handles.subject = subj_tmp; 
end

if ~handles.trig_arduino_connected
    handles.trigger_arduino = []; 
end

if ~handles.lick_arduino_connected 
    handles.lick_arduino = []; 
end

gaze_offset_x = str2num(get(handles.offset_x_edit, 'String')); 
gaze_offset_y = str2num(get(handles.offset_y_edit, 'String')); 
gaze_offset =[gaze_offset_x, gaze_offset_y]; 

repeats_per_loc = str2num(get(handles.repeats_edit, 'String')); 
response_time = str2num(get(handles.response_time_edit, 'String'));  
hold_time = str2num(get(handles.hold_time_edit, 'String'));  
trial_time = str2num(get(handles.trial_time_edit, 'String'));  

session_time = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

mouse_for_eye = get(handles.mouse_eye_check, 'Value'); 

require_fix_tr_init = get(handles.require_fix_check, 'Value');  
fixation_to_init = str2double(get(handles.fixation_edit, 'String')); 
time_out_trial_init_s  = str2double(get(handles.time_out_edit, 'String')); 

n_rsvp = str2double(get(handles.n_rsvp_edit, 'String')); 
break_after = str2double(get(handles.break_after_edit, 'String')); 

if isfield(handles, 'gui_lick_timer') && strcmp(handles.gui_lick_timer.Running, 'on')
   stop(handles.gui_lick_timer); 
end

% get selected reward string
reward_list = get(handles.reward_popup, 'String'); 
reward_idx = get(handles.reward_popup, 'Value'); 
handles.reward_selected = reward_list{reward_idx}; 

if ~isempty(ra)
    % check if reward type is selected
    if strcmp(handles.reward_selected, '[select reward type]')
        disp('Pump is connected - please select a reward type');
        return
    end
end
 

training_notes_str = get(handles.notes_edit, 'String'); 
handles.save_data_dir_extra = get(handles.data_save_local_dir_edit, 'String'); 
set(handles.status_text, 'String', sprintf('Session: calib_%s_%s.mat', handles.subject, session_time)); 
[~,calib,save_full] = EyeTracker_Calibrate_gui_fcn(ra, handles.subject,...
    handles.params_file, handles.calib_file, gaze_offset, repeats_per_loc, ...
    response_time, hold_time, trial_time, session_time, mouse_for_eye,...
    require_fix_tr_init, fixation_to_init, time_out_trial_init_s, handles.reward_today_txt,...
    handles.reward_vol/1e3, handles.punish_time , break_after, n_rsvp, ...
    handles.trigger_arduino, handles.lick_arduino, handles.reward_selected,...
    handles.setup_config, training_notes_str, handles.save_data_dir_extra);

if ~isempty(handles.subject_file)
    %handles.reward_today = str2double(char(regexp(get(handles.reward_today_txt, 'String'), '\d*\.\d*', 'match')));
    handles.reward_today = getRewardTodayFromTxt(handles.reward_today_txt); 
    updateSubjectLogTable(handles.subject_file, handles.reward_today, curr_date);
end

if isfield(handles, 'gui_lick_timer') && strcmp(handles.gui_lick_timer.Running, 'off')
    start(handles.gui_lick_timer);
end

guidata(hObject, handles);


n_pts_tot = calib.n_pts_x * calib.n_pts_y; 
if calib.n_completed > 3
    try
        if contains(handles.params.save_params_name,'center_point') || contains(handles.expt_params_name, 'center_point')
            questans = questdlg('Would you like to run center point estimation?', 'Center point', 'Yes', 'No', 'Yes');
            if strcmp(questans, 'Yes')
                [offset_x, offset_y] = get_mean_x_y_pts(save_full);
            else
                fprintf('Will not run center point estimation\nTo run manually use get_mean_x_y_pts.m\n');
            end
        elseif n_pts_tot > 1
            questans = questdlg('Would you like to run multi-point calibration?', 'Run calibration', 'Yes', 'No', 'Yes');
            if strcmp(questans, 'Yes')
                eyetracker_calibration_fcn(save_full)
            else
                fprintf('Will not run multi-point calibration\nTo run manually use eyetracker_calibration_fcn.m\n');
            end
        end
    catch me
        disp('Attempted to run calibration functions, but an error occurred');
    end
end

%f = parfeval(@EyeTracker_Calibrate_gui_fcn, 1, ra, handles.reward_pin, handles.subject,...
%    handles.params_file, handles.calib_file, gaze_offset, repeats_per_loc, ...
%    response_time, hold_time, trial_time, session_time);


function offset_x_edit_Callback(hObject, eventdata, handles)
% hObject    handle to offset_x_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset_x_edit as text
%        str2double(get(hObject,'String')) returns contents of offset_x_edit as a double


% --- Executes during object creation, after setting all properties.
function offset_x_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset_x_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offset_y_edit_Callback(hObject, eventdata, handles)
% hObject    handle to offset_y_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset_y_edit as text
%        str2double(get(hObject,'String')) returns contents of offset_y_edit as a double


% --- Executes during object creation, after setting all properties.
function offset_y_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset_y_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in select_popup.
function select_popup_Callback(hObject, eventdata, handles)
% hObject    handle to select_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_popup
curridx = get(gcbo, 'Value');
flist = get(gcbo, 'String');
handles.expt_params_name = flist{curridx}; 
if curridx == 1
    handles.params_file = [];
    set(handles.hold_time_edit, 'String', '');
    set(handles.response_time_edit, 'String', '');
    set(handles.mode_text, 'String', '');
    set(handles.stim_text, 'String', '');
    set(handles.n_x_y_string, 'String', '');
    set(handles.fixation_edit, 'String', '');
    set(handles.time_out_edit, 'String', '');
    set(handles.require_fix_check, 'Value', 0);
    set(handles.punish_time_edit, 'String', '');
    set(handles.n_rsvp_edit, 'String', '');
    set(handles.break_after_edit, 'String', '');
else
    handles.params_file = fullfile(handles.expt_params_dir, handles.expt_params_name);
    handles.params = load(handles.params_file);
    handles.repeats_per_loc = handles.params.repeats_per_stim;
    handles.nr_pts_x_y = [handles.params.n_pts_x, handles.params.n_pts_y];
    handles.hold_time = handles.params.time_to_reward;
    handles.response_time = handles.params.time_out_after;
    handles.trial_time = handles.params.presentation_time;
    handles.punish_time = handles.params.punish_length_ms; 
    
    if isfield(handles.params, 'require_fix_tr_init')
        handles.require_fix_tr_init = handles.params.require_fix_tr_init;
    else
        handles.require_fix_tr_init = 0;
    end
    if isfield(handles.params, 'fixation_to_init')
        handles.fixation_to_init = handles.params.fixation_to_init;
    else
        handles.fixation_to_init = [];
    end
    if isfield(handles.params, 'time_out_trial_init_s')
        handles.time_out_trial_init_s = handles.params.time_out_trial_init_s;
    else
        handles.time_out_trial_init_s = [];
    end
    
    if isfield(handles.params, 'n_rsvp')
        handles.n_rsvp = handles.params.n_rsvp; 
    else
        handles.n_rsvp = []; 
    end
    
    if isfield(handles.params, 'break_after')
        handles.break_after = handles.params.break_after; 
    else
        handles.break_after = []; 
    end
    
    %keyboard
    set(handles.trial_time_edit, 'String', sprintf('%d', handles.trial_time));
    set(handles.n_x_y_string, 'String', sprintf('[%d, %d]', handles.nr_pts_x_y(1), handles.nr_pts_x_y(2)));
    set(handles.repeats_edit, 'String', sprintf('%d', handles.repeats_per_loc));
    set(handles.hold_time_edit, 'String', sprintf('%d', handles.hold_time));
    set(handles.response_time_edit, 'String', sprintf('%d', handles.response_time));
    if (isfield(handles.params, 'n_rsvp') && handles.params.n_rsvp > 1) || (isfield(handles.params, 'fixation_exp_mode') && handles.params.fixation_exp_mode)
        handles.trial_mode = 'fix_exp_mode';
    else
        handles.trial_mode = handles.params.trial_mode;
    end
    set(handles.mode_text, 'String', handles.trial_mode);
    set(handles.stim_text, 'String', handles.params.stim_mode);
    set(handles.require_fix_check, 'Value', handles.require_fix_tr_init);
    set(handles.fixation_edit, 'String', num2str(handles.fixation_to_init));
    set(handles.time_out_edit, 'String', num2str(handles.time_out_trial_init_s));
    set(handles.punish_time_edit, 'String', num2str(handles.punish_time));
    set(handles.n_rsvp_edit, 'String', num2str(handles.n_rsvp)); 
    set(handles.break_after_edit, 'String', num2str(handles.break_after)); 
    
    if ~handles.require_fix_tr_init
        set(handles.fixation_edit, 'Enable', 'off');
        set(handles.time_out_edit, 'Enable', 'off');
    else
        set(handles.fixation_edit, 'Enable', 'on');
        set(handles.time_out_edit, 'Enable', 'on');
    end
    
    %keyboard
    
    if strcmp(handles.trial_mode, 'foraging')
        set(handles.trial_time_edit,'Enable', 'off');
        set(handles.hold_time_edit, 'Enable', 'on');
        set(handles.response_time_edit, 'Enable', 'on');
        set(handles.n_rsvp_edit, 'Enable', 'off'); 
        set(handles.break_after_edit, 'Enable', 'off'); 
    elseif strcmp(handles.trial_mode, 'trial')
        set(handles.trial_time_edit,'Enable', 'on');
        set(handles.hold_time_edit, 'Enable', 'off');
        set(handles.response_time_edit, 'Enable', 'off');
        set(handles.n_rsvp_edit, 'Enable', 'off'); 
        set(handles.break_after_edit, 'Enable', 'off'); 
    elseif strcmp(handles.trial_mode, 'fix_exp_mode')
        set(handles.trial_time_edit,'Enable', 'on');
        set(handles.hold_time_edit, 'Enable', 'off');
        set(handles.response_time_edit, 'Enable', 'off');
        set(handles.n_rsvp_edit, 'Enable', 'on'); 
        set(handles.break_after_edit, 'Enable', 'on'); 
    end
end

fprintf('loaded params: %s \n', flist{curridx}); 
disp(handles.params); 

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function select_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function expt_params_edit_Callback(hObject, eventdata, handles)
% hObject    handle to expt_params_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of expt_params_edit as text
%        str2double(get(hObject,'String')) returns contents of expt_params_edit as a double

handles.expt_params_dir = get(gcbo, 'String'); 

flist = dir([handles.expt_params_dir '\*.mat']); 
flist = flist(~[flist.isdir]); 

handles.expt_params_list = {flist.name}; 
popup_list_all = ['[select experiment params]', handles.expt_params_list]; 
set(handles.select_popup, 'String', popup_list_all); 

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function expt_params_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expt_params_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function subject_edit_Callback(hObject, eventdata, handles)
% hObject    handle to subject_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subject_edit as text
%        str2double(get(hObject,'String')) returns contents of subject_edit as a double

% load log 

curr_date = {datestr(now, 'yyyy-mm-dd')};
if ~isempty(handles.subject) && ~isempty(handles.subject_file)
    updateSubjectLogTable(handles.subject_file, handles.reward_today, handles.curr_date)
end

new_subject = get(gcbo, 'String'); 

handles.subject = new_subject; 
%handles.reward_total = 0; 
log_dir = fullfile(handles.base_dir, 'subject_logs');
if ~exist(log_dir)
    mkdir(log_dir); 
end
handles.subject_file = fullfile(log_dir, [handles.subject '.mat']); 
if ~exist(handles.subject_file, 'file')
    questans = questdlg(sprintf('No subject log exists for %s. Would you like to make one?', handles.subject),...
        'Create new file'); 
    if strcmp(questans, 'Yes')
        handles.reward_total = 0; 
        handles.subject_log_table = table(curr_date, handles.reward_total, 'VariableNames',...
            {'dates', 'reward_vol'}); 
        subject_log_table = handles.subject_log_table; 
        save(handles.subject_file, 'subject_log_table');
        handles.reward_today = 0; 
    else
        handles.subject_file = []; 
        handles.reward_today = 0; 
    end
else
    %keyboard
    sf = load(handles.subject_file);
    handles.subject_log_table = sf.subject_log_table;
    curridx = find(strcmp(handles.subject_log_table.dates, curr_date));
    if ~isempty(curridx)
        handles.reward_today = handles.subject_log_table.reward_vol(curridx);
    else
        handles.reward_today  = 0;
    end
   % keyboard 
end

handles.curr_date = curr_date; 
set(handles.reward_today_txt, 'String', sprintf('%0.3f mL', handles.reward_today)); 

guidata(hObject, handles); 

   

% --- Executes during object creation, after setting all properties.
function subject_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subject_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function calib_edit_Callback(hObject, eventdata, handles)
% hObject    handle to calib_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calib_edit as text
%        str2double(get(hObject,'String')) returns contents of calib_edit as a double
if isempty(get(gcbo, 'String'))
    set(handles.offset_y_edit, 'Enable', 'on');
    set(handles.offset_x_edit, 'Enable', 'on');
end
guidata(hObject, handles); 



% --- Executes during object creation, after setting all properties.
function calib_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calib_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in browse_calib_push.
function browse_calib_push_Callback(hObject, eventdata, handles)
% hObject    handle to browse_calib_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

start_calib_dir = [];  
if isfield(handles, 'subject') && ~isempty(handles.subject) 
    subj = handles.subject; 
    thisdatestr = datestr(now, 'yyyy-mm-dd'); 
    start_calib_dir = fullfile(handles.default_calib_dir, subj, thisdatestr, 'calibration', 'calib_coeffs'); 
    if ~exist(start_calib_dir, 'dir') 
        start_calib_dir = fullfile(handles.default_calib_dir, subj); 
    end
end

if exist(start_calib_dir, 'dir')
    [fname_tmp, pname_tmp]= uigetfile('*.mat', 'Load calibration file',start_calib_dir);
else
    [fname_tmp, pname_tmp]= uigetfile('*.mat', 'Load calibration file',handles.default_calib_dir);
end

%keyboard
if ~isempty(fname_tmp) && all(fname_tmp ~= 0)
    handles.calib_file = fullfile(pname_tmp, fname_tmp);
    handles.calibration_loaded = true;
    set(handles.calib_edit, 'String', handles.calib_file);
    set(handles.offset_y_edit, 'String', '[]');
    set(handles.offset_x_edit, 'String', '[]');
    set(handles.offset_y_edit, 'Enable', 'off');
    set(handles.offset_x_edit, 'Enable', 'off');
else
    set(handles.offset_y_edit, 'Enable', 'on');
    set(handles.offset_x_edit, 'Enable', 'on');
    handles.calibration_loaded = false;
    set(handles.calib_edit, 'String', ''); 
    handles.calib_file = [];
    fprintf('No calib loaded\n');
end

guidata(hObject, handles);


function repeats_edit_Callback(hObject, eventdata, handles)
% hObject    handle to repeats_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of repeats_edit as text
%        str2double(get(hObject,'String')) returns contents of repeats_edit as a double


% --- Executes during object creation, after setting all properties.
function repeats_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to repeats_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function response_time_edit_Callback(hObject, eventdata, handles)
% hObject    handle to response_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of response_time_edit as text
%        str2double(get(hObject,'String')) returns contents of response_time_edit as a double


% --- Executes during object creation, after setting all properties.
function response_time_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to response_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hold_time_edit_Callback(hObject, eventdata, handles)
% hObject    handle to hold_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hold_time_edit as text
%        str2double(get(hObject,'String')) returns contents of hold_time_edit as a double


% --- Executes during object creation, after setting all properties.
function hold_time_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hold_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function n_x_y_string_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_x_y_string (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trial_time_edit_Callback(hObject, eventdata, handles)
% hObject    handle to trial_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trial_time_edit as text
%        str2double(get(hObject,'String')) returns contents of trial_time_edit as a double


% --- Executes during object creation, after setting all properties.
function trial_time_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trial_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refresh_push.
function refresh_push_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.expt_params_dir = get(handles.expt_params_edit, 'String'); 

flist = dir([handles.expt_params_dir '\*.mat']); 
flist = flist(~[flist.isdir]); 

handles.expt_params_list = {flist.name}; 
popup_list_all = ['[select experiment params]', handles.expt_params_list]; 
set(handles.select_popup, 'String', popup_list_all); 

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in mouse_eye_check.
function mouse_eye_check_Callback(hObject, eventdata, handles)
% hObject    handle to mouse_eye_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mouse_eye_check


% --- Executes on button press in require_fix_check.
function require_fix_check_Callback(hObject, eventdata, handles)
% hObject    handle to require_fix_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of require_fix_check

handles.require_fix_tr_init = get(gcbo, 'Value');

if handles.require_fix_tr_init
    set(handles.fixation_edit, 'Enable', 'on');
    set(handles.time_out_edit, 'Enable', 'on');
else
    set(handles.fixation_edit, 'Enable', 'off');
    set(handles.time_out_edit, 'Enable', 'off');
end
guidata(hObject, handles);



function fixation_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fixation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixation_edit as text
%        str2double(get(hObject,'String')) returns contents of fixation_edit as a double


% --- Executes during object creation, after setting all properties.
function fixation_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function time_out_edit_Callback(hObject, eventdata, handles)
% hObject    handle to time_out_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_out_edit as text
%        str2double(get(hObject,'String')) returns contents of time_out_edit as a double


% --- Executes during object creation, after setting all properties.
function time_out_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_out_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function com_serial_edit_Callback(hObject, eventdata, handles)
% hObject    handle to com_serial_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of com_serial_edit as text
%        str2double(get(hObject,'String')) returns contents of com_serial_edit as a double


% --- Executes during object creation, after setting all properties.
function com_serial_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to com_serial_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in connect_serial_push.
function connect_serial_push_Callback(hObject, eventdata, handles)
% hObject    handle to connect_serial_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.pump_serial_comport = get(handles.com_serial_edit, 'String');
if ~handles.pump_serial_connected
    try
        handles.reward_serial = serialport(handles.pump_serial_comport, 19200); 
        flush(handles.reward_serial);
        %set(gcbo, 'Enable', 'off')
        handles.pump_serial_connected = 1;
        set(gcbo, 'String', 'Disconnect')
        % keyboard
        disp('Serial connected');
        
        sendParamsToPump(handles.reward_serial, handles.reward_vol, ...
            handles.reward_rate, handles.syringe_diam); 

    catch me
        disp(me);
        disp('Serial not connected!');
        handles.pump_serial_connected = 0;
    end
else
    %fid = instrfind('Port', handles.pump_serial_comport);
    %fclose(fid);
    delete(handles.reward_serial);
    handles.reward_serial = [];
    handles.pump_serial_connected = 0;
    set(gcbo, 'String', 'Connect Serial')
    disp('Pump via serial disconnected');
end

guidata(hObject, handles);



function vol_edit_Callback(hObject, eventdata, handles)
% hObject    handle to vol_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vol_edit as text
%        str2double(get(hObject,'String')) returns contents of vol_edit as a double


% --- Executes during object creation, after setting all properties.
function vol_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vol_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rate_edit_Callback(hObject, eventdata, handles)
% hObject    handle to rate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rate_edit as text
%        str2double(get(hObject,'String')) returns contents of rate_edit as a double


% --- Executes during object creation, after setting all properties.
function rate_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diam_edit_Callback(hObject, eventdata, handles)
% hObject    handle to diam_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diam_edit as text
%        str2double(get(hObject,'String')) returns contents of diam_edit as a double


% --- Executes during object creation, after setting all properties.
function diam_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diam_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in send_params_push.
function send_params_push_Callback(hObject, eventdata, handles)
% hObject    handle to send_params_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.pump_serial_connected
    handles.reward_vol = str2double(get(handles.vol_edit, 'String'));
    handles.reward_rate = str2double(get(handles.rate_edit, 'String'));
    handles.syringe_diam = str2double(get(handles.diam_edit, 'String'));
    
    reward_vol = handles.reward_vol;
    reward_rate = handles.reward_rate;
    syringe_diam = handles.syringe_diam;
    
    sendParamsToPump(handles.reward_serial, reward_vol, reward_rate, syringe_diam); 

    lastUpdate = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    
    save(handles.ppfile, 'reward_vol', 'reward_rate', 'syringe_diam', 'lastUpdate');
    disp('Updated pump params file');
    
else
    disp('No serial connection');
end

guidata(hObject, handles);

function sendParamsToPump(pumphand, reward_vol, reward_rate, syringe_diam)

pausetime = 0.1;
configureTerminator(pumphand,'CR')
diam_str = sprintf('DIA %0.2f', syringe_diam);
writeline(pumphand, diam_str);
WaitSecs(pausetime);
rat_str = sprintf('RAT %0.2f MM', reward_rate);
writeline(pumphand, rat_str);
WaitSecs(pausetime);
writeline(pumphand, 'DIR INF');
WaitSecs(pausetime);
vol_units_str = sprintf('VOL UL');
writeline(pumphand, vol_units_str);
WaitSecs(pausetime);
vol_str = sprintf('VOL %0.2f', reward_vol);
writeline(pumphand, vol_str);
WaitSecs(pausetime);
writeline(pumphand, 'LN 1');
WaitSecs(pausetime);
writeline(pumphand, 'TRG T2');
WaitSecs(pausetime);
disp('Sent pump parameters:');
disp(diam_str);
disp(rat_str);
disp(vol_units_str);
disp(vol_str);

function updateSubjectLogTable(subject_file, reward_today_vol, curr_date)

sf = load(subject_file);
subject_log_table = sf.subject_log_table;
curridx = find(strcmp(subject_log_table.dates, curr_date));

if ~isempty(curridx)
    subject_log_table.reward_vol(curridx) = reward_today_vol; 
else
    subject_log_table(end+1,:) = {curr_date, reward_today_vol}; 
end

%subject_log_table
save(subject_file, 'subject_log_table');
fprintf('Updated subject file %s\n', subject_file); 


function reward_today = getRewardTodayFromTxt(reward_today_txt_hand)
reward_today = str2double(char(regexp(get(reward_today_txt_hand, 'String'), '\d*\.\d*', 'match')));
    
function cleanUpFun(hObject, eventdata, handles)

try 
h = guidata(hObject); 
if ~isempty(h.subject_file)
    reward_today = getRewardTodayFromTxt(h.reward_today_txt); 
    updateSubjectLogTable(h.subject_file, reward_today, h.curr_date);
end

if h.pump_arduino_connected
%     fid = instrfind('Port', h.arduino_comport);
%     fclose(fid);
%     delete(h.reward_arduino);
    IOPort('Flush', h.reward_arduino.ahand);
    IOPort('Close',h.reward_arduino.ahand); 
    disp('Pump arduino disconnected');
end

if h.pump_serial_connected
    delete(h.reward_serial);
    disp('Pump via serial disconnected');
end


if h.trig_arduino_connected
    IOPort('Flush', h.trigger_arduino.ahand); 
    IOPort('Close', h.trigger_arduino.ahand);    
    disp('Trigger arduino disconnected');
end

if h.lick_arduino_connected
    IOPort('Flush', h.lick_arduino.ahand);
    IOPort('Close', h.lick_arduino.ahand);
    disp('Lickometer arduino disconnected');
    try
        stop(h.gui_lick_timer);
        delete(h.gui_lick_timer);
    catch
    end
end

catch me
end
delete(hObject);



function punish_time_edit_Callback(hObject, eventdata, handles)
% hObject    handle to punish_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of punish_time_edit as text
%        str2double(get(hObject,'String')) returns contents of punish_time_edit as a double


% --- Executes during object creation, after setting all properties.
function punish_time_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to punish_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n_rsvp_edit_Callback(hObject, eventdata, handles)
% hObject    handle to n_rsvp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_rsvp_edit as text
%        str2double(get(hObject,'String')) returns contents of n_rsvp_edit as a double


% --- Executes during object creation, after setting all properties.
function n_rsvp_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_rsvp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function break_after_edit_Callback(hObject, eventdata, handles)
% hObject    handle to break_after_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of break_after_edit as text
%        str2double(get(hObject,'String')) returns contents of break_after_edit as a double


% --- Executes during object creation, after setting all properties.
function break_after_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to break_after_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trig_arduino_com_edit_Callback(hObject, eventdata, handles)
% hObject    handle to trig_arduino_com_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trig_arduino_com_edit as text
%        str2double(get(hObject,'String')) returns contents of trig_arduino_com_edit as a double


% --- Executes during object creation, after setting all properties.
function trig_arduino_com_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trig_arduino_com_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in connect_trig_arduino_push.
function connect_trig_arduino_push_Callback(hObject, eventdata, handles)
% hObject    handle to connect_trig_arduino_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

port = get(handles.trig_arduino_com_edit, 'String');
if ~handles.trig_arduino_connected
    baudrate = 115200;
    readtimeout = 0.01; % 10ms 
    %[ahand, errmsg] = IOPort('OpenSerialPort', port, sprintf('BaudRate=%d ReceiveTimeout=0.1', baudrate));
    %configstr = ''; 
    configstr = sprintf('DataBits=8 Parity=None StopBits=1 DTR=1 RTS=1 FlowControl=Hardware HardwareBufferSizes=16,16 BaudRate=%d ReceiveTimeout=0.1', baudrate, readtimeout); 
    [ahand, errmsg] = IOPort('OpenSerialPort', port, configstr);
    if isempty(errmsg) % good, success
        fprintf('Trigger arduino on port %s is connected\n', port);
        trigger_arduino.ahand = ahand; 
        IOPort('Flush', ahand); 
        set(gcbo, 'String', 'Disconnect');
        set(handles.trig_arduino_com_edit, 'enable', 'off');
        handles.trig_arduino_connected = 1;
        trigger_arduino.session_pin = handles.sc.session_pin; 
        trigger_arduino.trial_pin = handles.sc.trial_pin; 
        trigger_arduino.stim_pin = handles.sc.stim_pin; 
        if isfield(handles.sc,'sampleCommand_pin') 
            trigger_arduino.sampleCommand_pin = handles.sc.sampleCommand_pin; 
        else
            trigger_arduino.sampleCommand_pin = []; 
        end
        assign_trigger_pins(trigger_arduino); 
        handles.trigger_arduino = trigger_arduino;
        handles.trigger_arduino_comport = port; 
    else
        fprintf('ERROR CONNECTING TO %s\n', port);
        fprintf('Check:\nis device plugged in?\nis device on %s?\nis device being used by another program?', port);
        handles.trig_arduino_connected = 0;
    end
else
    IOPort('Flush', handles.trigger_arduino.ahand); 
    IOPort('Close', handles.trigger_arduino.ahand);    
    set(handles.trig_arduino_com_edit, 'enable', 'on');
    set(gcbo, 'String', 'Connect');
    handles.trig_arduino_connected = 0;
    fprintf('Disconnecting from port %s\n', port);
    handles.trigger_arduino_comport = []; 
end

guidata(hObject, handles);



function lick_arduino_com_edit_Callback(hObject, eventdata, handles)
% hObject    handle to lick_arduino_com_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lick_arduino_com_edit as text
%        str2double(get(hObject,'String')) returns contents of lick_arduino_com_edit as a double


% --- Executes during object creation, after setting all properties.
function lick_arduino_com_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lick_arduino_com_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in connect_lick_arduino_push.
function connect_lick_arduino_push_Callback(hObject, eventdata, handles)
% hObject    handle to connect_lick_arduino_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

port = get(handles.lick_arduino_com_edit, 'String');
if ~handles.lick_arduino_connected
    baudrate = 115200;
    readtimeout = 0.01; % 10ms 
    configstr = sprintf('DataBits=8 Parity=None StopBits=1 DTR=1 RTS=1 FlowControl=Hardware HardwareBufferSizes=16,16 BaudRate=%d ReceiveTimeout=0.1', baudrate, readtimeout); 
    [ahand, errmsg] = IOPort('OpenSerialPort', port, configstr);
    if isempty(errmsg) % good, success
        fprintf('Lickometer arduino on port %s is connected\n', port);
        lick_arduino.ahand = ahand; 
        IOPort('Flush', ahand); 
        set(gcbo, 'String', 'Disconnect');
        set(handles.lick_arduino_com_edit, 'enable', 'off');
        handles.lick_arduino_connected = 1;
        % MAKE THESE PART OF SETUP SOON 
        lick_arduino.lick_pin = 3;
        assign_lickometer_pins(lick_arduino); 
        handles.lick_arduino = lick_arduino;
        handles.lick_arduino_comport = port;
        %%% set up timer for gui display of licks 
        read_lick_cmd = gen_lickometer_command(lick_arduino); 
        lick_box_hand = handles.lick_box; 
        timer_period = 0.05;
        t = timer('TimerFcn', {@checkLickometer_gui, ahand, read_lick_cmd, lick_box_hand}, 'Period', timer_period, 'ExecutionMode', 'fixedDelay',...
            'Name', 'lickometer_gui_timer', 'ErrorFcn', @timerErrorFcn_gui);        
        handles.gui_lick_timer = t;
        start(handles.gui_lick_timer); 
    else
        fprintf('ERROR CONNECTING TO %s\n', port);
        fprintf('Check:\nis device plugged in?\nis device on %s?\nis device being used by another program?', port);
        handles.lick_arduino_connected = 0;
    end
else
    try
        stop(handles.gui_lick_timer);
        delete(handles.gui_lick_timer);
        IOPort('Flush', handles.lick_arduino.ahand);
        IOPort('Close', handles.lick_arduino.ahand);
    catch
    end
    set(handles.lick_arduino_com_edit, 'enable', 'on');
    set(gcbo, 'String', 'Connect');
    handles.lick_arduino_connected = 0;
    fprintf('Disconnecting from port %s\n', port);
    handles.lick_arduino_comport = []; 
end

guidata(hObject, handles);


% --- Executes on selection change in reward_popup.
function reward_popup_Callback(hObject, eventdata, handles)
% hObject    handle to reward_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns reward_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from reward_popup




% --- Executes during object creation, after setting all properties.
function reward_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reward_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function notes_edit_Callback(hObject, eventdata, handles)
% hObject    handle to notes_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of notes_edit as text
%        str2double(get(hObject,'String')) returns contents of notes_edit as a double


% --- Executes during object creation, after setting all properties.
function notes_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to notes_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function reward_dur_edit_Callback(hObject, eventdata, handles)
% hObject    handle to reward_dur_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reward_dur_edit as text
%        str2double(get(hObject,'String')) returns contents of reward_dur_edit as a double


% --- Executes during object creation, after setting all properties.
function reward_dur_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reward_dur_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in update_reward_dur.
function update_reward_dur_Callback(hObject, eventdata, handles)
% hObject    handle to update_reward_dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.pump_arduino_connected
    handles.reward_dur = str2double(get(handles.reward_dur_edit, 'String'));
    fprintf('reward duration updated to %0.2f seconds\n', handles.reward_dur)
    handles.reward_vol = 100 *handles.reward_dur; % microliters
else
    disp('No arduino connection');
end

guidata(hObject, handles);

% --- Executes on button press in reward_arduino_push.
function reward_arduino_push_Callback(hObject, eventdata, handles)
% hObject    handle to reward_arduino_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.pump_arduino_connected
%     IOPort('Write',handles.reward_arduino.ahand,handles.pump_cmd.on, 1);
%     WaitSecs(handles.reward_dur);
%     IOPort('Write',handles.reward_arduino.ahand,handles.pump_cmd.off, 1);

    time = typecast(uint16(handles.reward_dur*1000),'uint8'); %convert to 2 uint8 bytes
    IOPort('Write',handles.reward_arduino.ahand,uint8([53 97+handles.reward_arduino.pump_pin time(1) time(2)]),1);
    reward_given = 1;
else
    disp('No connection to pump');
    reward_given = 0;
end

if reward_given
    curr_date = datestr(now, 'yyyy-mm-dd');
    if strcmp(handles.curr_date, curr_date)
        handles.reward_today = handles.reward_today + handles.reward_vol/1e3;
        set(handles.reward_today_txt, 'String', sprintf('%0.3f mL', handles.reward_today));
    else % leftover, so save amount and update
        if ~isempty(handles.subject_file)
            reward_today = getRewardTodayFromTxt(handles.reward_today_txt);
            updateSubjectLogTable(handles.subject_file, reward_today, handles.curr_date);
        end
        handles.reward_today = 0 + handles.reward_vol/1e3; % reset
        set(handles.reward_today_txt, 'String', sprintf('%0.3f mL', handles.reward_today));
        handles.curr_date = curr_date;
    end
end

guidata(hObject,handles);



function data_save_local_dir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to data_save_local_dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_save_local_dir_edit as text
%        str2double(get(hObject,'String')) returns contents of data_save_local_dir_edit as a double


% --- Executes during object creation, after setting all properties.
function data_save_local_dir_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_save_local_dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in update_data_save_local_dir.
function update_data_save_local_dir_Callback(hObject, eventdata, handles)
% hObject    handle to update_data_save_local_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles.setup_config, 'save_dir_local_extra')
    start_calib_dir = handles.setup_config.save_dir_local_extra; 
    if isfolder(start_calib_dir)
        dname = uigetdir(start_calib_dir);
    else
        dname = '';
    end
    
    set(handles.data_save_local_dir_edit, 'String', dname); 
end


guidata(hObject, handles);
