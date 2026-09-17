% run this once, or anytime configuration params change! 

%%%% Set these parameters: 
marmulator_base_dir = ''; % e.g. '/home/user/Documents/MATLAB/Marmulator' 
save_dir_local = ''; % e.g. '/home/user/Data/Marmulator'  
save_dir_remote = ''; % e.g. '/mnt/server/user/Marmulator_data'
log_dir = fullfile(save_dir_local, 'subject_logs'); % directory for saving subject logs to track reward volume over days

eyetracker_IP = ''; % 'xxx.xx.x.xxx' as a str
eyetracker_port = ''; % 'xxxxx'as a str 

% SETUP PERIPHERALS - ARDUINOS, PUMP
serial_pump_comport = ''; % e.g. '/dev/ttyUSB0'
arduino_lickometer_comport = ''; % e.g. '/dev/ttyUSB1'
arduino_triggers_comport = ''; % e.g. '/dev/ttyACM0'

% TRIGGER ARDUINO PIN CONFIGURATION 
session_pin = 4; 
trial_pin = 7; 
stim_pin = 10; 
sampleCommand_pin = []; % set to [] if you don't want it

% SETUP AUDIO DEVICE
audio_device_keyword = ''; % e.g. 'Steinberg' or 'Rubix22'


% SCREEN PARAMS 
default_gaze_center_adjust = [0,0]; % in pixels from center, negative moves 
screenid_stim = 1; 
screenid_ctrl = 0; 



% NOTE BELOW IS NOT USED - SCREEN SIZE ETC IS RETURNED IN ENGINE.M
% SEE window_rect IN EXPT PARAMS FILE
% IF window_rect = [] THEN WINDOW WILL BE FULL STIM SCREEN
deviceBrand = ''; % e.g. 'ASUS'
screenInches = []; % [27.7,15.5]; 
screenPhysicalPixels = [];  % e.g.  [1920, 1080]
% screenPixels and screenScale are set in the Windows Setting 
screenPixels = []; % % e.g. [1920, 1080], could be different from physical pixels
screenScale = 100; % percentage 
devicePixelRatio = 1; 
%if ~isempty(devicePixelRatio) && isempty(screenPixels)
%    viewportPPI = screenPhysicalPixels(1)/devicePixelRatio/(screenScale/100)/screenInches(1);
%else
%    viewportPPI = screenPixels(1)/(screenScale/100)/screenInches(1); 
%end
viewportPPI = [] % DEPRECATED
dist_to_screen = 13 /2.54; %cm to inches 
deg_to_inch_on_screen = tan(0.5*pi/180)*dist_to_screen*2; %inches
deg_to_pixel_on_screen = deg_to_inch_on_screen * viewportPPI; 
%%%%


reward_types = {''}; % cell of strings for reward


%%%% DO NOT MODIFY %%%%
setup_date = datestr(now, 'yyyy-mm-dd_HH-MM_SS'); 
setup_save_path = fullfile(marmulator_base_dir, 'setup_config_TEST.mat'); 

save(setup_save_path, 'marmulator_base_dir', 'save_dir_local',...
    'save_dir_remote', 'log_dir', 'eyetracker_IP', 'eyetracker_port', 'serial_pump_comport',...
    'session_pin', 'trial_pin', 'stim_pin', 'sampleCommand_pin', ...
    'arduino_lickometer_comport', 'arduino_triggers_comport',...
    'default_gaze_center_adjust', 'screenid_stim', 'screenid_ctrl', ...
    'setup_date', 'reward_types',...
    'deviceBrand','screenInches','screenPixels','screenPhysicalPixels','screenScale','viewportPPI','devicePixelRatio',...
    'dist_to_screen','deg_to_inch_on_screen','deg_to_pixel_on_screen', ...
    'audio_device_keyword');  

fprintf('saved setup config to %s\n', setup_save_path); 







