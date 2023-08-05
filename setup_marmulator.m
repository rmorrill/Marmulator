% run this once, or anytime configuration params change! 

%%%% Set these parameters: 
marmulator_base_dir = 'C:\Users\issalab\Documents\GitHub\Marmulator'; 
expt_params_dir = fullfile(marmulator_base_dir,'experiment_params'); 
save_dir_local = 'C:\Users\issalab\Documents\GitHub\Marmulator\calib_data';
save_dir_local_extra = 'D:\YJ'; 
save_dir_remote = 'Z:\Data'; %'\\locker-smb.engram.rc.zi.columbia.edu\Data'; 
log_dir = 'C:\Users\issalab\Documents\GitHub\Marmulator\subject_logs_session'; 
eyetracker_toolbox_dir = 'C:\Users\issalab\Desktop\EYETRACKER\Interfaces\3rdParty\Microsoft-Windows\MATLAB\ViewPoint_EyeTracker_Toolbox'; 

serial_pump_comport = ''; 
arduino_pump_comport = 'COM10'; 
arduino_lickometer_comport = ''; 
arduino_triggers_comport = 'COM5'; 

% arduino pin configuration
session_pin = 4; 
trial_pin = 7; 
stim_pin = 10; 
sampleCommand_pin = 12; % set to [] if you don't want it
reward_pin = 2; % 

default_gaze_center_adjust = [0,0];

screenid_stim = 0; 
screenid_ctrl = 0; % ResolutionTest
audio_deviceNum  = 1; % PsychPortAudio('GetDevices')

deviceBrand = 'UperfectMonitor'; 

screenInches = [10.9,6.1,12.5]; % [27.7,15.5]; 
screenPhysicalPixels = [3840,2560]; 
% screnPixels and screenScale are set in the Windows Setting 
screenPixels = [3840, 2560]; %
screenScale = 150; % percentage 
devicePixelRatio = 1; 
if ~isempty(devicePixelRatio) && isempty(screenPixels)
    viewportPPI = screenPhysicalPixels(1)/devicePixelRatio/(screenScale/100)/screenInches(1);
else
    viewportPPI = screenPixels(1)/(screenScale/100)/screenInches(1); 
end

dist_to_screen = 13 /2.54; %cm to inches
% PRINT PIXEL VALUE EQUIVALENT OF 1 visual degree 
deg_to_inch_on_screen = tan(0.5*pi/180)*dist_to_screen*2; %inches
deg_to_pixel_on_screen = deg_to_inch_on_screen * viewportPPI; 


reward_types = {'coconut_milk', 'apple_juice', 'water', 'ensure_diluted', 'evaporated_milk','ensure-yogurt_diluted'}; 


%%%% DO NOT MODIFY %%%%
setup_date = datestr(now, 'yyyy-mm-dd_HH-MM_SS'); 
setup_save_path = fullfile(marmulator_base_dir, 'setup_config.mat'); 

save(setup_save_path, 'marmulator_base_dir', 'save_dir_local','log_dir','expt_params_dir',...
    'save_dir_remote', 'save_dir_local_extra', 'eyetracker_toolbox_dir', 'serial_pump_comport',...
    'session_pin', 'trial_pin', 'stim_pin', 'sampleCommand_pin', ...
    'reward_pin', 'audio_deviceNum',...
    'arduino_pump_comport', 'arduino_lickometer_comport', 'arduino_triggers_comport',...
    'default_gaze_center_adjust', 'screenid_stim', 'screenid_ctrl', ...
    'setup_date', 'reward_types',...
    'deviceBrand','screenInches','screenPixels','screenPhysicalPixels','screenScale','viewportPPI','devicePixelRatio',...
    'dist_to_screen','deg_to_inch_on_screen','deg_to_pixel_on_screen');  

fprintf('saved setup config to %s\n', setup_save_path); 







