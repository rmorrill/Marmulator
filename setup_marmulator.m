% run this once, or anytime configuration params change! 

%%%% Set these parameters: 
marmulator_base_dir = 'C:\MATLAB\Marmulator'; 
save_dir_local = 'C:\MATLAB\eyetracker_calibration_071222\calib_data'; 
save_dir_remote = 'Z:\Data'; %'\\locker-smb.engram.rc.zi.columbia.edu\Data'; 
eyetracker_toolbox_dir = 'C:\Users\Hector\Desktop\eyetracker\Interfaces\3rdParty\Microsoft-Windows\MATLAB\ViewPoint_EyeTracker_Toolbox'; 

serial_pump_comport = 'COM9'; 
arduino_pump_comport = ''; 
arduino_lickometer_comport = 'COM13'; 
arduino_triggers_comport = 'COM8'; 

% arduino pin configuration
session_pin = 4; 
trial_pin = 7; 
stim_pin = 10; 
sampleCommand_pin = 12; % set to [] if you don't want it

default_gaze_center_adjust = [0, 210]; % [x, y]
default_gaze_center_adjust = [0,116]; %2022-10-25 yj
default_gaze_center_adjust = [0,70]; 
default_gaze_center_adjust = [0,116]; 
default_gaze_center_adjust = [0,90]; % 2022-12-03 yj 
default_gaze_center_adjust = [0,45]; % 2022-12-03 yj 

screenid_stim = 1; 
screenid_ctrl = 2; 

deviceBrand = 'UperfectMonitor'; 

screenInches = [10.9,6.1,12.5]; % [27.7,15.5]; 
screenPhysicalPixels = [3840,2160]; 
% screnPixels and screenScale are set in the Windows Setting 
screenPixels = [2560, 1440]; %
screenScale = 100; % percentage 
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

save(setup_save_path, 'marmulator_base_dir', 'save_dir_local',...
    'save_dir_remote', 'eyetracker_toolbox_dir', 'serial_pump_comport',...
    'session_pin', 'trial_pin', 'stim_pin', 'sampleCommand_pin', ...
    'arduino_pump_comport', 'arduino_lickometer_comport', 'arduino_triggers_comport',...
    'default_gaze_center_adjust', 'screenid_stim', 'screenid_ctrl', ...
    'setup_date', 'reward_types',...
    'deviceBrand','screenInches','screenPixels','screenPhysicalPixels','screenScale','viewportPPI','devicePixelRatio',...
    'dist_to_screen','deg_to_inch_on_screen','deg_to_pixel_on_screen');  

fprintf('saved setup config to %s\n', setup_save_path); 







