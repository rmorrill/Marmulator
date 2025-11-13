% run this once, or anytime configuration params change! 

%%%% Set these parameters: 
marmulator_base_dir = '/home/ryan/Documents/MATLAB/Marmulator'; 
save_dir_local = '/home/ryan/Data/Marmulator'; 
save_dir_remote = '/mnt/hopper_freiwald/rmorrill/Marmulator_data'; %'\\locker-smb.engram.rc.zi.columbia.edu\Data'; 
eyetracker_toolbox_dir = ''; 

serial_pump_comport = '/dev/ttyUSB0'; 
arduino_pump_comport = ''; 
arduino_lickometer_comport = '/dev/ttyACM1'; 
arduino_triggers_comport = '/dev/ttyACM0'; 

% arduino pin configuration
session_pin = 4; 
trial_pin = 7; 
stim_pin = 10; 
sampleCommand_pin = []; % set to [] if you don't want it

%default_gaze_center_adjust = [0, 210]; % [x, y]
default_gaze_center_adjust = [0,0]; % positive y moves down, negative y moves up

screenid_stim = 1; 
screenid_ctrl = 0; 

deviceBrand = 'ASUS'; 

screenInches = [10.9,6.1,12.5]; % [27.7,15.5]; 
screenPhysicalPixels = [1920,1080]; 
% screnPixels and screenScale are set in the Windows Setting 
screenPixels = [1920, 1080]; %
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

reward_types = {'NONE', 'maple_diluted', 'coconut_milk', 'apple_juice', 'water', 'ensure_diluted', 'evaporated_milk','ensure-yogurt_diluted'}; 


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







