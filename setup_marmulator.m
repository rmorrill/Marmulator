% run this once, or anytime configuration params change! 

%%%% Set these parameters: 
marmulator_base_dir = 'C:\MATLAB\Marmulator'; 
save_dir_local = 'C:\MATLAB\eyetracker_calibration_071222\calib_data'; 
save_dir_remote = 'Z:\Data'; %'\\locker-smb.engram.rc.zi.columbia.edu\Data'; 
eyetracker_toolbox_dir = 'C:\Users\Hector\Desktop\eyetracker\Interfaces\3rdParty\Microsoft-Windows\MATLAB\ViewPoint_EyeTracker_Toolbox'; 

serial_pump_comport = 'COM9'; 
arduino_pump_comport = ''; 
arduino_lickometer_comport = ''; 
arduino_triggers_comport = ''; 

default_gaze_center_adjust = [0, 210]; % [x, y]
default_gaze_center_adjust = [0,116]; %2022-10-25 yj

screenid_stim = 1; 
screenid_ctrl = 2; 

%%%% DO NOT MODIFY %%%%
setup_date = datestr(now, 'yyyy-mm-dd_HH-MM_SS'); 
setup_save_path = fullfile(marmulator_base_dir, 'setup_config.mat'); 

save(setup_save_path, 'marmulator_base_dir', 'save_dir_local',...
    'save_dir_remote', 'eyetracker_toolbox_dir', 'serial_pump_comport',...
    'arduino_pump_comport', 'arduino_lickometer_comport', 'arduino_triggers_comport',...
    'default_gaze_center_adjust', 'screenid_stim', 'screenid_ctrl', ...
    'setup_date');  

fprintf('saved setup config to %s\n', setup_save_path); 







