% convert Marmulator Data to mkTurk data 
% Marmulator can currently run calibration/fixation tasks

function convert_to_mkTurkData(save_dir,filename) 
load(filename) 
defaultFilePath = 'default.json'; % mkTurk datafile 
fid = fopen(defaultFilePath);
raw = fread(fid);
str = char(raw'); 
jsonfile = jsondecode(str);
fclose(fid)
%% TASK   
jsonfile.TASK.Agent = settings.subject; 
jsonfile.TASK.NRSVP = calib_settings.n_rsvp;
%SamplingStrategy": "uniform_without_replacement",
%"NStickyResponse": "",
jsonfile.TASK.NGridPointsX = calib.n_pts_x;
jsonfile.TASK.NGridPointsY = calib.n_pts_y; 
if calib.n_pts_x > 1 && calib.n_pts_y == 1
    spacingPixels = calib.pts(2,1) - calib.pts(1,1);
elseif calib.n_pts_x == 1 && calib.n_pts_y > 1
    spacingPixels = calib.pts(1,2) - calib.pts(1,1); 
elseif calib.n_pts_x > 1 && calib.n_pts_y > 1
    spacingPixels = calib.pts(2) - calib.pts(1);
else
    spacingPixels = 0;
end

jsonfile.TASK.GridSpacingInches = spacingPixels/ settings.viewportPPI; 
% for calibration task, Marmulator generates multiple grid points
if size(calib.pts,1) > 1
    jsonfile.TASK.FixationGridIndex = -1;
    jsonfile.TASK.SampleGridIndex = -1; 
else
    jsonfile.TASK.FixationGridIndex = 1; 
    jsonfile.TASK.SampleGridIndex = 1; 
end

%"TestGridIndex": [],
jsonfile.TASK.FixationSizeInches = 0;
jsonfile.TASK.FixationTimeOut = calib_settings.time_out_trial_init_s * 1000; 
jsonfile.TASK.ImageBagsSample = calib_settings.img_folder; 
%"SamplePRE": 0,
jsonfile.TASK.SampleOFF =  calib_settings.rsvp_iti_ms;
%"ImageBagsTest": [],
jsonfile.TASK.RewardStage =double(reward.give_rewards);
%"NRewardMax": 1,
%"NConsecutiveHitsforBonus": 20000,
jsonfile.TASK.PunishTimeOut = punish.length_ms; 
%"Species": "marmoset",
jsonfile.TASK.RewardDuration = reward.pulse_on_dur * 1000; 
%"Target": "gridwindow",
jsonfile.TASK.Photodiode = double(calib_settings.photodiode_on);
jsonfile.TASK.SampleOutsideGracePeriod =  calib_settings.rsvp_break_after_t; % same as blinkgraceperiod
%"FixationOutsideGracePeriod": 20000,
jsonfile.TASK.FixationWindowSizeInches = calib_settings.bounding_box_x /settings.viewportPPI; 
jsonfile.TASK.FixationDotSizeInches = calib_settings.stim_pre_dot_sz/settings.viewportPPI; 
jsonfile.TASK.FixationDuration = calib_settings.fixation_to_init; 
%"NFixations": 1,
%"FixationUsesSample": 0,
%"SameDifferent": 0,
%"VisualSearch": 0,
%"ObjectGridIndex": [],
%"TestOFF": "",
%"KeepSampleON": 0,
%"ChoiceGridIndex": [],
%"ChoiceSizeInches": "",
%"KeepTestON": "",
%"ChoiceTimeOut": "",
%"HideChoiceDistractors": "",
%"ChoiceOutsideGracePeriod": "",
%"NStimuliPerBagBlock": "",
jsonfile.TASK.BlinkGracePeriod = calib_settings.rsvp_break_after_t; 
jsonfile.TASK.NRSVPMax = calib_settings.n_rsvp;
%"Automator": 0,
%"AutomatorFilePath": "",
%"CurrentAutomatorStage": 0,
jsonfile.TASK.GridXOffsetInches = calib_settings.gaze_center_adj_x / settings.viewportPPI; 
jsonfile.TASK.GridYOffsetInches = calib_settings.gaze_center_adj_y / settings.viewportPPI; 
jsonfile.TASK.BackgroundColor2D = rgb2hex(calib_settings.bg_col *ones(1,3));
% "HeadsupDisplayFraction": 0,
% "THREEJScameraFOV": 45,
% "THREEJScameraZDist": 10,
% "THREEJSRenderRatio": 2,
% "SaveImagesResolution": 0,
% "DeviceConfig": "",
% "BQSaveDisplayTimes": "",
% "BQSaveEye": "",
% "BQSaveTouch": "",
jsonfile.TASK.CalibrateEye = 0; 
% "CalibrateEyeCrossTerms": "",
% "CheckRFID": "",
jsonfile.TASK.InterTrialInterval = calib_settings.inter_stim_interval; % Marmulator draws randomly from a set range of iti

%% ENV
%"ResearcherDisplayName": "",
%"ResearcherEmail": "",
%"ResearcherID": "",
jsonfile.ENV.USBDeviceType = 'microcontroller trigger';
jsonfile.ENV.USBDeviceName = 'Arduino Leonardo'; 
jsonfile.ENV.Subject = settings.subject;
%"AgentRFID": "",
jsonfile.ENV.CurrentDate = datestr(datenum(settings.time_start,'yyyy-mm-dd_HH-MM-SS'),'yyyy-mm-ddTHH:MM:SS'); 
if isfield(settings,'devicePixelRatio') && ~isempty(settings.devicePixelRatio)
    jsonfile.ENV.CanvasRatio = 1/settings.devicePixelRatio; 
    jsonfile.ENV.DevicePixelRatio = settings.devicePixelRatio; 
    jsonfile.ENV.ScreenSizePixels = settings.screenPixels /settings.devicePixelRatio; 
else
    jsonfile.ENV.ScreenSizePixels = settings.screenPixels;
end
%"BackingStoreRatio": 1,

%"THREEJStoInches": "",
jsonfile.ENV.FixationRadius = 0; 
jsonfile.ENV.FixationWindowRadius= 0; 
%"FixationColor": "",
jsonfile.ENV.FixationDotRadius = calib_settings.stim_pre_dot_sz; 
%"FixationDotColor": "white",
%"ChoiceRadius": "",
%"ChoiceColor": "",
if numel(calib.pts) == 2
    jsonfile.ENV.XGridCenter = [calib.pts(:,1)];
    jsonfile.ENV.YGridCenter = [calib.pts(:,2)];
else
    jsonfile.ENV.XGridCenter = calib.pts(:,1);
    jsonfile.ENV.YGridCenter = calib.pts(:,2);
end
jsonfile.ENV.RewardDuration = reward.pulse_on_dur *1000;
jsonfile.ENV.ParamFileName = calib_settings.expt_params; 
%"ParamFileRev": "",
%"ParamFileDate": "",
jsonfile.ENV.DataFileName =  settings.data_file_name; 
%"FirestoreDocRoot": "",
%"CurrentAutomatorStageName": "",
%"MinPercentCriterion": -1,
%"MinTrialsCriterion": -1,
%"StagePctCorrect": -1,
%"StageNTrials": -1,

jsonfile.ENV.NRSVPMin = calib_settings.n_rsvp;
jsonfile.ENV.NRSVPMax = calib_settings.n_rsvp;
%"WebBluetoothAvailable": "",
%"WebUSBAvailable": "",
%"BatteryAPIAvailable": "",
%"OffscreenCanvasAvailable": 0,
jsonfile.ENV.UserAgent = "Psychtoolbox"; 
%"WebAppUrl": "",
%"DeviceType": "desktop",
jsonfile.ENV.DeviceBrand = settings.deviceBrand;
jsonfile.ENV.DeviceName = settings.hostname; 
jsonfile.ENV.DeviceScreenWidth = calib_settings.disp_rect(3); 
jsonfile.ENV.DeviceScreenHeight = calib_settings.disp_rect(4); 
jsonfile.ENV.DeviceGPU= "NVIDIA GeForce RTX 2080 Ti"; 
%"DeviceBrowserName": "",
%"DeviceBrowserVersion": "",
jsonfile.ENV.DeviceOSName = settings.computer_type; 
%"DeviceOSCodename": "",
jsonfile.ENV.DeviceOSVersion = settings.osversion; 
%"DeviceTouchscreen": 0,
jsonfile.ENV.ScreenRatio = settings.screenScale/100; 
jsonfile.ENV.ScreenPhysicalPixels = settings.screenPhysicalPixels;
jsonfile.ENV.ScreenSizeInches = settings.screenInches; 
jsonfile.ENV.ViewportPixels = calib_settings.disp_rect(3:4);
jsonfile.ENV.ViewportPPI = settings.viewportPPI; 
jsonfile.ENV.PhysicalPPI = jsonfile.ENV.ScreenPhysicalPixels(1) / settings.screenInches(1); 
jsonfile.ENV.FrameRateDisplay = 1/settings.ifi_monitor; 
jsonfile.ENV.FrameRateMovie = 1/settings.ifi_monitor; 
%"PrimeScenes": "",
%"NumPrebufferTrials": "",
%"MaxTrialsPerFile": "",
%"Task": "RSVP",
%"MTurkWorkerId": "",
%"AssignmentId": "",
%"HITId": "",
%"StressTest": 0,
if strcmp(eyetrack.method,'mouse')
    jsonfile.ENV.Eye.TrackEye = 0;
else
    jsonfile.ENV.Eye.TrackEey = 1;
end

if size(calib.pts,1) > 1
    jsonfile.ENV.Eye.calibration = 1;
else
    jsonfile.ENV.Eye.calibration = 0; 
end

jsonfile.ENV.Eye.CalibXTransform = eyetrack.cX; 
jsonfile.ENV.Eye.CalibYTransform = eyetrack.cY; 
%"CalibType": "default",
%"NCalibPointsTrain": 0,
%"NCalibPointsTest": 0,
%"CalibTrainMSE": [],
%"CalibTestMSE": [],
%"CalibTestMSETarg": { "x": [], "y": [], "n": [] }
jsonfile.ENV.EffectorSaveJSONDataRelativetoFixationDotDisplayMS = 0; 
jsonfile.ENV.PhotodiodeSquareSizeInches = calib_settings.photodiode_sz / settings.viewportPPI;
jsonfile.ENV.PhotodiodeSquareX = calib_settings.photodiode_pos(1);
jsonfile.ENV.PhotodiodeSquareY = calib_settings.photodiode_pos(2); 
%"RewardSquareSizeInches": "",
jsonfile.ENV.PunishSquareSizeInches = settings.screenInches(1);
%"RewardSquareXY": [],
%"PunishSquareXY": [0, 0],
%"FixationSquareWidth": "",
%"FixationSquareColor": "",
jsonfile.ENV.PhotodiodeSquareWidth = calib_settings.photodiode_sz; 
%"RewardSquareWidth": "",
%"PunishSquareWidth": 1536

%% CANVAS
% 
% "CANVAS": {
% "sequenceblank": [],
% "tsequenceblank": [],
% "sequencepre": [],
% "tsequencepre": [],
% "sequencepost": [],
% "tsequencepost": [],
% "headsupfraction": 0,
% "offsetleft": 0,
% "offsettop": 0,
jsonfile.CANVAS.workspace = [0, 0,calib_settings.disp_rect(3),calib_settings.disp_rect(4)]; 
%% SCENEMETA

%  "SCENEMETA": {
%     "THREEJStoPixels": "",
%     "SampleXGridCenterTHREEJS": [],
%     "SampleYGridCenterTHREEJS": [],
%     "TestXGridCenterTHREEJS": [],
%     "TestYGridCenterTHREEJS": [],
%     "SampleImageSetDir": "",
%     "SampleNouns": [],
%     "SampleObjects": [],
%     "SampleBagNames": [],
jsonfile.SCENEMETA.SampleBagIdx = calib_settings.img_folder_idx-1; %0-index
%     "SampleImageIdx": [],
jsonfile.SCENEMETA.SampleImagePathList = calib_settings.img_list;  
%     "TestImageSetDir": "",
%     "TestNouns": [],
%     "TestObjects": [],
%     "TestBagNames": [],
%     "TestBagIdx": [],
%     "TestImageIdx": []
%   },


%% SCENES
%"SCENES": { "SampleScenes": {}, "TestScenes": {} },
EXT = {'.json'};
if iscell(calib_settings.img_folder)
    for i = 1:numel(calib_settings.img_folder)
        img_fnames_tmp = dir(calib_settings.img_folder{i});
        [~,~,e] = fileparts({img_fnames_tmp.name}); 
        img_fnames_tmp = img_fnames_tmp(matches(e,EXT,'IgnoreCase',1)); 
        if ~isempty(img_fnames_tmp)
            fid = fopen(fullfile(calib_settings.img_folder{i}, img_fnames_tmp.name));
            raw = fread(fid);
            str = char(raw'); 
            scenefile = jsondecode(str);
            jsonfile.SCENES.SampleScenes.(strcat('x',num2str(i-1))) = scenefile;
        end
    end
else
    img_fnames_tmp = dir(calib_settings.img_folder);
    [~,~,e] = fileparts({img_fnames_tmp.name}); 
    img_fnames_tmp = img_fnames_tmp(matches(e,EXT,'IgnoreCase',1)); 
    if ~isempty(img_fnames_tmp)
        fid = fopen(fullfile(calib_settings.img_folder, img_fnames_tmp.name));
        raw = fread(fid);
        str = char(raw'); 
        scenefile = jsondecode(str);
        jsonfile.SCENES.SampleScenes.x0 = scenefile;
    end
end

%% "TRIALEVENTS": {
for i = 1:size(calib_settings.img_seq,1)
    jsonfile.TRIALEVENTS.Sample.(strcat('x',num2str(i-1))) = calib_settings.img_seq(i,1:calib.n_completed)-1; 
end

%"Test": {},
jsonfile.TRIALEVENTS.CorrectItem = calib.sequence(1:calib.n_completed)-1; 
jsonfile.TRIALEVENTS.FixationGridIndex = calib.sequence(1:calib.n_completed)-1;
jsonfile.TRIALEVENTS.StartTime =  calib.stim_pre_start(1:calib.n_completed) *1000; % convert to ms  
jsonfile.TRIALEVENTS.ReinforcementTime = reward.reward_time(1:calib.n_completed) *1000; % convert to ms 
jsonfile.TRIALEVENTS.EndTime = calib.end_t(1:calib.n_completed)*1000; % convert to ms after reward or punishment is delivered 

jsonfile.TRIALEVENTS.FixationTouchEvent = cell(1,calib.n_completed);
for i = 1: calib.n_completed
    if calib.trial_init_timed_out(i) == 1
        jsonfile.TRIALEVENTS.FixationTouchEvent{i} = 'broke outside'; 
    else
        jsonfile.TRIALEVENTS.FixationTouchEvent{i} = 'held'; 
    end
end

jsonfile.TRIALEVENTS.SampleStartTime = calib.start_t(1:calib.n_completed)*1000; % convert to ms  ; 
timedout_tr_ind = find(calib.trial_init_timed_out ==1);
jsonfile.TRIALEVENTS.SampleStartTime(timedout_tr_ind) = -1; 

% most recent fixation 
jsonfile.TRIALEVENTS.FixationXYT.x0 = nan * ones(1,calib.n_completed); 
jsonfile.TRIALEVENTS.FixationXYT.x1= nan * ones(1,calib.n_completed); 
jsonfile.TRIALEVENTS.FixationXYT.x2 = nan * ones(1,calib.n_completed); 

for i = 1:calib.n_completed
    t_ind = find(eyetrack.time >= calib.stim_pre_start(i) & eyetrack.time <=calib.stim_pre_end(i));
    jsonfile.TRIALEVENTS.FixationXYT.x0(i) = eyetrack.x(t_ind(end));
    jsonfile.TRIALEVENTS.FixationXYT.x1(i) = eyetrack.y(t_ind(end));
    jsonfile.TRIALEVENTS.FixationXYT.x2(i) = eyetrack.time(t_ind(end)) *1000; % convert to ms  
end
%"Response": [],

for i = 1:size(calib_settings.clip_sequence,2)
    jsonfile.TRIALEVENTS.TSequenceDesiredClip.(strcat('x',num2str(i-1))) = []; 
    jsonfile.TRIALEVENTS.TSequenceActualClip.(strcat('x',num2str(i-1))) = (calib.clip_sequence_start_t(1:calib.n_completed,i) *1000)'; % convert to ms '; 
end

jsonfile.TRIALEVENTS.SampleGridIndex = ones(1,calib.n_completed);

jsonfile.TRIALEVENTS.SampleFixationTouchEvent = cell(1,calib.n_completed); 
for i = 1:calib.n_completed
    if reward.correct_trial(i) == 1
        jsonfile.TRIALEVENTS.SampleFixationTouchEvent{i} = 'held';
    else
        jsonfile.TRIALEVENTS.SampleFixationTouchEvent{i} = 'broke_early';
    end
end

% most recent fixation during sample presentation per trial 
% group eyetrack data by trial

jsonfile.TRIALEVENTS.SampleFixationXYT.x0 = nan * ones(1,calib.n_completed); 
jsonfile.TRIALEVENTS.SampleFixationXYT.x1= nan * ones(1,calib.n_completed); 
jsonfile.TRIALEVENTS.SampleFixationXYT.x2 = nan * ones(1,calib.n_completed); 

for i = 1:calib.n_completed
    t_ind = find(eyetrack.time >= calib.start_t(i) & eyetrack.time <=calib.end_t(i));
    if ~isempty(t_ind)
        jsonfile.TRIALEVENTS.SampleFixationXYT.x0(i) = eyetrack.x(t_ind(end));
        jsonfile.TRIALEVENTS.SampleFixationXYT.x1(i) = eyetrack.y(t_ind(end));
        jsonfile.TRIALEVENTS.SampleFixationXYT.x2(i) = eyetrack.time(t_ind(end)) *1000; % convert to ms 
    end
end

jsonfile.TRIALEVENTS.ResponseTouchEvent = jsonfile.TRIALEVENTS.SampleFixationTouchEvent; 
jsonfile.TRIALEVENTS.ResponseXYT.x2 = jsonfile.TRIALEVENTS.SampleFixationXYT.x2; 

jsonfile.TRIALEVENTS.NReward = double(reward.correct_trial(1:calib.n_completed)); 
 %"SampleCommandReturnTime": [],
 %"SampleCommandOffReturnTime": []

 %% TIMEEVENTS
 if ~strcmp(eyetrack.method,'mouse')
     jsonfile.TIMEEVENTS.EffectorXY.t = struct; 
     jsonfile.TIMEEVENTS.EffectorXY.x = struct; 
     jsonfile.TIMEEVENTS.EffectorXY.y = struct; 
     jsonfile.TIMEEVENTS.EffectorXY.w = struct; 
     jsonfile.TIMEEVENTS.EffectorXY.a = struct; 
     jsonfile.TIMEEVENTS.EffectorXY.q = struct; 
     for i = 1:calib.n_completed
        t_ind = find(eyetrack.time >= calib.start_t(i) & eyetrack.time <=calib.end_t(i));
        if ~isempty(t_ind)
            jsonfile.TIMEEVENTS.EffectorXY.t.(strcat('x',num2str(i-1))) = eyetrack.time(t_ind) * 1000; 
            jsonfile.TIMEEVENTS.EffectorXY.x.(strcat('x',num2str(i-1))) = eyetrack.x(t_ind); 
            jsonfile.TIMEEVENTS.EffectorXY.y.(strcat('x',num2str(i-1))) = eyetrack.y(t_ind); 
            if numel(eyetrack.pupil_size_x) >1
                jsonfile.TIMEEVENTS.EffectorXY.w.(strcat('x',num2str(i-1))) = eyetrack.pupil_size_x(t_ind); 
                jsonfile.TIMEEVENTS.EffectorXY.a.(strcat('x',num2str(i-1))) = eyetrack.pupil_size_x(t_ind)./eyetrack.pupil_size_y(t_ind); 
            end
            jsonfile.TIMEEVENTS.EffectorXY.q.(strcat('x',num2str(i-1))) = eyetrack.quality(t_ind);
        end
     end
 end
 
 jsonfile.TIMEEVENTS.DiplayTimes = struct;
 jsonfile.TIMEEVENTS.DisplayTimes.t_actual = calib.frame_st_t(1:calib.n_completed);
 jsonfile.TIMEEVENTS.DisplayTimes.frame_num = calib.frame_clip_num(1:calib.n_completed); 
 
%% save
save_file_date = strrep(jsonfile.ENV.CurrentDate,':','_'); 
save_file_name = append(save_file_date,'_', settings.subject,'.json'); 

jsonfile_to_save = jsonencode(jsonfile,PrettyPrint=true); 
fid = fopen(fullfile(save_dir,save_file_name), 'w');
fprintf(fid, '%s', jsonfile_to_save);
fclose(fid);
end
