
function [eyetrack, calib, save_full] = EyeTracker_Calibrate_gui_fcn(reward_pumphand, reward_arduino_pin,...
    subject, expt_params, calib_fname, gaze_offset, trs_per_location, time_out_after, ...
    time_to_reward, presentation_time, session_time, eye_method_mouse,...
    require_fix_tr_init, fixation_to_init, time_out_trial_init_s, ...
    reward_today_hand, reward_vol, punish_length_ms, rsvp_break_after_t, n_rsvp, ...
    trigger_arduino, lick_arduino, reward_type, setup_config, training_notes_str)

profile_memory = false; % flag for tracking memory usage
% if true, will place mem_used, avail_sys_mem, avail_phys_mem into base
% workspace for debugging

newPriority = 1;
oldPriority = Priority(newPriority);
fprintf('PTB old priority: %d, new priority %d\n', oldPriority, newPriority);
PsychJavaSwingCleanup;

% mkTurk_data_save_flag
mkTurk_data_save = 0;

if ~isempty(calib_fname)
    c = load(calib_fname);
    if isfield(c,'coeff_Proj')
        cProj = c.coeff_Proj; % 3 x 3 matrix
    else
        cProj = [];
    end
    cX = c.coeff_X;
    cY = c.coeff_Y;
    apply_calib = 1;
else
    cY = nan(1,3);
    cX = nan(1,3);
    cProj = nan(1,3);
    apply_calib = 0;
end

if ~isempty(gaze_offset)
    gaze_offset_x = gaze_offset(1);
    gaze_offset_y = gaze_offset(2);
else
    gaze_offset_x = [];
    gaze_offset_y = [];
end

% check reward "arduino" - serial or arduino
if ~isempty(reward_pumphand) && contains(class(reward_pumphand), 'Serialport')
    reward_serial = true;
elseif ~isempty(reward_pumphand) && contains(class(reward_pumphand), 'arduino')
    reward_serial = false;
end

if ~isempty(trigger_arduino)
    trig_hand = trigger_arduino.ahand;
    trig_flag = true;
    [sess_trig_cmd, trial_trig_cmd, stim_trig_cmd,sampleCommand_trig_cmd] = gen_trig_commands(trigger_arduino);
    % ensure that everything is off
    % turn off session trigger
    IOPort('Write', trig_hand, sess_trig_cmd.off, 1);
    WaitSecs(0.05);
    IOPort('Write', trig_hand, trial_trig_cmd.off, 1);
    WaitSecs(0.05);
    IOPort('Write', trig_hand, stim_trig_cmd.off, 1);
    WaitSecs(0.05);
    if ~isempty(sampleCommand_trig_cmd)
        IOPort('Write', trig_hand, sampleCommand_trig_cmd.off, 1);
        WaitSecs(0.05);
    end
else
    trig_flag = false;
    fprintf('WARNING: NO TRIGGER ARDUINO HANDLE PROVIDED!\nTriggers will not be recorded\n');
end

if ~isempty(lick_arduino)
    lick_hand = lick_arduino.ahand;
    lick_flag = true;
    rew_trs_since_lick = NaN;
    lick_times_all = nan(1,1e5);
    lick_idx = 0;
    Alphabet = 'abcdefghijklmnopqrstuvwxyz';
    read_lick_cmd = ['1' Alphabet(lick_arduino.lick_pin + 1)];
    lick_arduino_pin = lick_arduino.lick_pin;
else
    rew_trs_since_lick = NaN;
    lick_flag = false;
    lick_times_all = [];
    lick_arduino_pin = [];
    fprintf('WARNING: NO LICKOMETER ARDUINO HANDLE PROVIDED!\nLicks will not be recorded\n');
end

%% photodiode flash
photodiode_flash = true;
flash_rect_size = [90 90];


%% SETUP: load settings
s = load(expt_params);
%ctrl_screen = false;
ctrl_screen = true;
multiflip = 0;
multiflip2 = 0;
dontclear = 0;
dontclear2 = 0;
dontsync = 0;
dontsync2 = 1;

% unpack everything
save_params_name = s.save_params_name;
save_params_here = s.save_params_here;
eyetracker_toolbox_dir = setup_config.eyetracker_toolbox_dir;
window_rect = s.window_rect;
skip_sync_tests = s.skip_sync_tests;
screenid_stim = setup_config.screenid_stim;
screenid_ctrl = setup_config.screenid_ctrl;
calibration_win_len = s.calibration_win_len;
calibration_win_ht = s.calibration_win_ht;
n_pts_x = s.n_pts_x;
n_pts_y = s.n_pts_y;
inter_stim_interval = s.inter_stim_interval;
iti_random = s.iti_random;
position_order = s.order;
if isfield(s,'image_order')
    image_order = s.image_order;
else
    image_order = 'random';
end

bg_col = s.bg_col;
stim_rect_size_x = s.stim_rect_size_x;
stim_rect_size_y = s.stim_rect_size_y;
show_curr_bounding_box = s.show_curr_bounding_box;
show_all_bounding_boxes = s.show_all_bounding_boxes;
bounding_rect_size_x = s.bounding_rect_size_x;
bounding_rect_size_y = s.bounding_rect_size_y;
manual_bounding_boxes = s.manual_bounding_boxes;
if isfield(s,'bounding_rect_size_stim_pre_dot_x') && ~isempty(s.bounding_rect_size_stim_pre_dot_x)
    bounding_rect_size_stim_pre_dot_x = s.bounding_rect_size_stim_pre_dot_x;
else
    bounding_rect_size_stim_pre_dot_x = bounding_rect_size_x;
end
if isfield(s,'bounding_rect_size_stim_pre_dot_y') && ~isempty(s.bounding_rect_size_stim_pre_dot_y)
    bounding_rect_size_stim_pre_dot_y = s.bounding_rect_size_stim_pre_dot_y;
else
    bounding_rect_size_stim_pre_dot_y = bounding_rect_size_y;
end
if isfield(s,'manual_bounding_boxes stim_pre_dot') && ~isempty(s.manual_bounding_boxes_stim_pre_dot)
    manual_bounding_boxes_stim_pre_dot = s.manual_bounding_boxes_stim_pre_dot;
else
    manual_bounding_boxes_stim_pre_dot = manual_bounding_boxes;
end
trial_mode = s.trial_mode;
reward_on = s.reward_on;
location_require_quality = s.location_require_quality;
stim_mode = s.stim_mode;
img_folder = s.img_folder;
pulse_size = s.pulse_size;
movie_folder = s.movie_folder;
movie_rate = s.movie_rate;
scale_fact_move = s.scale_fact_move;
scale_fact_size = s.scale_fact_size;
mvdot_sz = s.mvdot_sz;
reverse_rate = s.reverse_rate;
reverse_rate_sz = s.reverse_rate_sz;
draw_gaze_pt = s.draw_gaze_pt;
n_gaze_pts_draw = s.n_gaze_pts_draw;
gaze_pt_dot_col = s.dot_col;
gaze_pt_sz = s.gaze_pt_sz;
draw_retain_bb_pts = s.draw_retain_bb_pts;
color_shift_feedback = s.color_shift_feedback;
rew_col_start = s.rew_col_start;
rew_col_end = s.rew_col_end;
play_reward_sound = s.play_reward_sound;
reward_sound_file = s.reward_sound_file;
give_punishments = s.give_punishments;
play_punish_sound = s.play_punish_sound;
punish_sound_file = s.punish_sound_file;
if eye_method_mouse
    fprintf('EYE METHOD MOUSE CHECKED - control position through mouse cursor\n');
    eye_method = 'mouse';
else
    eye_method = s.eye_method;
end
n_frames_blink = s.n_frames_blink;
give_rewards = s.give_rewards;
reward_thresh = s.reward_thresh;
reward_on_dur = s.reward_on_dur;
wait_after_reward = s.wait_after_reward;

if isfield(s, 'stimulus_pre_dot')
    % stimulus pre dot
    stimulus_pre_dot = s.stimulus_pre_dot;
    stimulus_pre_time = s.stimulus_pre_time; %ms
    stim_pre_dot_sz = s.stim_pre_dot_sz;
else
    stimulus_pre_dot = true;
    stimulus_pre_time = 500; %ms
    stim_pre_dot_sz = 10;
end

if isfield(s,'stimulus_pre_dot_disappear')
    stimulus_pre_dot_disappear = s.stimulus_pre_dot_disappear; % if 1, disappear when stimulus is presented
else
    stimulus_pre_dot_disappear = 0;
end

if isfield(s, 'stop_rewards_n_nolicks')
    stop_rewards_n_nolicks = s.stop_rewards_n_nolicks;
else
    stop_rewards_n_nolicks = 4;
end

if strcmp(stim_mode, 'match_to_sample')
    m2s_mode = true; % match to sample mode 
    %search_grid = s.search_grid;
    test_bb_size_x = s.test_bb_size_x; 
    test_bb_size_y = s.test_bb_size_y; 
    test_img_len = s.test_img_len; 
    test_img_ht = s.test_img_ht
    dots_on_test_grid = s.dots_on_test_grid; 
    task_folder = s.task_folder; 
    if isfield(s, 'task_params_file')
        task_params_file = s.task_params_file; 
    else
        task_params_file = ''; 
    end
    
    if isfield(s, 'punish_incorrect_fix') 
        punish_incorrect_fix = s.punish_incorrect_fix; 
    else
        punish_incorrect_fix = true;
    end
else
    m2s_mode = false; 
end


if isfield(s, 'draw_crosshairs') % draw white crosshairs on fixation point instead of standard dot 
    draw_crosshairs = s.draw_crosshairs; 
    crosshair_sz = s.crosshair_sz; 
    crosshair_lw = s.crosshair_lw; 
else
    draw_crosshairs = false; 
    crosshair_sz = []; 
    crosshair_lw = []; 
end
    
%% fixation mode/rsvp setup
% fix mode will use presentation time for each stimulus duration
if isfield(s, 'rsvp_iti_t')
    rsvp_iti_t = s.rsvp_iti_t;
else
    rsvp_iti_t = 0;
end

if isempty(n_rsvp) || isnan(n_rsvp)
    n_rsvp = 1;
end

if isfield(s, 'fixation_exp_mode')
    fix_exp_mode = s.fixation_exp_mode;
else
    if n_rsvp>1
        fix_exp_mode = true;
    else
        fix_exp_mode = false;
    end
end

if isfield(s, 'constant_bg_mode')
    constant_bg_mode = s.constant_bg_mode;
else
    constant_bg_mode = false;
end

if give_punishments
    punish_length_ms_draw = punish_length_ms;
else
    punish_length_ms_draw = NaN;
end
fprintf('punish_length_ms is %d ms\n', punish_length_ms); 

gaze_center_adj_x = setup_config.default_gaze_center_adjust(1);
gaze_center_adj_y = setup_config.default_gaze_center_adjust(2);
apply_gaze_center_adj = true;


bonus_reward_min_wait = 1; % give bonus rewards at most every 1s

start_rew_vol_str = get(reward_today_hand, 'String');
start_reward_vol = str2double(char(regexp(start_rew_vol_str, '\d*\.\d*', 'match')));

expt_type = 'calibration';
curr_date = datestr(now, 'yyyy-mm-dd');
save_base_local = setup_config.save_dir_local;
save_base_remote = setup_config.save_dir_remote;

save_data_dir = fullfile(save_base_local, subject, curr_date, expt_type);
if ~exist(save_data_dir, 'dir')
    mkdir(save_data_dir);
end

if ~isempty(save_base_remote)
    save_data_dir_remote = fullfile(save_base_remote, subject, curr_date, expt_type);
    if ~exist(save_data_dir_remote, 'dir')
        mkdir(save_data_dir_remote);
    end
end

%% SETUP: keystroke shortcuts

esc_key = KbName('ESC');
rew_end_key = KbName('r');
bonus_rew_key = KbName('b');
wake_up_movie_key = KbName('v');
wake_up_movie_terminate_key = KbName('t');

%% SETUP: eyetracker-related
eye_rect = [0 0 1 1] ;
eye = 0; % eye A

good_pts = [];
good_pts_cols = [];

if strcmp(eye_method, 'rand') | strcmp(eye_method, 'mouse')
    use_vpx = false;
else
    use_vpx = true;
end

if use_vpx
    % ensure that vpx eyetracker library is loaded
    if libisloaded('vpx') == 0
        addpath(genpath(eyetracker_toolbox_dir));
        vpx_Initialize;
    end
    
    switch eye_method
        case 'pupil'
            [eyepos_x, eyepos_y] = vpx_GetPupilPoint(eye);
        case 'pupil-glint'
            [eyepos_x, eyepos_y] = vpx_GetDiffVector(eye);
        case 'gaze_point'
            [eyepos_x, eyepos_y] = vpx_GetGazePoint(eye);
        case 'gaze_point_corrected'
            [eyepos_x, eyepos_y] = vpx_GetGazePointCorrected(eye);
    end
    calib_rect_cmd = sprintf('Calibration_RealRect %0.2f %0.2f %0.2f %0.2f', deal(eye_rect));
    vpx_SendCommandString( calib_rect_cmd )
    eyepos_raw = [eyepos_x, eyepos_y];
elseif strcmp(eye_method, 'rand')
    eyepos_x = rand();
    eyepos_y = rand();
    eyepos_raw = rand(1,2);
elseif strcmp(eye_method, 'mouse')
    [eyepos_x, eyepos_y] = GetMouse();
    eyepos_raw = rand(1,2);
    location_require_quality = false;
end

% initialize shared variables
eyetracker_time = [];
eyetracker_qual = [];
pupil_size_x = [];
pupil_size_y = [];

idx_all = 0;
gaze_pt_dot_cols = round(repmat(gaze_pt_dot_col', [1,n_gaze_pts_draw])*255);
alphas = linspace(0, 255, n_gaze_pts_draw);
rgba_cols = [gaze_pt_dot_cols; alphas];
dotsz = gaze_pt_sz;

%% SETUP: image presentation mode

if strcmp(stim_mode, 'images') || strcmp(stim_mode,'spinning') || strcmp(stim_mode, 'smooth pursuit')
    
    [img_fnames, img_folder_idx] = ret_image_fnames(img_folder);
    
    if constant_bg_mode % then exclude background stimuli from the list
        [~,imnames_tmp] = fileparts(img_fnames);
        rmidx = contains(imnames_tmp, 'background_');
        img_fnames = img_fnames(~rmidx);
        img_folder_idx = img_folder_idx(~rmidx);
    end
    
    [imgs, im_ht, im_wd, ar] = load_images(img_fnames);
    n_imgs = numel(imgs);
    
    [unique_img_folder,~,ic] = unique(img_folder_idx);
    n_imgs_per_folder = accumarray(ic,1);
else
    img_folder_idx = []; 
end

%% SETUP: constant background mode (new 2/2023)
if constant_bg_mode % read the backgrounds separately
    % go through folders and search for a 'background_*.*' images in each
    if iscell(img_folder)
        bg_img_fnames = [];
        for i = 1:numel(img_folder)
            bg_img_fnames_tmp = dir(fullfile(img_folder{i}, 'background_*'));
            if numel(bg_img_fnames_tmp)~= 1
                error('number of possible backgrounds in folder %s does not equal 1\nlabel background images as "background_*"', img_folder{i});
            end
            [~,~,e] = fileparts({bg_img_fnames_tmp.name});
            bg_img_fnames_tmp = bg_img_fnames_tmp(matches(e,EXT,'IgnoreCase',1));
            bg_img_fnames_tmp = fullfile(img_folder{i}, {bg_img_fnames_tmp.name});
            bg_img_fnames = [bg_img_fnames bg_img_fnames_tmp];
        end
    else
        bg_img_fnames_tmp = dir(fullfile(img_folder, 'background_*'));
        [~,~,e] = fileparts({bg_img_fnames_tmp.name});
        bg_img_fnames_tmp = bg_img_fnames_tmp(matches(e,EXT,'IgnoreCase',1));
        bg_img_fnames = fullfile(img_folder, {bg_img_fnames_tmp.name});
    end
    
    [bg_imgs, bg_im_ht, bg_im_wd, bg_ar] = load_images(bg_img_fnames);
    n_bg_imgs = numel(bg_imgs);
end

%% SETUP: match_to_sample mode
if m2s_mode
    
    n_rsvp = 2; % overwrite n_rsvp
    
    % load the task file
    if ~isempty(task_params_file)
        taskfile = dir(fullfile(task_folder, task_params_file));
    else
        assert(length(taskfile) == 1, sprintf('Number of possible taskfiles in %s does not equal 1', task_folder));
        taskfile = dir(fullfile(task_folder, '*task_params*.mat'));
    end
    
    tf_tmp = load(fullfile(taskfile.folder, taskfile.name));
    tf = tf_tmp.task;
    
    % load all sample_imgs
    %[sample_img_fnames, ~] = ret_image_fnames(tf.sample_img_folder)
    sample_img_fnames = fullfile(tf.sample_img_folder, tf.sample_img_names);
    [sample_imgs, sample_im_ht, sample_im_wd, sample_ar] = load_images(sample_img_fnames);
    
    % load all test imgs
    test_img_fnames = fullfile(tf.test_img_folder, tf.test_img_names);
    [test_imgs, test_im_ht, test_im_wd, test_ar] = load_images(test_img_fnames);
    n_test_imgs = numel(test_img_fnames);
    m2s_img_seq = zeros(n_rsvp, n_test_imgs);
    for i = 1:n_test_imgs
        m2s_img_seq(1,i) = tf.target_ID(i);
        m2s_img_seq(2,i) = i;
    end
    % convert the test grid (in fractional coordinates) to pix relative to im center
    test_img_x_cent = round(test_img_len/2);
    test_img_y_cent = round(test_img_ht/2);
    if iscell(tf.test_grid) % grid is specified in a per-test image basis
        test_pts = cell(numel(tf.test_grid),1);
        for tp = 1:numel(tf.test_grid)
            test_pts_tmp = tf.test_grid{tp} .* [test_img_len, test_img_ht];
            test_pts_tmp(:,1) = test_pts_tmp(:,1) - test_img_x_cent;
            test_pts_tmp(:,2) = test_pts_tmp(:,2) - test_img_y_cent;
            test_pts{tp} = test_pts_tmp;
        end
    else
        test_pts = tf.test_grid .* [test_img_len, test_img_ht];
        test_pts(:,1) = test_pts(:,1) - test_img_x_cent;
        test_pts(:,2) = test_pts(:,2) - test_img_y_cent;
    end
   % if apply_gaze_center_adj
    %    test_pts(:,2) = test_pts(:,2) + gaze_center_adj_y;
   %     test_pts(:,1) = test_pts(:,1) + gaze_center_adj_x;
    %end
end

%% SETUP: wake up images

if isfield(s, 'wake_up_trials')
    wake_up_trials = s.wake_up_trials;
    wake_up_every = s.wake_up_every;
    wake_up_tr_dur = s.wake_up_tr_dur;
    wake_up_img_fold = s.wake_up_img_fold;
    wake_up_stim_size_x = s.wake_up_stim_size_x;
    wake_up_stim_size_y = s.wake_up_stim_size_y;
    if isfield(s,'wake_up_reward')
        wake_up_reward = s.wake_up_reward;
    else
        wake_up_reward = 1;
    end
else
    wake_up_trials = false;
    wake_up_every = [];
    wake_up_tr_dur = [];
    wake_up_img_fold = '';
    wake_up_stim_size_x = [];
    wake_up_stim_size_y = [];
end

if wake_up_trials
    wu_image_fname_list = {};
    if iscell(wake_up_img_fold)
        wake_up_imgs = {};
        for j = 1:length(wake_up_img_fold)
            wu_img_d = dir(wake_up_img_fold{j});
            wu_img_d = wu_img_d(~[wu_img_d.isdir]);
            for i = 1:numel(wu_img_d)
                wu_image_fname_list = [wu_image_fname_list,fullfile(wake_up_img_fold{j}, wu_img_d(i).name)];
                [wu_img_tmp, ~,alpha] = imread(fullfile(wake_up_img_fold{j}, wu_img_d(i).name));
                if ~isempty(alpha)
                    wu_img_tmp(:,:,4) = alpha;
                end
                wake_up_imgs = [wake_up_imgs,wu_img_tmp];
            end
        end
        n_wu_imgs = numel(wake_up_imgs);
        %wu_img_idx = randi(n_wu_imgs);
        wu_idx = 1;
    else
        wu_img_d = dir(wake_up_img_fold);
        wu_img_d = wu_img_d(~[wu_img_d.isdir]);
        wake_up_imgs = cell(1,numel(wu_img_d));
        for i = 1:numel(wu_img_d)
            wu_image_fname_list{i} = fullfile(wake_up_img_fold, wu_img_d(i).name);
            [wu_img_tmp, ~,alpha] = imread(wu_image_fname_list{i});
            if ~isempty(alpha)
                wu_img_tmp(:,:,4) = alpha;
            end
            wake_up_imgs{i} = wu_img_tmp;
        end
        n_wu_imgs = numel(wake_up_imgs);
        %wu_img_idx = 1;
        wu_idx = 1;
    end
end

%% PTB open windows
PsychJavaSwingCleanup;
%InitializeMatlabOpenGL;
AssertOpenGL;
InitializePsychSound(1);

%Screen('Preference', 'Verbosity', 0);
Screen('Preference', 'Verbosity', 1);
%Screen('Preference', 'SkipSyncTests', skip_sync_tests)
Screen('Preference', 'SkipSyncTests', 0);
Screen('Preference', 'VisualDebugLevel', 0);
Screen('Preference', 'TextRenderer', 0);

% get colors
whitecol = WhiteIndex(screenid_stim);
blackcol = BlackIndex(screenid_stim);
whitecol_interrsvp = whitecol *0.8;
blackcol_interrsvp = whitecol * 0.2;
graycol = (whitecol + blackcol)/2;

%rect_col = [255 0 0];
rect_col = graycol;

switch bg_col
    case 'gray'
        bg_col_val = graycol;
    case 'white'
        bg_col_val = whitecol;
    case 'black'
        bg_col_val = blackcol;
end

%fprintf('GUI MODE\n\n')
%window_rect = [ 0   0 2000 864];
%screenid_stim = 1;
[win, win_rect] = Screen('OpenWindow', screenid_stim, bg_col_val, window_rect);
%ctrl_rect_debug = [0 0 1650 1200];
shrink_factor = 0.5;
ctrl_rect_debug = shrink_factor * win_rect;

screen_hz = Screen('NominalFrameRate', screenid_stim);
ifi=Screen('GetFlipInterval', win);
halfifi = 0.5*ifi;

if fix_exp_mode
    inter_rsvp_frames = round(rsvp_iti_t/1e3/ifi);
end

if ~isempty(lick_arduino)
    lick_rect = [ctrl_rect_debug(3)-90*shrink_factor, ctrl_rect_debug(4)-200, ctrl_rect_debug(3), ctrl_rect_debug(4)-200+90*shrink_factor];
end
if ctrl_screen; [win_ctrl, win_rect_ctrl] = Screen('OpenWindow', screenid_ctrl, bg_col_val, ctrl_rect_debug); end
%if ctrl_screen; [win_ctrl, win_rect_ctrl] = Screen('OpenWindow', screenid_ctrl, bg_col_val); end
%[win0, winRect0] = Screen('OpenWindow', screenId, bgcolor * 255, [0 0 300 300], [], [], [], [], [], kPsychGUIWindow);
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

HideCursor(screenid_stim);
HideCursor(screenid_ctrl);

if ctrl_screen
    Screen('BlendFunction', win_ctrl, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    ctrl_win_text_size = 18;
    Screen('TextSize', win_ctrl, ctrl_win_text_size);
end

[x_cent, y_cent] = RectCenter(win_rect);
if photodiode_flash
    flash_rect = [win_rect(3)-flash_rect_size(1), win_rect(4)-flash_rect_size(2),...
        win_rect(3), win_rect(4)];
end

if apply_gaze_center_adj
    y_cent = y_cent + gaze_center_adj_y;
    x_cent = x_cent + gaze_center_adj_x;
end

if m2s_mode
    % put test pts into screen space
    if iscell(test_pts)
        for tp = 1:numel(test_pts)
            test_pts{tp}(:,1) =  test_pts{tp}(:,1) + x_cent;
            test_pts{tp}(:,2) =  test_pts{tp}(:,2) + y_cent;
        end
    else
        test_pts(:,1) = test_pts(:,1) + x_cent;
        test_pts(:,2) = test_pts(:,2) + y_cent;
    end
end


%% TRIGGER: turn on session trigger
if trig_flag
    IOPort('Write', trig_hand, sess_trig_cmd.on, 1);
    WaitSecs(0.05);
    IOPort('Write', trig_hand, sess_trig_cmd.on, 1);
    
    % send file barcode to sample command line
    if ~isempty(sampleCommand_trig_cmd)
        time_str = replace(session_time(12:end),'-','');
        for i = 1:length(time_str)
            IOPort('Write',trig_hand,sampleCommand_trig_cmd.on,1);
            WaitSecs(0.01 * str2double(time_str(i)));
            IOPort('Write', trig_hand, sampleCommand_trig_cmd.off, 1);
            WaitSecs(0.025);
        end
    end
end

%% COLOR SHIFT FEEDBACK: DEPRECATED???
if color_shift_feedback
    n_frames_til_rew = round(time_to_reward/1e3/ifi);
    col_shift_change = nan(3,n_frames_til_rew);
    col_shift_change(1,:) = linspace(rew_col_start(1), rew_col_end(1), n_frames_til_rew);
    col_shift_change(2,:) = linspace(rew_col_start(2), rew_col_end(2), n_frames_til_rew);
    col_shift_change(3,:) = linspace(rew_col_start(3), rew_col_end(3), n_frames_til_rew);
    col_shift_change(4,:) = ones(1,n_frames_til_rew)*0.5;
end

%% LOAD: preload all textures
% load all image textures
if strcmp(stim_mode, 'images') || strcmp(stim_mode,'spinning') || strcmp(stim_mode, 'smooth pursuit')
    imgs_texture = zeros(1,length(imgs));
    imgs_texture_ctrl = zeros(1,length(imgs));
    for i = 1:length(imgs)
        imgs_texture(i) = Screen('MakeTexture', win,imgs{i});
        Screen('PreloadTextures', win,imgs_texture(i));
        if ctrl_screen
            imgs_texture_ctrl(i) = Screen('MakeTexture',win_ctrl,imgs{i});
            Screen('PreloadTextures',win_ctrl,imgs_texture_ctrl(i));
        end
    end
    
    
    if constant_bg_mode
        bg_imgs_texture = zeros(1,length(bg_imgs));
        bg_imgs_texture_ctrl = zeros(1,length(bg_imgs));
        for i = 1:length(bg_imgs)
            bg_imgs_texture(i) = Screen('MakeTexture', win,bg_imgs{i});
            Screen('PreloadTextures', win,bg_imgs_texture(i));
            if ctrl_screen
                bg_imgs_texture_ctrl(i) = Screen('MakeTexture',win_ctrl,bg_imgs{i});
                Screen('PreloadTextures',win_ctrl,bg_imgs_texture_ctrl(i));
            end
        end
    end
elseif m2s_mode
    % load SAMPLE images 
    sample_imgs_texture = zeros(1,length(sample_imgs));
    sample_imgs_texture_ctrl = zeros(1,length(sample_imgs));
    for i = 1:length(sample_imgs)
        sample_imgs_texture(i) = Screen('MakeTexture', win,sample_imgs{i});
        Screen('PreloadTextures', win,sample_imgs_texture(i));
        if ctrl_screen
            sample_imgs_texture_ctrl(i) = Screen('MakeTexture',win_ctrl,sample_imgs{i});
            Screen('PreloadTextures',win_ctrl,sample_imgs_texture_ctrl(i));
        end
    end
    
    % load TEST images
    test_imgs_texture = zeros(1,length(test_imgs));
    test_imgs_texture_ctrl = zeros(1,length(test_imgs));
    for i = 1:length(test_imgs)
        test_imgs_texture(i) = Screen('MakeTexture', win,test_imgs{i});
        Screen('PreloadTextures', win,test_imgs_texture(i));
        if ctrl_screen
            test_imgs_texture_ctrl(i) = Screen('MakeTexture',win_ctrl,test_imgs{i});
            Screen('PreloadTextures',win_ctrl,test_imgs_texture_ctrl(i));
        end
    end
    
end


%% LOAD: load all wake up images
if wake_up_trials
    %wu_imgs_texture = cell(1,length(wake_up_imgs));
    %wu_imgs_texture_ctrl = cell(1,length(wake_up_imgs));
    wu_imgs_texture = zeros(1,length(wake_up_imgs));
    wu_imgs_texture_ctrl = zeros(1,length(wake_up_imgs));
    for i = 1:length(wake_up_imgs)
        %wu_imgs_texture{i} = Screen('MakeTexture',win,wake_up_imgs{i});
        wu_imgs_texture(i) = Screen('MakeTexture',win,wake_up_imgs{i});
        if ctrl_screen
            %wu_imgs_texture_ctrl{i} = Screen('MakeTexture',win_ctrl,wake_up_imgs{i});
            wu_imgs_texture_ctrl(i) = Screen('MakeTexture',win_ctrl,wake_up_imgs{i});
        end
    end
end

% wake_up movies (keyboard v)
if isfield(s,'wake_up_movie')
    moviefiles_all = dir(s.wake_up_movie);
    movienames =  {moviefiles_all(~[moviefiles_all.isdir]).name};
    moviePtr_wu = Screen('OpenMovie', win, fullfile(s.wake_up_movie, movienames{1}), [], 1, 64);
    Screen('PlayMovie', moviePtr_wu, movie_rate, 1);
    
    wu_movie_stim_rect = CenterRectOnPoint([0, 0, wake_up_stim_size_x, wake_up_stim_size_y], x_cent,y_cent);
end

if strcmp(stim_mode, 'movie')
    moviefiles_all = dir(movie_folder);
    movienames =  {moviefiles_all(~[moviefiles_all.isdir]).name};
    moviePtr = Screen('OpenMovie', win, fullfile(movie_folder, movienames{1}), [], 1, 64);
    Screen('PlayMovie', moviePtr, movie_rate, 1);
    WaitSecs(1);
end

%% SETUP: set up stimulus locations on screen

calib_rect = [0 0 calibration_win_len, calibration_win_ht];
calib_rect_cent = CenterRectOnPoint(calib_rect, x_cent, y_cent);

% calculate x locs, y locs for stim
if n_pts_x == 1
    x_pts = round(mean([calib_rect_cent(1), calib_rect_cent(3)]));
else
    x_pts = round(linspace(calib_rect_cent(1), calib_rect_cent(3), n_pts_x));
end

if n_pts_y == 1
    y_pts = round(mean([calib_rect_cent(2), calib_rect_cent(4)]));
else
    y_pts = round(linspace(calib_rect_cent(2), calib_rect_cent(4), n_pts_y));
end

[A,B] = meshgrid(x_pts,y_pts);
c=cat(2,A',B');
all_pts =reshape(c,[],2);
n_calib_pts = size(all_pts, 1);

% if add center and gridpoints don't always contain center points
if isfield(s,'add_center')
    if s.add_center == 1
        if (find(all_pts(:,1) ==x_cent) ~= find(all_pts(:,2) == y_cent))
            n_calib_pts = n_calib_pts + 1;
            n_pts_x = n_pts_x + 1;
            n_pts_y = n_pts_y + 1;
            all_pts = [all_pts; x_cent,y_cent];
        elseif (isempty(find(all_pts(:,1) == x_cent, 1)) && isempty(find(all_pts(:,2) == y_cent, 1)))
            n_calib_pts = n_calib_pts + 1;
            n_pts_x = n_pts_x + 1;
            n_pts_y = n_pts_y + 1;
            all_pts = [all_pts; x_cent,y_cent];
        end
    end
end

%% SETUP: spinning: DEPRECATED???
if strcmp(stim_mode,'spinning')
    center_coord = [win_rect(3)/2, win_rect(4)/2];
    
    per_circle = diam_circle * pi; % circle perimeter
    
    % calculate how long it takes for a full turn
    
    presentation_time = per_circle/speed; % second
    presentation_time = presentation_time * 1000; % millisecond
    
    % calculate frames per stim, iti
    stim_frames = round(presentation_time/1e3/ifi);
    n_pts_x = stim_frames;
    n_pts_y = stim_frames;
    
    % pick stim_frames number of points on this trajectory
    all_angle = linspace(0,pi*2,stim_frames);
    
    x_pts = zeros(1,length(all_angle));
    y_pts = zeros(1,length(all_angle));
    for i = 1:length(x_pts)
        x_pts(i) = diam_circle/2 * cos(all_angle(i)) + center_coord(1);
        y_pts(i) = diam_circle/2 * sin(all_angle(i)) + center_coord(2);
    end
    
    c=cat(2,x_pts',y_pts');
    all_pts =reshape(c,[],2);
    n_calib_pts = size(all_pts, 1);
end

if strcmp(stim_mode,'smooth pursuit')
    center_coord = [win_rect(3)/2, win_rect(4)/2];
    % calculate frames per stim, iti
    stim_frames = round(presentation_time/1e3/ifi);
    n_pts_x = stim_frames;
    n_pts_y = stim_frames;
    
    all_angle = zeros(size(angle_endpts));
    for q = 1:length(angle_endpts)
        all_angle(q) = angle_endpts(q) * pi/180;
    end
    
    x_pts = zeros(1,length(all_angle));
    y_pts = zeros(1,length(all_angle));
    for i = 1:length(x_pts)
        x_pts(i) = diam_circle/2 * cos(all_angle(i)) + center_coord(1);
        y_pts(i) = diam_circle/2 * sin(all_angle(i)) + center_coord(2);
    end
    c=cat(2,x_pts',y_pts');
    
    pt_start = repmat(center_coord,[size(c,1),1]);
    pt_end = c;
    
    all_unique_pts = cell(size(c,1),1);
    for q = 1:size(c,1)
        all_unique_pts{q} = zeros(stim_frames,2);
        all_unique_pts{q}(:,1) = linspace(pt_start(q,1),pt_end(q,1), stim_frames);
        all_unique_pts{q}(:,2) = linspace(pt_start(q,2),pt_end(q,2),stim_frames);
    end
    
    n_calib_pts = stim_frames;
    
end

%% SETUP: img sequence
if strcmp(stim_mode, 'images')
    n_trs_requested = n_calib_pts*trs_per_location;
    if constant_bg_mode
        n_reps = ceil(n_trs_requested*n_rsvp/sum(n_imgs_per_folder));
        img_seq = [];
        bg_seq = [];
        for i = 1:length(unique_img_folder)
            stidx = sum(n_imgs_per_folder(1:i-1))+1;
            img_seq_tmp = repmat(stidx:stidx+sum(n_imgs_per_folder(i))-1, [1, n_reps]);
            n_trs_per_folder(i) = floor(length(img_seq_tmp)/n_rsvp);
            n_stim_per_folder(i) = n_trs_per_folder(i) * n_rsvp;
            img_seq_tmp = img_seq_tmp(1:n_stim_per_folder(i));
            ridx_tmp = randperm(length(img_seq_tmp));
            img_seq_rnd = img_seq_tmp(ridx_tmp);
            img_seq = [img_seq; reshape(img_seq_rnd, [numel(img_seq_rnd)/n_rsvp, n_rsvp])];
            bg_seq = [bg_seq; ones(n_trs_per_folder(i), 1)*i];
        end
        
        % trial randomization
        tr_ridx = randperm(numel(bg_seq));
        bg_seq = bg_seq(tr_ridx, :);
        img_seq = img_seq(tr_ridx, :);
        img_seq = img_seq';
        
    elseif strcmp(image_order, 'random')
        n_rsvps = n_trs_requested*n_rsvp;
        x = floor(n_rsvps/n_imgs);
        if x>0
            img_seq = [];
            for q = 1:x
                img_seq = [img_seq randperm(n_imgs)];
            end
            img_seq = [img_seq randperm(n_imgs, mod(n_rsvps, n_imgs))];
        else
            img_seq = randperm(n_imgs, n_rsvps);
        end
        
    elseif strcmp(image_order,'block')
        % randomly draw from each image folder, but present folders
        % sequentially as listed in the img_folder
        % number of trials specified in the Marmulator GUI will now apply
        % to each image folder
        img_seq = [];
        for i = 1:length(unique_img_folder)
            n_rsvps = n_trs_requested/length(unique_img_folder)*n_rsvp;
            x = floor(n_rsvps/n_imgs_per_folder(i));
            if i > 1
                idx_start = n_imgs_per_folder(i-1);
            else
                idx_start = 0;
            end
            if x>0
                for q = 1:x
                    img_seq = [img_seq randperm(n_imgs_per_folder(i))+idx_start];
                end
                img_seq = [img_seq randperm(n_imgs_per_folder(i),mod(n_rsvps,n_imgs_per_folder(i)))+idx_start];
            else
                img_seq = [img_seq randperm(n_imgs_per_folder(i), n_rsvps)+idx_start];
            end
        end
        
    else % sequential presentation
        n_rsvps = n_trs_requested*n_rsvp;
        x = floor(n_rsvps/n_imgs);
        img_seq = [repmat(1:n_imgs, 1, x) 1:mod(n_rsvps, x*n_imgs)];
    end
    
    if ~constant_bg_mode
        img_seq = reshape(img_seq, [n_rsvp, n_trs_requested]);
    end
    
    n_trs_curr = size(img_seq, 2);
elseif m2s_mode
    n_trs_requested = trs_per_location;
    x = floor(n_trs_requested/n_test_imgs);
    if x>0
        img_seq = [];
        for q = 1:x
            img_seq = [img_seq randperm(n_test_imgs)];
        end
        img_seq = [img_seq randperm(n_test_imgs, mod(n_trs_requested, n_test_imgs))];
    else
        img_seq = randperm(n_test_imgs, n_trs_requested);
    end
    img_seq = m2s_img_seq(:,img_seq);
    n_trs_curr = n_trs_requested;

    target_img_seq = img_seq(2,:); 
    target_idx = tf.target_location(target_img_seq); % which location is the correct one? 
else
    n_trs_requested = n_calib_pts*trs_per_location;
    img_seq = [];
    image_displayed = [];
    n_trs_curr = n_trs_requested;
end

%% SETUP: trial sequence
if strcmp(image_order,'block')
    trs_per_location = trs_per_location * length(unique_img_folder);
end

tr_seq = [];
if strcmp(position_order, 'random')
    for q = 1:ceil(n_trs_curr/n_calib_pts)
        tr_seq = [tr_seq randperm(n_calib_pts)];
    end
else
    tr_seq = repmat(1:n_calib_pts,1,ceil(n_trs_curr/n_calib_pts));
end

tr_seq = tr_seq(1:n_trs_curr);

%% SETUP: wake_up images
if wake_up_trials && ~m2s_mode % for now, m2s doesn't work with wakeup images
    %y = nan(1,n_trs_requested);
    y = nan(1,n_trs_curr);
    z = [tr_seq; y];
    %n_wake_up_trs = ceil(n_trs_requested/mean(wake_up_every))+1;
    n_wake_up_trs = ceil(n_trs_curr/mean(wake_up_every))+1;
    wue = wake_up_every -1;
    insert_every = round(rand(1,n_wake_up_trs)*diff(wue)+wue(1));
    insert_idx = cumsum(insert_every);
    %insert_idx = insert_idx(insert_idx<=n_trs_requested);
    insert_idx = insert_idx(insert_idx<=n_trs_curr);
    z(2,insert_idx) = 0;
    z = z(~isnan(z))';
    tr_seq = z;
    n_trs_tot = numel(tr_seq);
    n_wake_up_trs = sum(tr_seq == 0);
    wake_up_stim_frames = round(wake_up_tr_dur/ifi);
    
    % wake up image sequence randomized
    wu_seq = repmat(linspace(1,n_wu_imgs,n_wu_imgs),[1,ceil(n_wake_up_trs/n_wu_imgs)]);
    wu_seq = wu_seq(randperm(length(wu_seq)));
    
    % modify img sequence as well
    if strcmp(stim_mode, 'images')
        wu_tr_idx = find(tr_seq == 0);
        for i = 1:numel(wu_tr_idx)
            wu_tr_idx_curr = wu_tr_idx(i);
            img_seq = [img_seq(:,1:wu_tr_idx_curr-1) zeros(n_rsvp, 1) img_seq(:,wu_tr_idx_curr:end)];
            if constant_bg_mode
                bg_seq = [bg_seq(1:wu_tr_idx_curr-1); 0; bg_seq(wu_tr_idx_curr:end)];
            end
        end
    end
    
else
    n_wake_up_trs = 0;
    n_trs_tot = n_trs_requested;
end

image_displayed = cell(n_rsvp, n_trs_tot); % will be empty for all except 'images' sessions

if strcmp(stim_mode,'spinning') || strcmp(stim_mode,'smooth pursuit')
    n_trs_tot = trs_per_location;
end

% desired trial clip sequence
if stimulus_pre_dot
    clip_sequence = {'fixation pt','stimulus'};
    if require_fix_tr_init == 1
        clip_sequence_t = [0,fixation_to_init];
    else
        clip_sequence_t = [0, stimulus_pre_time];
    end
    
else
    clip_sequence = {'stimulus'};
    clip_sequence_t = [0];
end


if fix_exp_mode
    for i = 1:n_rsvp
        if i == n_rsvp
            clip_sequence_t = [clip_sequence_t,clip_sequence_t(2) + rsvp_iti_t*(i-1) + presentation_time*i];
            clip_sequence = [clip_sequence,'trial_end'];
        else
            clip_sequence_t = [clip_sequence_t,clip_sequence_t(2) + rsvp_iti_t*(i-1) + presentation_time*i,...
                clip_sequence_t(2) + rsvp_iti_t*i + presentation_time*i];
            clip_sequence = [clip_sequence,'blank','stimulus'];
        end
    end
else
    clip_sequence_t = [clip_sequence_t,presentation_time];
    clip_sequence = [clip_sequence,'trial_end'];
end

if lick_flag
    lick_trial = false(1,n_trs_tot);
else
    lick_trial = [];
end

% calculate frames per stim, iti
stim_frames = round(presentation_time/1e3/ifi);
pulse_rate = (rand(1, n_trs_tot)+0.5)*3;

if strcmp(trial_mode, 'trial')
    t_sin = 0:ifi:(presentation_time/1e3);
elseif strcmp(trial_mode, 'foraging')
    t_sin = 0:ifi:(time_out_after+time_to_reward/1e3);
end

%% SETUP: bounding boxes

% set up bounding boxes for display on the control screen
cols = round(distinguishable_colors(n_calib_pts+1)*255);

bounding_rect_tmp = [0,0,bounding_rect_size_x, bounding_rect_size_y];

if ~isempty(manual_bounding_boxes)
    if apply_gaze_center_adj
        manual_bounding_boxes = OffsetRect(manual_bounding_boxes, gaze_center_adj_x, gaze_center_adj_y);
    end
    bounding_rects = manual_bounding_boxes;
else
    bounding_rects = CenterRectOnPoint(bounding_rect_tmp, all_pts(:,1), all_pts(:,2)); % all_pts is already gaze center adjusted
end


bounding_rect_stim_pre_dot_tmp = [0,0,bounding_rect_size_stim_pre_dot_x, bounding_rect_size_stim_pre_dot_y];
if ~isempty(manual_bounding_boxes_stim_pre_dot)
    if apply_gaze_center_adj
        manual_bounding_boxes_stim_pre_dot = OffsetRect(manual_bounding_boxes_stim_pre_dot,gaze_center_adj_x, gaze_center_adj_y);
    end
    bounding_rects_stim_pre_dot = manual_bounding_boxes_stim_pre_dot;
else
    bounding_rects_stim_pre_dot = CenterRectOnPoint(bounding_rect_stim_pre_dot_tmp,all_pts(:,1),all_pts(:,2));
end

if m2s_mode
    % define test bounding boxes
    test_bounding_rect_tmp = [0,0,test_bb_size_x, test_bb_size_y];
    if iscell(test_pts) % test points are defined per-test image
        test_bounding_rects = cell(numel(test_pts), 1);
        test_rects_idx = cell(numel(test_pts), 1);
        for tp = 1:numel(test_pts)
            test_bounding_rects{tp} = CenterRectOnPoint(test_bounding_rect_tmp, test_pts{tp}(:,1), test_pts{tp}(:,2));
            test_rects_idx{tp} = 1:size(test_pts{tp},1);
        end
    else
        test_bounding_rects = CenterRectOnPoint(test_bounding_rect_tmp, test_pts(:,1), test_pts(:,2));
        test_rects_idx = 1:size(test_bounding_rects, 1);
    end
end


%% SETUP: Reward and Punish sound handles using PsychPortAudio
audio_handle = [];
% reward
[aud_y, aud_fs] = psychwavread(reward_sound_file);
[samplecount,ninchannels] = size(aud_y);
aud_y = repmat(aud_y',2/ninchannels,1);
suggestedLat = []; % PsychPortAudio('GetDevices') LowOutputLatency
audio_handle(1) = PsychPortAudio('Open', 1, [], 1, aud_fs,2, [],suggestedLat);
PsychPortAudio('FillBuffer', audio_handle(1), aud_y);

% punish
[aud_pun_y, aud_pun_fs] = psychwavread(punish_sound_file);
[samplecount,ninchannels] = size(aud_pun_y);
aud_pun_y = repmat(aud_pun_y',2/ninchannels,1);
audio_handle(2) = PsychPortAudio('Open', 1, [], 1, aud_pun_fs,2,[],suggestedLat);
PsychPortAudio('FillBuffer', audio_handle(2), aud_pun_y);

%% SETUP: ITI, timing 

if iti_random
    iti_frames = round((rand(1,n_trs_tot)*diff(inter_stim_interval)+inter_stim_interval(1))/1e3/ifi);
else
    iti_frames = ones(1,n_trs_tot)*round(inter_stim_interval(1)/1e3/ifi);
end

stimulus_pre_frames = round(stimulus_pre_time/1e3/ifi); % strimulus_pre_time in ms
wait_after_rew_frames = round(wait_after_reward/ifi); % wait_after_reward in s

if m2s_mode 
    foraging_max_t = time_out_after + presentation_time + rsvp_iti_t; 
else
    foraging_max_t = time_out_after; 
end

if strcmp(stim_mode,'spinning')
    tr_seq = nan(n_calib_pts,n_trs_tot);
    
    for q = 1:n_trs_tot
        tr_seq(1,q) = randi(n_calib_pts);
    end
    if strcmp(position_order, 'random')
        for q = 1: n_trs_tot
            if rand(1) > 0.5 %ccw
                tr_seq(2:length(tr_seq),q)...
                    = cat(2,linspace(tr_seq(1,q)+1,n_calib_pts,abs(n_calib_pts-tr_seq(1,q)-1+1)),...
                    linspace(1,tr_seq(1,q)-1,abs(tr_seq(1,q)-1-1+1)));
            else
                tr_seq(2:length(tr_seq),q)...
                    = cat(2,linspace(tr_seq(1,q)-1,1,abs(1-tr_seq(1,q)-1+1)),...
                    linspace(n_calib_pts,tr_seq(1,q)+1,abs(tr_seq(1,q)-n_calib_pts)));
            end
        end
    else % always cw
        for q = 1:n_trs_tot
            tr_seq(2:length(tr_seq),q)...
                = cat(2,linspace(tr_seq(1,q)+1,n_calib_pts,abs(n_calib_pts-tr_seq(1,q)-1+1)),...
                linspace(1,tr_seq(1,q)-1,abs(tr_seq(1,q)-1-1+1)));
        end
    end
end

if strcmp(stim_mode,'smooth_pursuit')
    tr_seq = randi(length(all_pts),[n_trs_tot,1]);
end

%% SETUP: initialize state vars and empty vectors for run

reward_trial = false(1,n_trs_tot); % for recording whether reward was administered
correct_trial = false(1,n_trs_tot); % for recording whether a trial was correct - not necesssarily same as reward_trial if using stop_rewards_n_nolicks
punish_trial = false(1,n_trs_tot); % whether a trial was punished
reward_time = nan(1,n_trs_tot); % for recording the time of a reward
time_to_bb = nan(1,n_trs_tot);
trial_aborted = false(1, n_trs_tot);
manual_reward_t = nan(1,n_trs_tot);
stim_pre_start = nan(1,n_trs_tot);
stim_pre_end = nan(1,n_trs_tot);
trial_init_timed_out = false(1,n_trs_tot);
calib_st_t = nan(1,n_trs_tot);
calib_end_t = nan(1,n_trs_tot);
calib_t_clip = nan * ones(n_trs_tot,length(clip_sequence_t));
calib_frame_st_t = cell(1,n_trs_tot);
calib_frame_clip_num = cell(1,n_trs_tot);
wake_up_image_displayed = cell(1, n_trs_tot);
incorrect_choice_m2s = zeros(1, n_trs_tot); 
if n_rsvp>1
    rsvp_start_t = nan(n_rsvp, n_trs_tot);
    rsvp_end_t = nan(n_rsvp, n_trs_tot);
end
reward_ct = 0;
man_reward_ct = 0;
bonus_reward_ct = 0;
bonus_reward_t = [];
wu_mov_start_t = [];
wu_mov_end_t = [];
exit_flag = 0;
mov_dur = [];
fix_pre_fr_ctr = 0;

if profile_memory
    mem_used = nan(1,n_trs_tot);
    avail_sys_mem = nan(1,n_trs_tot);
    avail_phys_mem = nan(1,n_trs_tot);
end

%% RUN: start the stimulus loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Screen('FillRect', win, bg_col_val)
seqidx = 0;
vbl = Screen('Flip', win);
t_start_sec = GetSecs();

for i = 1:n_trs_tot
    
    % set variables for each trial
    curr_tr = i;
    seqidx = tr_seq(i);
    stfridx = 0;
    rsvpfridx = 0;
    rsvp_ctr = 1; 
    interrsvpfridx = 0;
    loop_brk_ctr = 0;
    punish_loop_brk_ctr = 0; 
    pre_stim_timer = 0;
    manual_reward_flag = false; % has a manual reward been given on the trial yet?
    entered_bb = false; % did eye go into bb?
    skip_reward = 0;
    frame_st_t = [];
    frame_clip_num = [];
    t_dur = [];
    
    if seqidx ~= 0
        xcurr = all_pts(seqidx,1);
        ycurr = all_pts(seqidx,2);
        stim_rect = [xcurr-stim_rect_size_x/2 ycurr-stim_rect_size_y/2 xcurr+stim_rect_size_x/2 ycurr+stim_rect_size_y/2];
        if m2s_mode % we need sample and test bounding boxes
            sample_bb = bounding_rects(seqidx,:);
            if iscell(test_bounding_rects)
                test_bounding_rects_curr = test_bounding_rects{img_seq(end,curr_tr)};
                test_rects_idx_curr = test_rects_idx{img_seq(end,curr_tr)}; 
            else
                test_bounding_rects_curr = test_bounding_rects;
                test_rects_idx_curr = test_rects_idx; 
            end

            test_bb = test_bounding_rects_curr(target_idx(curr_tr),:);
            punish_rects = test_rects_idx_curr(test_rects_idx_curr ~= target_idx(curr_tr));
            punish_bbs = test_bounding_rects_curr(punish_rects,:);
            n_punish_bbs = size(punish_bbs, 1);
            
            curr_bb = sample_bb;
        else
            curr_bb = bounding_rects(seqidx,:);
        end
        
        curr_bb_stim_pre_dot = bounding_rects_stim_pre_dot(seqidx,:);
    else
        % put stim rect around center
        xcurr = x_cent;
        ycurr = y_cent;
        stim_rect = CenterRectOnPoint([0, 0, wake_up_stim_size_x, wake_up_stim_size_y], xcurr, ycurr);
        curr_bb = stim_rect;
        curr_bb_stim_pre_dot = stim_rect;
        
        if wake_up_reward == 0
            skip_reward = 1;
            reward_this_trial = 0;
        end
    end
    
    idx_all = idx_all +1;
    
    drawInfoText();
    drawBoundingBoxes();
    if draw_retain_bb_pts
        drawGoodEyePts(0);
    end
    get_eyetracker_draw_dots();
    checkLick();
    
    Screen('FillRect', win, bg_col_val)
    vbl = Screen('Flip', win);
    
    if ctrl_screen
        Screen('FillRect', win_ctrl, bg_col_val)
        vbl2 = Screen('Flip', win_ctrl);
    end
    
    for j = 1:iti_frames(i)
        idx_all = idx_all +1;
        drawInfoText();
        
        drawBoundingBoxes();
        if draw_retain_bb_pts
            drawGoodEyePts(0);
        end
        checkLick();
        
        vbl = Screen('Flip', win, vbl + halfifi, dontclear, dontsync, multiflip);
        if ctrl_screen; vbl2 = Screen('Flip', win_ctrl, vbl2 + halfifi, dontclear2, dontsync2, multiflip2); end
        get_eyetracker_draw_dots();
    end
    
    end_stim_pre = 0;
    j = 0;
    fix_pre_fr_ctr = 0;
    blink_ctr = 0;
    trial_trig_hi = 0;
    sampleCommand_trig_hi = 0;
    clip_ctr = 1;
    
    % ENTER STIMULUS PRE
    if stimulus_pre_dot
        stps = GetSecs();
        stim_pre_start(i) = stps - t_start_sec;
        while ~end_stim_pre
            
            j = j+1;
            idx_all = idx_all +1;
            
            % draw background if in constant background mode
            if seqidx~= 0 && constant_bg_mode
                bgidx = bg_seq(i);
                rsx_curr = stim_rect_size_x;
                rsy_curr = round(rsx_curr/bg_ar(bgidx));
                bg_stim_rect_curr = [xcurr-rsx_curr/2 ycurr-rsy_curr/2 xcurr+rsx_curr/2 ycurr+rsy_curr/2];
                
                Screen('DrawTexture', win, bg_imgs_texture(bgidx), [], bg_stim_rect_curr)
                if ctrl_screen
                    Screen('DrawTexture', win_ctrl, bg_imgs_texture_ctrl(bgidx), [], bg_stim_rect_curr*shrink_factor)
                end
            end
            
            drawInfoText();
            checkLick();
            
            if seqidx ~= 0
                if draw_crosshairs
                    drawCrosshairs(all_pts(seqidx,:));
                else
                    if ctrl_screen; Screen('DrawDots', win_ctrl, all_pts(seqidx,:)*shrink_factor, stim_pre_dot_sz*shrink_factor, whitecol, [], 1); end
                    Screen('DrawDots', win, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                end
            else
                if draw_crosshairs
                    drawCrosshairs([xcurr, ycurr]); 
                else
                    if ctrl_screen; Screen('DrawDots', win_ctrl, [xcurr, ycurr]*shrink_factor, stim_pre_dot_sz*shrink_factor, whitecol, [], 1); end
                    Screen('DrawDots', win, [xcurr, ycurr], stim_pre_dot_sz, whitecol, [], 1);
                end
            end
            
            [eyeposx_cur, eyeposy_cur] = get_eyetracker_draw_dots();
            drawBoundingBoxes(curr_bb_stim_pre_dot);
            
            vbl = Screen('Flip', win, vbl + halfifi, dontclear, dontsync, multiflip);
            if ctrl_screen; vbl2 = Screen('Flip', win_ctrl, vbl2 + halfifi, dontclear2, dontsync2, multiflip2); end
            
            if trig_flag && ~trial_trig_hi
                IOPort('Write', trig_hand, trial_trig_cmd.on, 1);
                trial_trig_hi = 1;
            end
            
            %[eyeposx_cur, eyeposy_cur] = get_eyetracker_draw_dots();
            curr_in_bb = IsInRect(eyeposx_cur, eyeposy_cur, curr_bb_stim_pre_dot);
            
            pre_stim_timer = GetSecs() - stps;
            [~,~,kCode] = KbCheck(0);
            if find(kCode) == esc_key
                disp('ESC key recognized, exiting');
                exit_flag = 1;
                trial_aborted(i) = true;
                break
            elseif find(kCode) == bonus_rew_key
                if isempty(bonus_reward_t) || GetSecs() - bonus_reward_t(end) - t_start_sec > bonus_reward_min_wait
                    bonus_reward_t(end+1) = GetSecs() - t_start_sec;
                    %man_reward_ct = man_reward_ct + 1;
                    bonus_reward_ct = bonus_reward_ct + 1;
                    pumpReward_updateGUI();
                    fprintf('bonus reward, time = %0.3fs\n', bonus_reward_t(end));
                end
            elseif any(find(kCode) == wake_up_movie_key) && isfield(s,'wake_up_movie')
                wake_up_movie = 1;
                wu_mov_start_t(end+1) = GetSecs()-t_start_sec;
                while wake_up_movie
                    disp('wake_up_movie is playing')
                    drawInfoText();
                    checkLick();
                    Screen('PlayMovie', moviePtr_wu, movie_rate, 1);
                    [eyeposx_cur, eyeposy_cur] = get_eyetracker_draw_dots();
                    [keyIsDown, ~, keyCode] = KbCheck(-1);
                    if (keyIsDown==1 && keyCode(wake_up_movie_terminate_key))
                        Screen('PlayMovie',moviePtr_wu,0);
                        wu_mov_end_t(end+1) = GetSecs()-t_start_sec;
                        break;
                    end
                    
                    waitforimage = 1;
                    movtex_new_wu = Screen('GetMovieImage', win, moviePtr_wu, waitforimage);
                    %Screen('PlayMovie', movie, rate, 1, 1.0);
                    if movtex_new_wu>0
                        movtex_wu = movtex_new_wu;
                    end
                    
                    if ctrl_screen
                        Screen('DrawTexture', win_ctrl, movtex_wu, [],  wu_movie_stim_rect*shrink_factor);
                        Screen('DrawDots', win_ctrl, [eyeposx_cur, eyeposy_cur]*shrink_factor,...
                            dotsz, rgba_cols(:,end), [], 1);
                    end
                    Screen('DrawTexture', win, movtex_wu, [],  wu_movie_stim_rect)
                    vbl = Screen('flip',win,vbl+halfifi,dontclear,dontsync,multiflip);
                    if ctrl_screen; vbl2 = Screen('Flip', win_ctrl,vbl2+halfifi, dontclear2,dontsync2,multiflip2);end
                end
                wake_up_movie = 0;
            end
            
            if require_fix_tr_init && seqidx ~= 0
                if location_require_quality
                    qual_check = eyetracker_qual(idx_all) < 2;
                else
                    qual_check = true;
                end
                
                if curr_in_bb  && qual_check
                    fix_pre_fr_ctr = fix_pre_fr_ctr + 1;
                    blink_ctr = 0;
                    
                    if fix_pre_fr_ctr == 3
                        calib_t_clip(i,clip_ctr) = GetSecs() - t_start_sec;
                        clip_ctr = clip_ctr + 1;
                    end
                elseif (~curr_in_bb || ~qual_check) && fix_pre_fr_ctr>3 && blink_ctr<n_frames_blink
                    blink_ctr = blink_ctr + 1;
                else
                    blink_ctr = 0;
                    fix_pre_fr_ctr = 0;
                end
                
                if fix_pre_fr_ctr >= round(fixation_to_init/1e3/ifi)
                    end_stim_pre = 1;
                elseif pre_stim_timer>time_out_trial_init_s
                    end_stim_pre = 1;
                    trial_init_timed_out(i) = 1;
                end
                
            elseif j >=stimulus_pre_frames
                end_stim_pre = 1;
            end
        end
        stim_pre_end(i) = GetSecs() - t_start_sec;
    end
    
    if strcmp(stim_mode, 'moving_dot')
        curr_rect = stim_rect;
        curr_x = round(rand * [curr_rect(3) - curr_rect(1)]) + curr_rect(1);
        curr_y = round(rand * [curr_rect(4) - curr_rect(2)]) + curr_rect(2);
        dx_sign = sign(randn);
        dy_sign = sign(randn);
        size_curr = mvdot_sz;
        dsize_sign = 1;
    end
    
    eye_data_qual_curr = nan(1,stim_frames);
    eye_in_bb_curr = nan(1,stim_frames);
    
    if trial_init_timed_out(i)
        end_stim = 1;
    else
        end_stim = 0;
    end
    blink_ctr = 0;
    %rsvp_ctr = 1;
    rsvp_break_ctr = 0;
    inter_rsvp_fr_ctr = 1;
    pulse_t_idx = 0; % for size-pulsed imaged
    stim_trig_hi = 0;
    
    % ENTER STIMULUS
    while ~end_stim && ~exit_flag % loops for every frame
        t1_frame = GetSecs();
        if fix_exp_mode && rsvp_ctr > 1 && inter_rsvp_fr_ctr <= inter_rsvp_frames && seqidx ~= 0
            % go into a break
            inter_rsvp = true;
            interrsvpfridx = interrsvpfridx + 1;
        else
            inter_rsvp = false;
            rsvpfridx = rsvpfridx + 1;
            interrsvpfridx = 0;
            %rsvpfridx
        end
        
        if (strcmp(stim_mode, 'images') || m2s_mode)&& seqidx ~=0
            img_idx = img_seq(rsvp_ctr, i);
        end
        
        if m2s_mode
            if rsvp_ctr == 1 
                curr_bb = sample_bb; 
            else
                curr_bb = test_bb; 
            end
        end
        
        stfridx = stfridx + 1; % stim frame idx, counts every frame from stim start (rsvp_1) through trial stim end (rsvp_n)
        idx_all = idx_all +1; % idx for all frames recorded
        
        if ~trial_init_timed_out(i)
            % draw gray background:
            Screen('FillRect', win, bg_col_val);
            if ctrl_screen; Screen('FillRect', win_ctrl, bg_col_val); end
            
            if seqidx ~= 0 && strcmp(trial_mode, 'foraging') && color_shift_feedback && loop_brk_ctr>0
                Screen('FillRect', win, col_shift_change(:,loop_brk_ctr), curr_bb);
                if ctrl_screen; Screen('FillRect', win_ctrl, col_shift_change(:,loop_brk_ctr), curr_bb); end
            end
            
            [eyeposx_cur, eyeposy_cur] = get_eyetracker_draw_dots();
            eye_data_qual_curr(stfridx) = eyetracker_qual(idx_all);

            curr_in_bb = IsInRect(eyeposx_cur, eyeposy_cur, curr_bb);

            %t_pb1 = GetSecs(); 
            if m2s_mode && rsvp_ctr>1 && ~inter_rsvp
                in_punish_bb_tmp = false(1,n_punish_bbs); 
                for p = 1:n_punish_bbs
                    in_punish_bb_tmp(p) = IsInRect(eyeposx_cur,eyeposy_cur,punish_bbs(p,:)); 
                end
                if any(in_punish_bb_tmp) 
                    curr_in_bb_punish = true; 
                    in_punish_bb = find(in_punish_bb_tmp); 
                else
                    curr_in_bb_punish = false; 
                end
                %t_pb2 = GetSecs();            
            end
            
            %if m2s_mode && rsvp_ctr>1 && curr_in_bb
                %disp('curr in bb')
                %fprintf('curr in bb: %d\n', loop_brk_ctr); 
                %loop_brk_ctr
            %end
            
            eye_in_bb_curr(stfridx) = curr_in_bb;
            
            if draw_retain_bb_pts
                drawGoodEyePts(1);
            end
            
            if seqidx ~= 0  && ~inter_rsvp
                switch stim_mode
                    case 'spinning'
                        seqidx = tr_seq(stfridx,i);
                        xcurr = all_pts(seqidx,1);
                        ycurr = all_pts(seqidx,2);
                        stim_rect_curr = [xcurr-stim_rect_size_x/2 ycurr-stim_rect_size_y/2 xcurr+stim_rect_size_x/2 ycurr+stim_rect_size_y/2];
                        %Screen('DrawTexture', win, imgs_texture{img_idx}, [], stim_rect_curr)
                        Screen('DrawTexture', win, imgs_texture(img_idx), [], stim_rect_curr)
                        if ctrl_screen
                            Screen('DrawTexture', win_ctrl, imgs_texture_ctrl{img_idx}, [], stim_rect_curr*shrink_factor)
                        end
                    case 'smooth pursuit'
                        xcurr = all_pts{seqidx}(stfridx,1);
                        ycurr = all_pts{seqidx}(stfridx,2);
                        stim_rect_curr = [xcurr-stim_rect_size_x/2 ycurr-stim_rect_size_y/2 xcurr+stim_rect_size_x/2 ycurr+stim_rect_size_y/2];
                        
                        Screen('DrawTexture', win, imgs_texture(img_idx), [], stim_rect_curr)
                        if ctrl_screen
                            Screen('DrawTexture', win_ctrl, imgs_texture_ctrl{img_idx}, [], stim_rect_curr*shrink_factor)
                        end
                        %                 if stimulus_pre_dot
                        %                     Screen('DrawDots', win_ctrl, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                        %                     Screen('DrawDots', win, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                        %                 end
                        
                    case 'images'
                        if isempty(image_displayed{rsvp_ctr,i})
                            image_displayed{rsvp_ctr, i} = img_fnames{img_idx};
                        end
                        
                        if pulse_size
                            pulse_sin = (cos(2*pi*pulse_rate(i)*t_sin) +1)/2;
                            pulse_sin = (pulse_sin + 1)/2; % 50% modulation depth
                            pulse_t_idx = pulse_t_idx + 1;
                            if pulse_t_idx>length(t_sin)
                                pulse_t_idx =1;
                            end
                            rsx_curr = stim_rect_size_x*pulse_sin(pulse_t_idx);
                            rsy_curr = round(rsx_curr/ar(img_idx));
                        else
                            rsx_curr = stim_rect_size_x;
                            rsy_curr = round(rsx_curr/ar(img_idx));
                        end
                        
                        stim_rect_curr = [xcurr-rsx_curr/2 ycurr-rsy_curr/2 xcurr+rsx_curr/2 ycurr+rsy_curr/2];
                        Screen('DrawTexture', win, imgs_texture(img_idx), [], stim_rect_curr)
                        if ctrl_screen
                            Screen('DrawTexture', win_ctrl, imgs_texture_ctrl(img_idx), [], stim_rect_curr*shrink_factor)
                            %Screen('DrawDots', win_ctrl, [eyeposx_cur, eyeposy_cur]*shrink_factor,...
                            %    dotsz, rgba_cols(:,end), [], 1);
                        end
                    case 'moving_dot'
                        t_mov_start = GetSecs();
                        if curr_x>=curr_rect(3)
                            dx_sign = -1;
                        elseif curr_x<=curr_rect(1)
                            dx_sign = 1;
                        else
                            dx_sign = sign(rand-reverse_rate) * dx_sign;
                        end
                        dx_sign = sign(rand-reverse_rate) * dx_sign;
                        
                        if curr_y>=curr_rect(4)
                            dy_sign = -1;
                        elseif curr_y<=curr_rect(2)
                            dy_sign = 1;
                        else
                            dy_sign = sign(rand-reverse_rate) * dy_sign;
                        end
                        dy_sign = sign(rand-reverse_rate) * dy_sign;
                        if size_curr>8 && size_curr<30
                            dsize_sign = sign(rand-reverse_rate_sz) * dsize_sign;
                        elseif size_curr<=8
                            dsize_sign = 1;
                        elseif size_curr>60
                            dsize_sign = -1;
                        end
                        
                        dx = round(rand*scale_fact_move)*dx_sign;
                        dy = round(rand*scale_fact_move)*dy_sign;
                        dsize = round(rand*scale_fact_size)*dsize_sign;
                        mvdot_coords = [curr_x + dx, curr_y + dy];
                        size_curr = size_curr + dsize;
                        Screen('DrawDots', win, mvdot_coords,...
                            size_curr, [0 0 0], [], 1);
                        if ctrl_screen; Screen('DrawDots', win_ctrl, mvdot_coords*shrink_factor,...
                                size_curr*shrink_factor, [0 0 0], [], 1); end
                        
                        curr_x = mvdot_coords(1);
                        curr_y = mvdot_coords(2);
                    case 'movie'
                        t_mov_start = GetSecs();
                        waitforimage = 0;
                        movtex_new = Screen('GetMovieImage', win, moviePtr, waitforimage);
                        if movtex_new>0
                            movtex = movtex_new;
                        end
                        if ctrl_screen
                            Screen('DrawTexture', win_ctrl, movtex, [],  stim_rect*shrink_factor);
                            %Screen('DrawDots', win_ctrl, [eyeposx_cur, eyeposy_cur]*shrink_factor,...
                            %    dotsz, rgba_cols(:,end), [], 1);
                        end
                        Screen('DrawTexture', win, movtex, [],  stim_rect)
                        % for movies, re-draw dots
                        %                         if stimulus_pre_dot && ~stimulus_pre_dot_disappear
                        %                             Screen('DrawDots', win_ctrl, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                        %                             Screen('DrawDots', win, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                        %                         end
                        %mov_dur(end+1) = GetSecs()-t_mov_start;
                    case 'rectangle'
                        Screen('FillRect', win, rect_col, stim_rect);
                        %Screen('DrawDots', win, all_pts(seqidx,:), 5, blackcol, [], 1);
                        if ctrl_screen
                            Screen('FillRect', win_ctrl, rect_col, stim_rect*shrink_factor);
                            Screen('DrawDots', win_ctrl, [eyeposx_cur, eyeposy_cur]*shrink_factor,...
                                dotsz, rgba_cols(:,end), [], 1);
                        end
                    case 'match_to_sample'

                        if rsvp_ctr == 1 % draw a sample
                            %fprintf('sample frame: %d\n', rsvpfridx); 
                            rsx_curr = stim_rect_size_x;
                            rsy_curr = round(rsx_curr/sample_ar(img_idx));
                            stim_rect_curr = [xcurr-rsx_curr/2 ycurr-rsy_curr/2 xcurr+rsx_curr/2 ycurr+rsy_curr/2];
                            if isempty(image_displayed{rsvp_ctr,i})
                                image_displayed{rsvp_ctr, i} = sample_img_fnames{img_idx};
                            end
                            Screen('DrawTexture', win, sample_imgs_texture(img_idx), [], stim_rect_curr)
                            if ctrl_screen
                                Screen('DrawTexture', win_ctrl, sample_imgs_texture_ctrl(img_idx), [], stim_rect_curr*shrink_factor)
                            end

                        else % draw a test
                            rsx_curr = test_img_len;
                            %rsy_curr = round(rsx_curr/test_ar(img_idx));
                            rsy_curr = test_img_ht; 
                            stim_rect_curr = [xcurr-rsx_curr/2 ycurr-rsy_curr/2 xcurr+rsx_curr/2 ycurr+rsy_curr/2];
                            
                            if isempty(image_displayed{rsvp_ctr,i})
                                image_displayed{rsvp_ctr, i} = test_img_fnames{img_idx};
                            end
                            Screen('DrawTexture', win, test_imgs_texture(img_idx), [], stim_rect_curr)
                            if ctrl_screen
                                Screen('DrawTexture', win_ctrl, test_imgs_texture_ctrl(img_idx), [], stim_rect_curr*shrink_factor)
                            end
                        end
                end
                
                
                if ctrl_screen % draw eye position dots
                    Screen('DrawDots', win_ctrl, [eyeposx_cur, eyeposy_cur]*shrink_factor,...
                        dotsz, rgba_cols(:,end), [], 1);
                end
                
                
                if stimulus_pre_dot && ~stimulus_pre_dot_disappear
                    if ~m2s_mode || (m2s_mode && rsvp_ctr == 1)
                        if draw_crosshairs 
                            drawCrosshairs(all_pts(seqidx,:));  
                        else
                            if ctrl_screen; Screen('DrawDots', win_ctrl, all_pts(seqidx,:)*shrink_factor, stim_pre_dot_sz*shrink_factor, whitecol, [], 1); end
                            Screen('DrawDots', win, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                        end
                    elseif m2s_mode && rsvp_ctr > 1
                        
                        if iscell(test_pts)
                            test_pts_curr = test_pts{img_seq(end,curr_tr)}; 
                        else
                            test_pts_curr = test_pts; 
                        end
                        
                        if dots_on_test_grid
                            if draw_crosshairs
                                for tpi = 1:size(test_pts_curr,1)
                                    drawCrosshairs(test_pts_curr(tpi,:));
                                end
                            else
                                if ctrl_screen; Screen('DrawDots', win_ctrl, test_pts_curr'*shrink_factor, stim_pre_dot_sz*shrink_factor, whitecol, [], 1); end
                                Screen('DrawDots', win, test_pts_curr', stim_pre_dot_sz, whitecol, [], 1);
                            end
                        end
                    end
                end
                
            elseif seqidx == 0 % draw wakeup trial
                wu_img_idx = wu_seq(wu_idx);
                if isempty(wake_up_image_displayed{i})
                    wake_up_image_displayed{i} = wu_image_fname_list{wu_img_idx};
                end
                
                % randomly draw from a list of
                Screen('DrawTexture', win, wu_imgs_texture(wu_img_idx), [], stim_rect);
                
                if ctrl_screen
                    Screen('DrawTexture', win_ctrl, wu_imgs_texture_ctrl(wu_img_idx), [], stim_rect*shrink_factor);
                    Screen('DrawDots', win_ctrl, [eyeposx_cur, eyeposy_cur]*shrink_factor,...
                        dotsz, rgba_cols(:,end), [], 1);
                end
                
            elseif inter_rsvp
               %fprintf('inter_rsvp: %d\n', inter_rsvp_fr_ctr); 
                if constant_bg_mode
                    Screen('DrawTexture', win, bg_imgs_texture(bgidx), [], bg_stim_rect_curr)
                    if ctrl_screen
                        Screen('DrawTexture', win_ctrl, bg_imgs_texture_ctrl(bgidx), [], bg_stim_rect_curr*shrink_factor)
                    end
                end
                
                if stimulus_pre_dot && ~stimulus_pre_dot_disappear && ~m2s_mode
                    if draw_crosshairs
                        drawCrosshairs(all_pts(seqidx,:));
                    else
                        if ctrl_screen; Screen('DrawDots', win_ctrl, all_pts(seqidx,:)*shrink_factor, stim_pre_dot_sz*shrink_factor, whitecol, [], 1); end
                        Screen('DrawDots', win, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                    end
                end
                
                if ctrl_screen; Screen('DrawDots', win_ctrl, [eyeposx_cur, eyeposy_cur]*shrink_factor,...
                        dotsz, rgba_cols(:,end), [], 1); end
                inter_rsvp_fr_ctr = inter_rsvp_fr_ctr + 1;
                
            end
            
            checkLick();
            
            if strcmp(reward_on, 'location')
                if location_require_quality
                    qual_check = eyetracker_qual(idx_all)<2;
                else
                    qual_check = true;
                end
                
                if (fix_exp_mode && ~m2s_mode)|| (m2s_mode && rsvp_ctr == 1)
                    if curr_in_bb && qual_check
                        rsvp_break_ctr = 0;
                        % sample command trigger starts
                    elseif (~curr_in_bb || ~qual_check)
                        rsvp_break_ctr = rsvp_break_ctr + 1;
                        %rsvp_break_ctr
                    end
                else
                    if curr_in_bb && qual_check  && ~inter_rsvp
                        entered_bb = true;
                        loop_brk_ctr = loop_brk_ctr + 1;
                        blink_ctr = 0;
                        if isnan(time_to_bb(i))
                            time_to_bb(i) = GetSecs() - t_start_sec;
                        end
                        
                    elseif (~curr_in_bb || ~qual_check) && loop_brk_ctr>3 && blink_ctr<n_frames_blink
                        blink_ctr = blink_ctr + 1;
                        
                    elseif m2s_mode && ~inter_rsvp && punish_incorrect_fix && curr_in_bb_punish 
                        punish_loop_brk_ctr = punish_loop_brk_ctr + 1; 
                    else
                        blink_ctr = 0;
                        loop_brk_ctr = 0;
                        punish_loop_brk_ctr = 0; 
                    end
                end
                
            elseif strcmp(reward_on, 'quality')
                if eyetracker_qual(idx_all)<2
                    loop_brk_ctr = loop_brk_ctr + 1;
                else
                    loop_brk_ctr = 0;
                end
            end
            
            if photodiode_flash % per elias' request, show photodiode square during inter_rsvp
                whichflash = mod(stfridx, 2);
                
                if whichflash
                    if ~inter_rsvp
                        Screen('FillRect', win, whitecol, flash_rect);
                        if ctrl_screen; Screen('FillRect', win_ctrl, whitecol, flash_rect*shrink_factor); end
                    else
                        Screen('FillRect', win, whitecol_interrsvp, flash_rect);
                        if ctrl_screen;Screen('FillRect', win_ctrl, whitecol_interrsvp, flash_rect*shrink_factor); end
                    end
                else
                    if ~inter_rsvp
                        Screen('FillRect', win, blackcol, flash_rect);
                        if ctrl_screen;Screen('FillRect', win_ctrl, blackcol, flash_rect*shrink_factor); end
                    else
                        Screen('FillRect', win, blackcol_interrsvp, flash_rect);
                        if ctrl_screen;Screen('FillRect', win_ctrl, blackcol_interrsvp, flash_rect*shrink_factor);end
                    end
                end
            end
            
            
            % draw informational text into control window:
            drawInfoText();
            drawBoundingBoxes();
            
            %t2_frame = GetSecs();
            %%%%% STIMULUS FRAME FLIP
            vbl = Screen('Flip', win, vbl + halfifi,dontclear, dontsync, multiflip);
            if ctrl_screen; vbl2 = Screen('Flip', win_ctrl, vbl2 + halfifi, dontclear2, dontsync2, multiflip2); end
            %t_dur(stfridx) = t1_frame - t2_frame;
            
            if stfridx == 1
                calib_st_t(i) = vbl-t_start_sec;
            end
            
            if n_rsvp>1 && rsvpfridx == 1 
                rsvp_start_t(rsvp_ctr, i) = vbl-t_start_sec; 
            end
            
            if trig_flag && ~stim_trig_hi && ~inter_rsvp
                %disp('stim trig on');
                IOPort('Write', trig_hand, stim_trig_cmd.on, 1);
                stim_trig_hi = 1;
            elseif trig_flag && stim_trig_hi && inter_rsvp
                %  disp('stim trig off');
                IOPort('Write', trig_hand, stim_trig_cmd.off, 1);
                stim_trig_hi = 0;
            end
            
            % sample command trigger starts
            if trig_flag && ~isempty(sampleCommand_trig_cmd) && ~sampleCommand_trig_hi && stfridx == 1
                IOPort('Write', trig_hand, sampleCommand_trig_cmd.on, 1);
                sampleCommand_trig_hi = 1;
            end

            %clip time
            frame_t = vbl-t_start_sec;
            frame_st_t = [frame_st_t,frame_t];
            frame_clip_num = [frame_clip_num,rsvp_ctr];
            
            if rsvpfridx == 1 || interrsvpfridx == 1
                calib_t_clip(i,clip_ctr) = frame_t;
                clip_ctr = clip_ctr+1;
            end
            
        end
        
        %if fix_exp_mode && rsvpfridx >= stim_frames
        if ((fix_exp_mode && ~m2s_mode) || (m2s_mode && rsvp_ctr == 1)) && rsvpfridx >= stim_frames
            rsvp_end_t(rsvp_ctr, i) = vbl-t_start_sec; 
            rsvp_ctr = rsvp_ctr + 1;
            inter_rsvp_fr_ctr = 1;
            rsvpfridx_old = rsvpfridx;
            rsvpfridx = 0;
        end
        
        % check if we need to end
        if trial_init_timed_out(i)
            end_stim = 1;
        elseif strcmp(trial_mode, 'trial') || seqidx == 0
            if ~fix_exp_mode && seqidx ~= 0 && stfridx >= stim_frames % finished successfully
                end_stim = 1;
            elseif fix_exp_mode && seqidx ~= 0 && stfridx >= stim_frames && rsvp_ctr > n_rsvp % finished successfully
                end_stim = 1;
            elseif seqidx == 0 && stfridx >= wake_up_stim_frames
                end_stim = 1;
            end
        elseif strcmp(trial_mode, 'foraging')
            if (fix_exp_mode && ~m2s_mode)|| (m2s_mode && rsvp_ctr == 1)
                if rsvp_ctr > n_rsvp
                    end_stim = 1;
                    calib_t_clip(i,clip_ctr) = GetSecs() - t_start_sec; % inter_rsvp
                elseif rsvp_break_ctr >= round(rsvp_break_after_t/1e3/ifi)
                    end_stim = 1;
                    fprintf('TRIAL BREAK: FIXATION BROKEN, %0.1f s \n', GetSecs() - calib_st_t(i) - t_start_sec);
                end
            else
                %if m2s_mode && rsvp_ctr > 1 stfridx >= round(time_out_after/1e3/ifi)
                if loop_brk_ctr >= round(time_to_reward/1e3/ifi)
                    end_stim = 1;
                    calib_t_clip(i,clip_ctr) = GetSecs() - t_start_sec;
                elseif m2s_mode && punish_incorrect_fix && (punish_loop_brk_ctr >= round(time_to_reward/1e3/ifi))
                    end_stim = 1; 
                    incorrect_choice_m2s(i) = in_punish_bb; 
                    fprintf('TRIAL BREAK: FIXATION IN INCORRECT BOX: %d\n', in_punish_bb); 
                elseif ~(curr_in_bb && qual_check) && stfridx >= round(foraging_max_t/1e3/ifi) && ~entered_bb
                    end_stim = 1;
                    fprintf('TRIAL BREAK: TIMED OUT, %0.1f s \n', GetSecs() - calib_st_t(i) - t_start_sec);
                elseif ~(curr_in_bb && qual_check) && stfridx >= round(foraging_max_t/1e3/ifi) && entered_bb
                    end_stim = 1;
                    fprintf('TRIAL BREAK: FIXATION BROKEN, %0.1f s \n', GetSecs() - calib_st_t(i) - t_start_sec);
                end
            end
        end
        
        [~,~,kCode] = KbCheck(0);
        if find(kCode) == esc_key
            disp('ESC key recognized, exiting');
            exit_flag = 1;
            calib_end_t(i) =  GetSecs()-t_start_sec;
            trial_aborted(i) = true;
            break
        elseif find(kCode) == rew_end_key
            manual_reward_t(i) =  GetSecs()-t_start_sec;
            man_reward_ct = man_reward_ct + 1;
            fprintf('manual reward with trial end! time = %0.3fs\n', manual_reward_t(i));
            end_stim = 1;
            manual_reward_flag = true;
        elseif find(kCode) == bonus_rew_key
            if  isempty(bonus_reward_t) || GetSecs() - bonus_reward_t(end) - t_start_sec > bonus_reward_min_wait
                bonus_reward_t(end+1) = GetSecs() - t_start_sec;
                %man_reward_ct = man_reward_ct + 1;
                bonus_reward_ct = bonus_reward_ct + 1;
                pumpReward_updateGUI();
                fprintf('bonus reward, time = %0.3fs\n', bonus_reward_t(end));
            end
        end
        
        if end_stim
            Screen('FillRect', win, bg_col_val);
            
            if ctrl_screen; Screen('FillRect', win_ctrl, bg_col_val); end
            
            if seqidx ~= 0 && give_rewards && color_shift_feedback && loop_brk_ctr > 0
                Screen('FillRect', win, col_shift_change(:,loop_brk_ctr), curr_bb);
                if ctrl_screen; Screen('FillRect', win_ctrl, col_shift_change(:,loop_brk_ctr), curr_bb); end
            end
            drawInfoText();
            drawBoundingBoxes();
            if draw_retain_bb_pts
                drawGoodEyePts(0);
            end
            vbl = Screen('Flip', win, vbl + halfifi, dontclear, dontsync, multiflip);
            if ctrl_screen; vbl2 = Screen('Flip', win_ctrl, vbl2 + halfifi, dontclear2, dontsync2, multiflip2); end
            
            % if ending, note the time:
            %calib_end_t(i) = GetSecs()-t_start_sec;
            calib_end_t(i) = vbl - t_start_sec;         
            if n_rsvp>1
                rsvp_end_t(rsvp_ctr, i) = vbl-t_start_sec;
            end
            fprintf('trial dur: %0.3f\n', calib_end_t(i)-calib_st_t(i)); 
            
            if trig_flag
                %disp('stim trig off');
                IOPort('Write', trig_hand, stim_trig_cmd.off, 1);
                stim_trig_hi = 0;
                % turn off trial trigger
                IOPort('Write', trig_hand, trial_trig_cmd.off, 1);
            end
            
            if give_rewards && ~skip_reward
                if strcmp(trial_mode, 'trial') || seqidx == 0
                    if strcmp(reward_on, 'quality') || (seqidx == 0 && ~eye_method_mouse)
                        frac_good = sum(eye_data_qual_curr < 2)/numel(eye_data_qual_curr);
                    elseif strcmp(reward_on, 'location') || (seqidx == 0 && eye_method_mouse)
                        frac_good = sum(eye_in_bb_curr)/numel(eye_in_bb_curr);
                    end
                    
                    if frac_good > reward_thresh
                        reward_this_trial = true;
                    else
                        reward_this_trial = false;
                    end
                    
                elseif strcmp(trial_mode, 'foraging')
                    if (fix_exp_mode && ~m2s_mode) || (m2s_mode && rsvp_ctr == 1)
                        if rsvp_break_ctr >= round(rsvp_break_after_t/1e3/ifi)
                            reward_this_trial = false;
                        else
                            reward_this_trial = true;
                        end
                        
                    else
                        if loop_brk_ctr >= round(time_to_reward/1e3/ifi)
                            reward_this_trial = true;
                        else
                            reward_this_trial = false;
                        end
                    end
                end
                
                if manual_reward_flag
                    reward_this_trial = true;
                end
                
                if reward_this_trial
                    correct_trial(i) = true;
                    if lick_flag
                        %%% check trials since last lick
                        if isnan(rew_trs_since_lick) && ~lick_trial(i)
                            %rew_trs_since_lick  = i;
                            rew_trs_since_lick = sum(correct_trial);
                        elseif lick_trial(i)
                            rew_trs_since_lick = 0;
                        elseif i == 1 && ~lick_trial(i)
                            rew_trs_since_lick = 1;
                        else
                            rew_trs_since_lick = rew_trs_since_lick + 1;
                        end
                    end
                    
                    if ~manual_reward_flag
                        reward_ct = reward_ct + 1;
                    end
                    % fprintf('****REWARD TRIAL %d, %0.2f good\n', i, frac_good);
                    
                    %fprintf('****REWARD TRIAL %d\n', i);
                    
                    if play_reward_sound
                        %sound(aud_y, aud_fs);
                        t1 = PsychPortAudio('Start',audio_handle(1), [], 0,1);
                        %disp('sound played');
                        %PsychPortAudio('Stop', audio_handle(1));
                    end
                    
                    if trig_flag && sampleCommand_trig_hi && ~isempty(sampleCommand_trig_cmd)
                        IOPort('Write', trig_hand, sampleCommand_trig_cmd.off, 1);
                        sampleCommand_trig_hi = 0;
                    end
                    
                    if ~isempty(reward_pumphand)
                        if ~lick_flag || (lick_flag && i == 1) || (lick_flag && rew_trs_since_lick < stop_rewards_n_nolicks) || manual_reward_flag
                            reward_time(i) = GetSecs()-t_start_sec;
                            reward_trial(i) = true;
                            if reward_serial
                                pumpReward_updateGUI();
                                %writeline(reward_pumphand, 'RUN');
                                
                                %WaitSecs(reward_on_dur);
                            else
                                reward_pumphand.digitalWrite(reward_arduino_pin, 1);
                                WaitSecs(reward_on_dur);
                                reward_pumphand.digitalWrite(10, 0);
                            end
                            %t1_rew = GetSecs();
                            %set(reward_today_hand, 'String', sprintf('%0.3f mL', (reward_ct + man_reward_ct)*reward_vol + start_reward_vol));
                            %drawnow;
                            t2_rew = GetSecs();
                            %pumpReward_updateGUI();
                        elseif lick_flag && rew_trs_since_lick >= stop_rewards_n_nolicks
                            fprintf('No reward given: %d trials since lick\n', rew_trs_since_lick);
                        end
                        
                    end
                    
                    for w = 1:wait_after_rew_frames
                        drawInfoText();
                        drawBoundingBoxes();
                        if draw_retain_bb_pts
                            drawGoodEyePts(0);
                        end
                        checkLick();
                        
                        vbl = Screen('Flip', win, vbl + halfifi, dontclear, dontsync, multiflip);
                        if ctrl_screen; vbl2 = Screen('Flip', win_ctrl, vbl2 + halfifi, dontclear2, dontsync2, multiflip2); end
                    end
                    %WaitSecs(wait_after_reward);
                    
                else
                    fprintf('no reward: trial %d\n', i);
                end
            end
            
            if give_punishments && (give_rewards && ~reward_this_trial) && ~skip_reward
                punish_trial(i) = true;
                %fprintf('Punish on\n');
                if play_punish_sound
                    %sound(aud_pun_y, aud_pun_fs);
                    t1 = PsychPortAudio('Start', audio_handle(2), 1, 0,1);
                    %disp('sound played');
                    % PsychPortAudio('Stop', audio_handle(2));
                end
                n_pun_frames = round(punish_length_ms/1e3/ifi); 
                for p = 1:n_pun_frames
                    Screen('FillRect', win, blackcol);
                    if ctrl_screen; Screen('FillRect', win_ctrl, blackcol); end
                    drawInfoText(1);
                    drawBoundingBoxes();
                    checkLick();
                    if draw_retain_bb_pts
                        drawGoodEyePts(0);
                    end
                    vbl = Screen('Flip', win, vbl + halfifi, dontclear, dontsync, multiflip);
                    if ctrl_screen; vbl2 = Screen('Flip', win_ctrl, vbl2 + halfifi, dontclear2, dontsync2, multiflip2); end
                    %fprintf('pun frame %d of %d\n', p, n_pun_frames); 
                end
                
                if trig_flag && sampleCommand_trig_hi && ~isempty(sampleCommand_trig_cmd)
                    IOPort('Write', trig_hand, sampleCommand_trig_cmd.off, 1);
                    sampleCommand_trig_hi = 0;
                end
            end
        end
    end
    
    if seqidx == 0
        wu_idx = wu_idx + 1;
        if wu_idx > n_wu_imgs
            wu_idx = 1;
        end
    end
    
    if exit_flag
        break
    end
    
    if profile_memory
        [u,m] = memory;
        mem_used(i) = u.MemUsedMATLAB;
        avail_sys_mem(i) = m.PhysicalMemory.Available;
        avail_phys_mem(i) = m.SystemMemory.Available;
    end
    
    calib_frame_st_t{i} = frame_st_t;
    calib_frame_clip_num{i} = frame_clip_num;
end

Screen('FillRect', win, bg_col_val);
if ctrl_screen; Screen('FillRect', win_ctrl, bg_col_val); end
for j = 1:iti_frames*2
    vbl = Screen('Flip', win);
    if ctrl_screen; vbl2 = Screen('Flip', win_ctrl); end
end

t2 = GetSecs;
Screen('CloseAll');
PsychPortAudio('Close');
sca

if trig_flag
    % turn off session trigger
    IOPort('Write', trig_hand, sess_trig_cmd.off, 1);
    WaitSecs(0.05);
    IOPort('Write', trig_hand, sess_trig_cmd.off, 1);
    % ensure everything else is off
    IOPort('Write', trig_hand, trial_trig_cmd.off, 1);
    IOPort('Write', trig_hand, stim_trig_cmd.off, 1);
    if ~isempty(sampleCommand_trig_cmd)
        IOPort('Write', trig_hand, sampleCommand_trig_cmd.off,1);
    end
end

if profile_memory
    assignin('base', 'mem_used', mem_used);
    assignin('base', 'avail_sys_mem', avail_sys_mem);
    assignin('base', 'avail_phys_mem', avail_phys_mem);
end

%% SAVE: gather data for save into structures
calib.n_pts_x = n_pts_x;
calib.n_pts_y = n_pts_y;
calib.trs_per_location = trs_per_location;
if any(trial_aborted)
    calib.n_completed = find(trial_aborted) - 1;
else
    calib.n_completed = i;
end
calib.pts = all_pts;
% calib_st_t(isnan(calib_st_t)) = [];
% calib_end_t(isnan(calib_end_t)) = [];
calib.start_t = calib_st_t;
calib.end_t = calib_end_t;
if n_rsvp>1
    calib.rsvp_start_t = rsvp_start_t; 
    calib.rsvp_end_t = rsvp_end_t; 
end
calib.time_to_bounding_box = time_to_bb;
calib.sequence = tr_seq;
calib.clip_sequence_start_t = calib_t_clip;
calib.frame_st_t = calib_frame_st_t;
calib.frame_clip_num = calib_frame_clip_num;
calib.trial_mode = trial_mode;
calib.reward_on = reward_on;
calib.time_to_reward = time_to_reward;
calib.time_out_after = time_out_after; % 9-17 YJ
calib.gaze_offset = [gaze_offset_x,gaze_offset_y]; % 9-17 YJ
calib.trial_aborted = trial_aborted;
calib.stim_pre_start = stim_pre_start;
calib.stim_pre_end = stim_pre_end;
calib.trial_init_timed_out = trial_init_timed_out;
calib.n_wake_up_trs = n_wake_up_trs;
calib.image_displayed = image_displayed;
calib.wake_up_image_displayed = wake_up_image_displayed;
calib.wake_up_movie_start_t = wu_mov_start_t;
calib.wake_up_movie_end_t = wu_mov_end_t;
calib.training_notes_objectives = training_notes_str;
calib.incorrect_choice_m2s = incorrect_choice_m2s; 

calib_settings.disp_rect = win_rect;
calib_settings.presentation_time = presentation_time;
calib_settings.rsvp_break_after_t = rsvp_break_after_t;
calib_settings.inter_stim_interval = inter_stim_interval;
calib_settings.iti_random = iti_random;
calib_settings.iti_frames = iti_frames;
calib_settings.position_order = position_order;
calib_settings.bg_col = bg_col_val;
calib_settings.stim_rect_size_x = stim_rect_size_x;
calib_settings.stim_rect_size_y = stim_rect_size_y;
calib_settings.stim_mode = stim_mode;
calib_settings.img_folder = img_folder;
calib_settings.img_folder_idx = img_folder_idx;
calib_settings.pulsed_img_size = pulse_size;
if strcmp(stim_mode, 'images') || strcmp(stim_mode,'spinning') ||strcmp(stim_mode,'smooth pursuit')
    calib_settings.img_list = img_fnames;
else
    calib_settings.img_list = [];
end
calib_settings.img_seq = img_seq;
calib_settings.n_rsvp = n_rsvp;
calib_settings.rsvp_iti_ms = rsvp_iti_t;

if strcmp(stim_mode,'spinning')  ||strcmp(stim_mode,'smooth pursuit')
    calib_settings.diam_circle = diam_circle;
    calib_settings.speed = speed;
end

calib_settings.bounding_box_x = unique(bounding_rects(:,3)-bounding_rects(:,1));
calib_settings.bounding_box_y =unique(bounding_rects(:,4)-bounding_rects(:,2));
calib_settings.bounding_box_x_stim_pre_dot = unique(bounding_rects_stim_pre_dot(:,3)-bounding_rects_stim_pre_dot(:,1));
calib_settings.bounding_box_y_stim_pre_dot =unique(bounding_rects_stim_pre_dot(:,4)-bounding_rects_stim_pre_dot(:,2));
calib_settings.draw_gaze_pt = draw_gaze_pt;
calib_settings.n_gaze_pts_draw = n_gaze_pts_draw;
calib_settings.gaze_pt_dot_col = gaze_pt_dot_col;
calib_settings.gaze_pt_dot_sz = gaze_pt_sz;
calib_settings.stimulus_pre_dot = stimulus_pre_dot;
calib_settings.stimulus_pre_time = stimulus_pre_time;
if draw_crosshairs 
    calib_settings.stim_pre_dot_sz = [];
else
    calib_settings.stim_pre_dot_sz = stim_pre_dot_sz;
end
calib_settings.draw_crosshairs = draw_crosshairs; 
calib_settings.crosshair_sz = crosshair_sz; 
calib_settings.crosshair_lw = crosshair_lw; 

calib_settings.stimulus_pre_dot_disappear = stimulus_pre_dot_disappear;

calib_settings.expt_params = expt_params;
calib_settings.wake_up_trials = wake_up_trials;
calib_settings.wake_up_every = wake_up_every;
calib_settings.wake_up_img_dir = wake_up_img_fold;
calib_settings.wake_up_stim_size_x = wake_up_stim_size_x;
calib_settings.wake_up_stim_size_x = wake_up_stim_size_y;
calib_settings.wake_up_tr_dur = wake_up_tr_dur;
calib_settings.photodiode_on = photodiode_flash;
calib_settings.photodiode_sz = flash_rect_size;
calib_settings.photodiode_pos = [flash_rect(1),flash_rect(2)];

if isfield(s,'wake_up_movie')
    calib_settings.wake_up_movie = s.wake_up_movie;
end

% offset in y
calib_settings.gaze_center_adj_y = gaze_center_adj_y;
calib_settings.gaze_center_adj_x = gaze_center_adj_x;
calib_settings.apply_gaze_center_adj = apply_gaze_center_adj;

calib_settings.require_fix_tr_init = require_fix_tr_init;
calib_settings.fixation_to_init = fixation_to_init;
calib_settings.time_out_trial_init_s = time_out_trial_init_s;
calib_settings.clip_sequence = clip_sequence;
calib_settings.clip_sequence_t = clip_sequence_t;


if m2s_mode
    calib_settings.test_bb_size_x = test_bb_size_x;
    calib_settings.test_bb_size_y = test_bb_size_y;
    calib_settings.test_img_len = test_img_len;
    calib_settings.test_img_ht = test_img_ht;
    calib_settings.dots_on_test_grid = dots_on_test_grid;
    calib_settings.task_folder = task_folder;
    calib_settings.task_params = tf; 
    calib_settings.taskfile = taskfile; 
    calib_settings.test_bounding_rects = test_bounding_rects; 
end

eyetrack.time = eyetracker_time;
eyetrack.x = eyepos_x;
eyetrack.y = eyepos_y;
eyetrack.method = eye_method;
eyetrack.eyepos_raw = eyepos_raw;
eyetrack.calib_applied = apply_calib;

if ~isempty(cProj)
    eyetrack.cProj = cProj;
end
eyetrack.cX = cX;
eyetrack.cY = cY;
eyetrack.quality = eyetracker_qual;
eyetrack.pupil_size_x =  pupil_size_x;
eyetrack.pupil_size_y =  pupil_size_y;

settings.eyetracker_toolbox_dir = eyetracker_toolbox_dir;
settings.save_data_dir = save_data_dir;
settings.subject = subject;
settings.window_rect = window_rect;
settings.skip_sync_tests = skip_sync_tests ;
settings.computer_type = computer;
[~, hname] = system('hostname');
settings.hostname = strtrim(hname);
[~, WindowsVersion] = system('ver');
settings.osversion = strtrim(WindowsVersion);
settings.time_save = datestr(now, 'yyyy-dd-mm_HH-MM-SS');
settings.time_start = session_time;
settings.run_time = GetSecs() - t_start_sec;
settings.ifi_monitor = ifi;
settings.screenScale = setup_config.screenScale;
settings.screenPixels = setup_config.screenPixels;
settings.screenPhysicalPixels = setup_config.screenPhysicalPixels;
settings.screenInches = setup_config.screenInches;
settings.viewportPPI = setup_config.viewportPPI;
settings.deviceBrand = setup_config.deviceBrand;
settings.devicePixelRatio = setup_config.devicePixelRatio;

reward.give_rewards = true;
reward.reward_type = reward_type;
reward.n_rewards_given = sum(reward_trial);
reward.reward_sequence = reward_trial;
reward.correct_trial = correct_trial;
reward.reward_vol = reward_vol;
reward.arduino_pin_reward = reward_arduino_pin;
reward.reward_thresh_frac_good_glint = reward_thresh;
reward.pulse_on_dur = reward_on_dur;
reward.wait_after_reward = wait_after_reward;
reward.reward_time = reward_time;
reward.play_sound = play_reward_sound;
reward.sound_file = reward_sound_file;
reward.manual_reward_t = manual_reward_t;
reward.man_reward_ct = man_reward_ct;
reward.bonus_reward_ct = bonus_reward_ct;
reward.bonus_reward_t = bonus_reward_t;
reward.lickometer = lick_flag;
reward.lick_times_all = lick_times_all(~isnan(lick_times_all));
reward.lick_arduino_pin = lick_arduino_pin;
reward.lick_trial = lick_trial;

punish.give_punishments = give_punishments;
punish.n_puns_given = sum(punish_trial);
punish.pun_sequence = punish_trial;
punish.length_ms = punish_length_ms;
punish.play_sound = play_punish_sound;
punish.sound_file = punish_sound_file;

%% SAVE: save data
savefname = sprintf('calib_%s_%s.mat', subject, session_time);
settings.data_file_name = savefname;
save_full = fullfile(save_data_dir, savefname); % output variable
save(fullfile(save_data_dir, savefname), 'calib', 'calib_settings', 'eyetrack', 'settings', 'reward', 'punish');
fprintf('session data saved to %s\n', fullfile(save_data_dir, savefname));

% try to save to remote
try
    save(fullfile(save_data_dir_remote, savefname), 'calib', 'calib_settings', 'eyetrack', 'settings', 'reward', 'punish');
    fprintf('remote save: session data saved to %s\n', fullfile(save_data_dir_remote, savefname));
    % convert to mkTurk style data and save
    if mkTurk_data_save
        convert_to_mkTurkData(save_data_dir_remote,fullfile(save_data_dir, savefname));
    end
catch me
    disp(me);
    disp('failed to save session data to remote');
end

%% SAVE: add to subject log table automatically
log_dir = 'C:/MATLAB/Marmulator/subject_logs_bysessions';
logData_bysession(log_dir,fullfile(save_data_dir, savefname));

%% SUPPORT FUNCTIONS: nested
    function [eyeposx_cur, eyeposy_cur, eye_data_qual] = get_eyetracker_draw_dots()
        switch eye_method
            case 'pupil-glint'
                [eyepos_x_tmp, eyepos_y_tmp] = vpx_GetDiffVector(eye);
                eye_data_qual = vpx_GetDataQuality(eye);
                [pupil_size_x(idx_all), pupil_size_y(idx_all)] = vpx_GetPupilSize(eye);
            case 'pupil'
                [eyepos_x_tmp, eyepos_y_tmp] = vpx_GetPupilPoint(eye);
                eye_data_qual = vpx_GetDataQuality(eye);
                [pupil_size_x(idx_all), pupil_size_y(idx_all)] = vpx_GetPupilSize(eye);
            case 'gaze_point'
                [eyepos_x_tmp, eyepos_y_tmp] = vpx_GetGazePoint(eye);
                eye_data_qual = vpx_GetDataQuality(eye);
                [pupil_size_x(idx_all), pupil_size_y(idx_all)] = vpx_GetPupilSize(eye);
            case 'gaze_point_corrected'
                [eyepos_x_tmp, eyepos_y_tmp] = vpx_GetGazePointCorrected(eye);
                eye_data_qual = vpx_GetDataQuality(eye);
                [pupil_size_x(idx_all), pupil_size_y(idx_all)] = vpx_GetPupilSize(eye);
            case 'rand'
                eyepos_x_tmp = rand();
                eyepos_y_tmp = rand();
                eye_data_qual = NaN;
            case 'mouse'
                [eyepos_x_tmp, eyepos_y_tmp] = GetMouse();
                eye_data_qual = NaN;
        end
        eyetracker_qual(idx_all) = eye_data_qual;
        
        if apply_calib
            if ~isempty(cProj) % PROJECTIVE
                raw_eye = [eyepos_x_tmp,eyepos_y_tmp];
                T = [cX,cY,cProj];
                transformed_eye = homography_transform(raw_eye,T);
                eyepos_x(idx_all) = transformed_eye(:,1);
                eyepos_y(idx_all) = transformed_eye(:,2);
            else
                if length(cX) == 4
                    eyepos_x(idx_all) = eyepos_x_tmp * cX(2) + eyepos_y_tmp * cX(3) + eyepos_x_tmp.*eyepos_y_tmp * cX(4) + cX(1);
                elseif length(cX) == 3
                    eyepos_x(idx_all) = eyepos_x_tmp*cX(2) + eyepos_y_tmp*cX(3) + cX(1);
                elseif length(cX) == 2
                    eyepos_x(idx_all) = eyepos_x_tmp*cX(2) + cX(1);
                end
                
                if length(cY) ==4
                    eyepos_y(idx_all) = eyepos_y_tmp * cY(2) + eyepos_x_tmp * cY(3) + eyepos_x_tmp.*eyepos_y_tmp * cY(4) + cY(1);
                elseif length(cY) == 3
                    eyepos_y(idx_all) = eyepos_y_tmp*cY(2) + eyepos_x_tmp*cY(3) + cY(1);
                elseif length(cY) == 2
                    eyepos_y(idx_all) = eyepos_y_tmp*cY(2) + cY(1);
                end
            end
        else
            if ~isempty(gaze_offset_x)
                eyepos_x(idx_all) = eyepos_x_tmp + gaze_offset_x;
            else
                eyepos_x(idx_all) = eyepos_x_tmp;
            end
            if ~isempty(gaze_offset_y)
                eyepos_y(idx_all) = eyepos_y_tmp + gaze_offset_y;
            else
                eyepos_y(idx_all) = eyepos_y_tmp;
            end
        end
        
        eyepos_raw(idx_all,:) = [eyepos_x_tmp, eyepos_y_tmp];
        eyetracker_time(idx_all) = GetSecs()-t_start_sec;
        
        draw_idx1 = max(idx_all - n_gaze_pts_draw + 1, 1);
        draw_idx2 = idx_all;
        if apply_calib || strcmp(eye_method, 'mouse')
            draw_coords = [eyepos_x(draw_idx1:draw_idx2); eyepos_y(draw_idx1:draw_idx2)];
        else
            draw_coords = ret_dots_in_rect(eyepos_x(draw_idx1:draw_idx2),eyepos_y(draw_idx1:draw_idx2),win_rect, eye_rect);
        end
        
        if draw_gaze_pt && ctrl_screen
            Screen('DrawDots', win_ctrl, draw_coords*shrink_factor,...
                dotsz, rgba_cols(:,(end-min(idx_all, n_gaze_pts_draw)+1):end), [], 1);
        end
        eyeposx_cur = draw_coords(1,end);
        eyeposy_cur = draw_coords(2,end);
    end

    function drawInfoText(duringPun)
        if ctrl_screen
            if nargin == 0
                duringPun = 0;
                textcol = [0 0 0];
            end
            if duringPun
                textcol = [255 255 255];
            else
                textcol = [0 0 0];
            end
            
            if strcmp(trial_mode, 'foraging')
                superstr = sprintf('mode: %s, %s; tr %d of %d, tr time: %dms\n%d auto reward, %d manual (%s), %d trs since lick\neyetracking: %s\nforaging: time out %dms, reward %dms, correct time: %dms', ...
                    stim_mode, trial_mode, curr_tr, n_trs_tot, round(stfridx*ifi*1e3), reward_ct, man_reward_ct, reward_on, rew_trs_since_lick, ...
                    eye_method, time_out_after, time_to_reward, round(loop_brk_ctr*ifi *1e3));
            else
                superstr = sprintf('mode: %s, %s; tr %d of %d, tr time: %dms\n%d auto reward, %d manual (%s), %d trs since lick\neyetracking: %s', ...
                    stim_mode, trial_mode, curr_tr, n_trs_tot, round(stfridx*ifi*1e3), reward_ct, man_reward_ct, reward_on, rew_trs_since_lick, ...
                    eye_method);
            end
            
            DrawFormattedText(win_ctrl, superstr, 80, win_rect_ctrl(4)-150, textcol);
            Screen('FrameRect', win_ctrl, [0 0 0], win_rect_ctrl, 1);
            
            if ~isempty(gaze_offset_x)
                calib_superstr = sprintf('gaze offset x, y: %0.3f, %0.3f', gaze_offset_x, gaze_offset_y);
            else
                if length(cX) ==4
                    calib_superstr = sprintf('calib x; y: %0.1f, %0.1f, %0.1f %0.1f; %0.1f, %0.1f, %0.1f %0.1f', ...
                        cX(1), cX(2), cX(3),cX(4), cY(1), cY(2), cY(3),cY(4));
                elseif length(cX) == 3
                    if ~isempty(cProj)
                        calib_superstr = sprintf('calib x; y: %0.1f, %0.1f, %0.1f; %0.1f, %0.1f, %0.1f\ncalib proj: %0.1f, %0.1f', ...
                            cX(1), cX(2), cX(3), cY(1), cY(2), cY(3), cProj(1), cProj(2));
                    else
                        calib_superstr = sprintf('calib x; y: %0.1f, %0.1f, %0.1f; %0.1f, %0.1f, %0.1f', ...
                            cX(1), cX(2), cX(3), cY(1), cY(2), cY(3));
                    end
                elseif length(cX) == 2
                    calib_superstr = sprintf('calib x; y: %0.1f, %0.1f; %0.1f, %0.1f', ...
                        cX(1), cX(2), cY(1), cY(2));
                end
            end
            
            DrawFormattedText(win_ctrl, calib_superstr, 80, win_rect_ctrl(4)-180, textcol);
            
            punish_superstr = sprintf('give time-out punish: %d, dur: %dms, nr punish: %d\ntr start requires fixation: %d, pre-stim time in box: %dms\ntime pre-stim total: %dms\n', give_punishments, punish_length_ms_draw,...
                sum(punish_trial), require_fix_tr_init,  round(fix_pre_fr_ctr*ifi*1e3), round(pre_stim_timer*1e3));
            DrawFormattedText(win_ctrl, punish_superstr, 700, win_rect_ctrl(4) - 150, textcol);
            
        end
    end

    function drawBoundingBoxes(curr_bb_stim_pre_dot)
        if ctrl_screen
            if show_curr_bounding_box
                if nargin>0 && seqidx ~= 0
                    Screen('FrameRect',win_ctrl,cols(seqidx+1,:), curr_bb_stim_pre_dot*shrink_factor,5);
                elseif nargin>0 && seqidx == 0
                    Screen('FrameRect',win_ctrl,whitecol,curr_bb_stim_pre_dot*shrink_factor,5);
                elseif nargin ==0 && seqidx ~= 0
                    Screen('FrameRect', win_ctrl, cols(seqidx+1,:), curr_bb*shrink_factor, 5);
                else
                    Screen('FrameRect', win_ctrl, whitecol, curr_bb*shrink_factor, 5);
                end
            end
            if show_all_bounding_boxes && seqidx ~= 0
                if  m2s_mode
                    if rsvp_ctr>1
                        Screen('FrameRect', win_ctrl, cols(2:end,:)', test_bounding_rects_curr'*shrink_factor, 1);
                    end
                else
                    Screen('FrameRect', win_ctrl, cols(2:end,:)', bounding_rects'*shrink_factor, 1);
                    Screen('FrameRect',win_ctrl,cols(2:end,:)',bounding_rects_stim_pre_dot'*shrink_factor,1);
                end
            end
        end
    end

    function drawGoodEyePts(during_stim)
        
        if during_stim % then add the last coordinate if it's in a bounding box
            last_coord = ret_dots_in_rect( eyepos_x(idx_all), eyepos_y(idx_all),win_rect,eye_rect);
            last_coord = last_coord';
            if IsInRect(last_coord(1), last_coord(2), curr_bb)
                good_pts(:,end+1) = last_coord;
                good_pts_cols(:,end+1) = [cols(seqidx+1,:) 0.3]';
            end
        end
        if ~isempty(good_pts) && ctrl_screen
            Screen('DrawDots', win_ctrl, good_pts*shrink_factor, dotsz-5, good_pts_cols, [], 1);
        end
    end

    function lick = checkLick()
        if lick_flag
            % check for licks
            IOPort('Write', lick_hand, read_lick_cmd, 2);
            lickout = IOPort('Read', lick_hand,1,3);
            lick = lickout(1)-48;
            
            if ctrl_screen
                if lick
                    Screen('FillRect', win_ctrl, [0, 142, 204], lick_rect)
                else
                    Screen('FillRect', win_ctrl, [30, 30, 30], lick_rect)
                end
            end
            
            if lick
                lick_idx = lick_idx + 1;
                lick_times_all(lick_idx) = GetSecs()- t_start_sec;
                if ~lick_trial(curr_tr)
                    fprintf('lick tr %d\n', curr_tr);
                end
                lick_trial(curr_tr) = true;
                
                if correct_trial(curr_tr)
                    rew_trs_since_lick = 0;
                end
                
            end
        end
    end

    function pumpReward_updateGUI()
        writeline(reward_pumphand, 'RUN');
        set(reward_today_hand, 'String', sprintf('%0.3f mL', (reward_ct + man_reward_ct + bonus_reward_ct)*reward_vol + start_reward_vol));
        drawnow;
    end

    function drawCrosshairs(pt)
        fromH = pt(1)-crosshair_sz;
        toH = pt(1)+crosshair_sz;
        fromV = pt(2)-crosshair_sz;
        toV = pt(2)+crosshair_sz;
        Screen('DrawLine', win, whitecol, pt(1), fromV, pt(1), toV, crosshair_lw);
        Screen('DrawLine', win, whitecol, fromH, pt(2), toH, pt(2), crosshair_lw);
        if ctrl_screen
            Screen('DrawLine', win_ctrl, whitecol, pt(1)*shrink_factor, fromV*shrink_factor, pt(1)*shrink_factor, toV*shrink_factor, crosshair_lw);
            Screen('DrawLine', win_ctrl, whitecol, fromH*shrink_factor, pt(2)*shrink_factor, toH*shrink_factor, pt(2)*shrink_factor, crosshair_lw);
        end
    end
            

end

%% SUPPORT FUNCTIONs: subfunctions
function draw_coords = ret_dots_in_rect(x,y,disp_rect,eye_rect)

disp_ht = disp_rect(4)-disp_rect(2);
disp_len = disp_rect(3)-disp_rect(1);

eyewin_xoffset = eye_rect(1) * disp_len;
eyewin_yoffset = eye_rect(2) * disp_ht;
eyewin_ht = (eye_rect(4) - eye_rect(2))*disp_ht;
eyewin_len = (eye_rect(3) - eye_rect(1))*disp_len;

% scale x
draw_coords(1,:) = round((x*eyewin_len) + eyewin_xoffset);

% scale y
draw_coords(2,:) = round(((y)*eyewin_ht) + eyewin_yoffset);
end

function [imgs, im_ht, im_wd, ar] = load_images(img_fnames)
% takes a list of images, reads them, returns image matrices, heights,
% widths, aspect ratios
% used in image mode and match_to_sample mode
fprintf('Reading %d stim images\n', numel(img_fnames));
imgs = cell(1,numel(img_fnames));
t_pre_imread = GetSecs();
for i = 1:numel(img_fnames)
    [img_tmp, ~,alpha] = imread(img_fnames{i});
    fprintf('.')
    if ~isempty(alpha)
        img_tmp(:,:,4) = alpha;
    end
    imgs{i} = img_tmp;
    im_ht(i) = size(imgs{i}, 1);
    im_wd(i) = size(imgs{i}, 2);
    ar(i) = im_wd(i)/im_ht(i);
    if mod(i, 100) == 0
        fprintf('\n');
    end
end
t_post_imread = GetSecs();
fprintf('\nThat took %0.2f s\n', t_post_imread - t_pre_imread);
end

function [img_fnames, img_folder_idx] = ret_image_fnames(img_folder)

EXT = {'.tiff','.png','.jpg','.jpeg'}; % Acceptable extensions
if iscell(img_folder)
    img_fnames = [];
    img_folder_idx = [];
    for i = 1:numel(img_folder)
        img_fnames_tmp = dir(img_folder{i});
        [~,~,e] = fileparts({img_fnames_tmp.name});
        img_fnames_tmp = img_fnames_tmp(matches(e,EXT,'IgnoreCase',1));
        img_fnames_tmp = fullfile(img_folder{i}, {img_fnames_tmp.name});
        img_fnames = [img_fnames img_fnames_tmp];
        img_folder_idx = [img_folder_idx, i* ones(1,length(img_fnames_tmp))];
    end
else
    img_fnames_tmp = dir(img_folder);
    [~,~,e] = fileparts({img_fnames_tmp.name});
    img_fnames_tmp = img_fnames_tmp(matches(e,EXT,'IgnoreCase',1));
    img_fnames = fullfile(img_folder, {img_fnames_tmp.name});
    img_folder_idx = ones(1,length(img_fnames));
end

end
