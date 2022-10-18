function [eyetrack, calib] = EyeTracker_Calibrate_gui_fcn(reward_pumphand, reward_arduino_pin,...
    subject, expt_params, calib_fname, gaze_offset, repeats_per_stim, time_out_after, ...
    time_to_reward, presentation_time, session_time, eye_method_mouse,...
    require_fix_tr_init, fixation_to_init, time_out_trial_init_s, ...
    reward_today_hand, reward_vol, punish_length_ms, setup_config)

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


PsychJavaSwingCleanup;

%% SETTINGS:
s = load(expt_params);

% unpack everything
save_params_name = s.save_params_name;
save_params_here = s.save_params_here;
eyetracker_toolbox_dir = setup_config.eyetracker_toolbox_dir;
window_rect = s.window_rect;
skip_sync_tests = s.skip_sync_tests;
screen_dist_cm = s.screen_dist_cm;
screen_width_cm = s.screen_width_cm;
screen_ht_cm = s.screen_ht_cm;
screenid_stim = setup_config.screenid_stim;
screenid_ctrl = setup_config.screenid_ctrl;
calibration_win_len = s.calibration_win_len;
calibration_win_ht = s.calibration_win_ht;
n_pts_x = s.n_pts_x;
n_pts_y = s.n_pts_y;
%repeats_per_stim = s.repeats_per_stim;
%presentation_time = s.presentation_time;
inter_stim_interval = s.inter_stim_interval;
iti_random = s.iti_random;
position_order = s.order;
bg_col = s.bg_col;
stim_rect_size_x = s.stim_rect_size_x;
stim_rect_size_y = s.stim_rect_size_y;
show_curr_bounding_box = s.show_curr_bounding_box;
show_all_bounding_boxes = s.show_all_bounding_boxes;
bounding_rect_size_x = s.bounding_rect_size_x;
bounding_rect_size_y = s.bounding_rect_size_y;
manual_bounding_boxes = s.manual_bounding_boxes;
trial_mode = s.trial_mode;
reward_on = s.reward_on;
%time_to_reward = s.time_to_reward;
%time_out_after = s.time_out_after;
location_require_quality = s.location_require_quality;
%location_require_quality = false; % RJM
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
dot_col = s.dot_col;
gaze_pt_sz = s.gaze_pt_sz;
draw_retain_bb_pts = s.draw_retain_bb_pts;
color_shift_feedback = s.color_shift_feedback;
rew_col_start = s.rew_col_start;
rew_col_end = s.rew_col_end;
play_reward_sound = s.play_reward_sound;
reward_sound_file = s.reward_sound_file;
give_punishments = s.give_punishments;
%punish_length_ms = s.punish_length_ms;
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
%reward_vol = s.reward_vol;
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

gaze_center_adj_x = setup_config.default_gaze_center_adjust(1); 
gaze_center_adj_y = setup_config.default_gaze_center_adjust(2); 
apply_gaze_center_adj = true;

%% rsvp setup
% rsvp will use presentation time for each stimulus duration
n_rsvp = 3; 
rsvp_iti_t = 300; 
if n_rsvp>1
    rsvp_mode = true; 
else
    rsvp_mode = false; 
end

image_order = 'random'; 

% if isfield(s, 'apply_gaze_center_adj')
%     gaze_center_adj_y = s.gaze_center_adj_y;
%     gaze_center_adj_x = s.gaze_center_adj_x;
%     apply_gaze_center_adj = s.apply_gaze_center_adj;
% else
%     % offset in y
%     gaze_center_adj_y = 116;
%     gaze_center_adj_x = 0;
%     apply_gaze_center_adj = true;
% end
% gaze_center_adj_y = 210;

if isfield(s, 'wake_up_trials')
    wake_up_trials = s.wake_up_trials;
    wake_up_every = s.wake_up_every;
    wake_up_tr_dur = s.wake_up_tr_dur;
    wake_up_img_fold = s.wake_up_img_fold;
    wake_up_stim_size_x = s.wake_up_stim_size_x;
    wake_up_stim_size_y = s.wake_up_stim_size_y;
else
    wake_up_trials = false;
    wake_up_every = [];
    wake_up_tr_dur = [];
    wake_up_img_fold = '';
    wake_up_stim_size_x = [];
    wake_up_stim_size_y = [];
end

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

%eye_method = 'mouse';
%%
eye_rect = [0 0 1 1] ;
eye = 0; % eye A
rect_col = [255 0 0];

good_pts = [];
good_pts_cols = [];

%session_time = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
%%
PsychJavaSwingCleanup;
%InitializeMatlabOpenGL;
AssertOpenGL;
InitializePsychSound(1);
esc_key = KbName('ESC');
rew_key = KbName('r');

% setup eyetracker-related
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
dot_cols = round(repmat(dot_col', [1,n_gaze_pts_draw])*255);
alphas = linspace(0, 255, n_gaze_pts_draw);
rgba_cols = [dot_cols; alphas];
dotsz = gaze_pt_sz;

% setup for images
if strcmp(stim_mode, 'images') || strcmp(stim_mode,'spinning') || strcmp(stim_mode, 'smooth pursuit')
    if iscell(img_folder)
        img_fnames = []; 
        for i = 1:numel(img_folder)
            img_fnames_tmp = dir(img_folder{i});
            img_fnames_tmp = img_fnames_tmp(~[img_fnames_tmp.isdir]);
            img_fnames_tmp = fullfile(img_folder{i}, {img_fnames_tmp.name}); 
            img_fnames = [img_fnames img_fnames_tmp]; 
        end
    else
        img_fnames_tmp = dir(img_folder);
        img_fnames_tmp = {img_fnames_tmp(~[img_fnames_tmp.isdir]).name};
        img_fnames = fullfile(img_folder, img_fnames_tmp); 
    end
    
    imgs = cell(1,numel(img_fnames));
    for i = 1:numel(img_fnames)
        %image_fname_list{i} = fullfile(img_folder, img_d(i).name);
        [img_tmp, ~,alpha] = imread(img_fnames{i});
        if ~isempty(alpha)
            img_tmp(:,:,4) = alpha;
        end
        imgs{i} = img_tmp;
        im_ht(i) = size(imgs{i}, 1);
        im_wd(i) = size(imgs{i}, 2);
        ar(i) = im_wd(i)/im_ht(i);
    end
    nr_imgs = numel(imgs);
    img_idx = 1;
end


%%
Screen('Preference', 'Verbosity', 0);
Screen('Preference', 'SkipSyncTests', skip_sync_tests)
Screen('Preference', 'VisualDebugLevel', 0)

%screenid_stim = max(Screen('Screens'));
%screenid_ctrl = max(Screen('Screens'))-1;

% get colors
whitecol = WhiteIndex(screenid_stim);
blackcol = BlackIndex(screenid_stim);
graycol = (whitecol + blackcol)/2;

switch bg_col
    case 'gray'
        bg_col_val = graycol;
    case 'white'
        bg_col_val = whitecol;
    case 'black'
        bg_col_val = blackcol;
end


% OPEN WINDOWS
%fprintf('GUI MODE\n\n')
[win, win_rect] = Screen('OpenWindow', screenid_stim, bg_col_val, window_rect);
ctrl_rect_debug = [0 0 1650 1200];
[win_ctrl, win_rect_ctrl] = Screen('OpenWindow', screenid_ctrl, bg_col_val, ctrl_rect_debug);
%[win0, winRect0] = Screen('OpenWindow', screenId, bgcolor * 255, [0 0 300 300], [], [], [], [], [], kPsychGUIWindow);
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('BlendFunction', win_ctrl, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[x_cent, y_cent] = RectCenter(win_rect);


if apply_gaze_center_adj
    y_cent = y_cent + gaze_center_adj_y;
    x_cent = x_cent + gaze_center_adj_x;
end

ctrl_win_text_size = 18;
Screen('TextSize', win_ctrl, ctrl_win_text_size);

screen_hz = Screen('NominalFrameRate', screenid_stim);
ifi=Screen('GetFlipInterval', win);
halfifi = 0.5*ifi;
if rsvp_mode
inter_rsvp_frames = round(rsvp_iti_t/1e3/ifi); 
end

if color_shift_feedback
    n_frames_til_rew = round(time_to_reward/1e3/ifi);
    col_shift_change = nan(3,n_frames_til_rew);
    col_shift_change(1,:) = linspace(rew_col_start(1), rew_col_end(1), n_frames_til_rew);
    col_shift_change(2,:) = linspace(rew_col_start(2), rew_col_end(2), n_frames_til_rew);
    col_shift_change(3,:) = linspace(rew_col_start(3), rew_col_end(3), n_frames_til_rew);
    col_shift_change(4,:) = ones(1,n_frames_til_rew)*0.5;
end


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


if strcmp(stim_mode, 'movie')
    moviefiles_all = dir(movie_folder);
    movienames =  {moviefiles_all(~[moviefiles_all.isdir]).name};
    %moviePtr_ctrl = Screen('OpenMovie', win_ctrl, fullfile(movie_folder, movienames{1}));
    moviePtr = Screen('OpenMovie', win, fullfile(movie_folder, movienames{1}), [], 1, 64);
    Screen('PlayMovie', moviePtr, movie_rate, 1);
end

% set up bounding boxes for display on the control screen
% first off, give them unique colors
cols = round(distinguishable_colors(n_calib_pts+1)*255);
bounding_rect_tmp = [0,0,bounding_rect_size_x, bounding_rect_size_y];
if ~isempty(manual_bounding_boxes)
    if apply_gaze_center_adj
        manual_bounding_boxes = OffsetRect(manual_bounding_boxes, gaze_center_adj_x, gaze_center_adj_y);
    end
    
    bounding_rects = manual_bounding_boxes;
else
    bounding_rects = CenterRectOnPoint(bounding_rect_tmp, all_pts(:,1), all_pts(:,2));
end

n_trs_tot = n_calib_pts*repeats_per_stim;


tr_seq = []; % trial sequence - note that this determines the position of stimuli only. 
% img_seq determines the sequence of images
if strcmp(position_order, 'random')
    for q = 1:repeats_per_stim
        tr_seq = [tr_seq randperm(n_calib_pts)];
    end
else
    tr_seq = repmat(1:n_calib_pts, 1, repeats_per_stim);
end

% set up img sequence 
n_rsvps = n_trs_tot*n_rsvp; 
x = floor(270/nr_imgs); 
img_seq = [repmat(1:nr_imgs, 1, x) 1:mod(n_rsvps, x*nr_imgs)]; 
if strcmp(image_order, 'random')
    img_seq =  img_seq(randperm(numel(img_seq))); 
end
img_seq = reshape(img_seq, [n_rsvp, n_trs_tot]); 


if wake_up_trials
    y = nan(1,n_trs_tot); 
    z = [tr_seq; y];
    n_wake_up_trs = ceil(n_trs_tot/mean(wake_up_every))+1;
    wue = wake_up_every -1;
    insert_every = round(rand(1,n_wake_up_trs)*diff(wue)+wue(1));
    insert_idx = cumsum(insert_every); 
    insert_idx = insert_idx(insert_idx<=n_trs_tot);
    z(2,insert_idx) = 0;
    z = z(~isnan(z))';
    tr_seq = z;
    n_trs_tot = numel(tr_seq);
    n_wake_up_trs = sum(tr_seq == 0);
    wake_up_stim_frames = round(wake_up_tr_dur/ifi);
    % read images
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
    nr_wu_imgs = numel(wake_up_imgs);
    wu_img_idx = 1;
else
    n_wake_up_trs = 0; 
end


if strcmp(stim_mode,'spinning') || strcmp(stim_mode,'smooth pursuit')
    n_trs_tot = repeats_per_stim;
end

image_displayed = cell(n_rsvp, n_trs_tot); % will be empty for all except 'images' sessions
wake_up_image_displayed = cell(1, n_trs_tot); 

% calculate frames per stim, iti
stim_frames = round(presentation_time/1e3/ifi);
pulse_rate = (rand(1, n_trs_tot)+0.5)*3;
if strcmp(trial_mode, 'trial')
    t_sin = 0:ifi:(presentation_time/1e3);
elseif strcmp(trial_mode, 'foraging')
    t_sin = 0:ifi:(time_out_after+time_to_reward/1e3);
end

% set up audio feedback
if play_reward_sound
    [aud_y, aud_fs] = audioread(reward_sound_file);
end

if play_punish_sound
    [aud_pun_y, aud_pun_fs] = audioread(punish_sound_file);
end

%iti_frames = round(inter_stim_interval(1)/1e3/ifi);
if iti_random
    iti_frames = round((rand(1,n_trs_tot)*diff(inter_stim_interval)+inter_stim_interval(1))/1e3/ifi);
else
    iti_frames = ones(1,n_trs_tot)*round(inter_stim_interval(1)/1e3/ifi);
end

stimulus_pre_frames = round(stimulus_pre_time/1e3/ifi);

Screen('FillRect', win, bg_col_val)
%tex1 = Screen('MakeTexture', win, img);
seqidx = 0;
vbl = Screen('Flip', win);
t_start_sec = GetSecs();


reward_trial = false(1,n_trs_tot); % for recording whether trial was rewarded
punish_trial = false(1,n_trs_tot);
reward_time = nan(1,n_trs_tot); % for recording the time of a reward
time_to_bb = nan(1,n_trs_tot);
trial_aborted = false(1, n_trs_tot);
manual_reward_t = nan(1,n_trs_tot);
stim_pre_start = nan(1,n_trs_tot);
stim_pre_end = nan(1,n_trs_tot);
trial_init_timed_out = false(1,n_trs_tot);
calib_st_t = nan(1,n_trs_tot); 
calib_end_t = nan(1,n_trs_tot); 


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

reward_ct = 0;
man_reward_ct = 0;
exit_flag = 0;
mov_dur = [];
fix_pre_fr_ctr = 0;

% start the stimulus loop
for i = 1:n_trs_tot
    
    curr_tr = i;
    seqidx = tr_seq(i);
    stfridx = 0;
    rsvpfridx = 0; 
    rsvp_ctr = 0; 
    loop_brk_ctr = 0;
    pre_stim_timer = 0;
    manual_reward_flag = false;
    entered_bb = false; % did eye go into bb?
    
    Screen('FillRect', win, bg_col_val)
    vbl = Screen('Flip', win);
    
    Screen('FillRect', win_ctrl, bg_col_val)
    vbl2 = Screen('Flip', win_ctrl);
    
    if seqidx ~= 0
        xcurr = all_pts(seqidx,1);
        ycurr = all_pts(seqidx,2);
        stim_rect = [xcurr-stim_rect_size_x/2 ycurr-stim_rect_size_y/2 xcurr+stim_rect_size_x/2 ycurr+stim_rect_size_y/2];
        curr_bb = bounding_rects(seqidx,:);
    else
        % put stim rect around center
        xcurr = x_cent;
        ycurr = y_cent;
        stim_rect = CenterRectOnPoint([0, 0, wake_up_stim_size_x, wake_up_stim_size_y], xcurr, ycurr);
        curr_bb = stim_rect;
    end
    
    drawInfoText();
    drawBoundingBoxes();
    if draw_retain_bb_pts
        drawGoodEyePts(0);
    end
    
    for j = 1:iti_frames(i)
        idx_all = idx_all +1;
        drawInfoText();
        drawBoundingBoxes();
        if draw_retain_bb_pts
            drawGoodEyePts(0);
        end
        vbl = Screen('Flip', win, vbl + halfifi);
        vbl2 = Screen('Flip', win_ctrl, vbl2 + halfifi);
        get_eyetracker_draw_dots();
        
    end

    end_stim_pre = 0;
    j = 0;
    fix_pre_fr_ctr = 0;
    blink_ctr = 0;
    % ENTER STIMULUS PRE
    if stimulus_pre_dot
        stps = GetSecs();
        stim_pre_start(i) = stps - t_start_sec;
        while ~end_stim_pre
            
            j = j+1;
            %for j = 1:stimulus_pre_frames
            idx_all = idx_all +1;
            drawInfoText();
            drawBoundingBoxes();
            if seqidx ~= 0
                Screen('DrawDots', win_ctrl, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                Screen('DrawDots', win, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
            else
                Screen('DrawDots', win_ctrl, [xcurr, ycurr], stim_pre_dot_sz, whitecol, [], 1);
                Screen('DrawDots', win, [xcurr, ycurr], stim_pre_dot_sz, whitecol, [], 1);
            end
            
            vbl = Screen('Flip', win, vbl + halfifi);
            vbl2 = Screen('Flip', win_ctrl, vbl2 + halfifi);
            
            [eyeposx_cur, eyeposy_cur] = get_eyetracker_draw_dots();
            curr_in_bb = IsInRect(eyeposx_cur, eyeposy_cur, curr_bb);
            
            pre_stim_timer = GetSecs() - stps;
            [~,~,kCode] = KbCheck(0);
            if find(kCode) == esc_key
                disp('ESC key recognized, exiting');
                exit_flag = 1;
                % calib_end_t(i) =  GetSecs()-t_start_sec;
                trial_aborted(i) = true;
                break
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
                elseif (~curr_in_bb || ~qual_check) && fix_pre_fr_ctr>3 && blink_ctr<n_frames_blink
                    blink_ctr = blink_ctr + 1;
                else
                    blink_ctr = 0;
                    fix_pre_fr_ctr = 0;
                end
                
                %sca
                %keyboard
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
    
    end_stim = 0;
    blink_ctr = 0;
    rsvp_ctr = 1; 
    rsvp_break_ctr = 0; 
    inter_rsvp_fr_ctr = 1; 
    
    % ENTER STIMULUS 
    while ~end_stim && ~exit_flag % loops for every frame 
        
        img_idx = img_seq(rsvp_ctr, i); 
        
        if rsvp_mode && rsvp_ctr > 1 && inter_rsvp_fr_ctr < inter_rsvp_frames
            % go into a break
            %sca
            %keyboard
            inter_rsvp = true;
        else
            inter_rsvp = false;
            rsvpfridx = rsvpfridx + 1;
        end
        %rsvpfridx
        
        stfridx = stfridx + 1; % stim frame idx, counts every frame from stim start (rsvp_1) through trial stim end (rsvp_n)
        idx_all = idx_all +1; % idx for all frames recorded
        
        if ~trial_init_timed_out(i)
            % draw gray background:
            Screen('FillRect', win, bg_col_val);
            Screen('FillRect', win_ctrl, bg_col_val);
            
            % draw informational text into control window:
            drawInfoText();
            drawBoundingBoxes();
            
            if seqidx ~= 0 && strcmp(trial_mode, 'foraging') && color_shift_feedback && loop_brk_ctr>0
                Screen('FillRect', win, col_shift_change(:,loop_brk_ctr), curr_bb);
                Screen('FillRect', win_ctrl, col_shift_change(:,loop_brk_ctr), curr_bb);
            end
            
            [eyeposx_cur, eyeposy_cur] = get_eyetracker_draw_dots();
            eye_data_qual_curr(stfridx) = eyetracker_qual(idx_all);
            
            %curr_in_bb = IsInRect(eyepos_x(idx_all), eyepos_y(idx_all), curr_bb);
            curr_in_bb = IsInRect(eyeposx_cur, eyeposy_cur, curr_bb);
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
                        tex1 = Screen('MakeTexture', win, imgs{img_idx});
                        Screen('DrawTexture', win, tex1, [], stim_rect_curr)
                        tex2 = Screen('MakeTexture', win_ctrl, imgs{img_idx});
                        Screen('DrawTexture', win_ctrl, tex2, [], stim_rect_curr)
                        
                    case 'smooth pursuit'
                        xcurr = all_pts{seqidx}(stfridx,1);
                        ycurr = all_pts{seqidx}(stfridx,2);
                        stim_rect_curr = [xcurr-stim_rect_size_x/2 ycurr-stim_rect_size_y/2 xcurr+stim_rect_size_x/2 ycurr+stim_rect_size_y/2];
                        tex1 = Screen('MakeTexture', win, imgs{img_idx});
                        Screen('DrawTexture', win, tex1, [], stim_rect_curr)
                        tex2 = Screen('MakeTexture', win_ctrl, imgs{img_idx});
                        Screen('DrawTexture', win_ctrl, tex2, [], stim_rect_curr)
                        %                 if stimulus_pre_dot
                        %                     Screen('DrawDots', win_ctrl, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                        %                     Screen('DrawDots', win, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                        %                 end
                        
                    case 'images'
                        tex1 = Screen('MakeTexture', win, imgs{img_idx});
                        
                        if isempty(image_displayed{i})
                            image_displayed{rsvp_ctr, i} = img_fnames{img_idx};
                        end
                        
                        if pulse_size
                            pulse_sin = (cos(2*pi*pulse_rate(i)*t_sin) +1)/2;
                            pulse_sin = (pulse_sin + 1)/2; % 50% modulation depth
                            rsx_curr = stim_rect_size_x*pulse_sin(stfridx);
                            rsy_curr = round(rsx_curr/ar(img_idx));
                        else
                            rsx_curr = stim_rect_size_x;
                            rsy_curr = round(rsx_curr/ar(img_idx));
                        end
                        
                        stim_rect_curr = [xcurr-rsx_curr/2 ycurr-rsy_curr/2 xcurr+rsx_curr/2 ycurr+rsy_curr/2];
                        Screen('DrawTexture', win, tex1, [], stim_rect_curr)
                        tex2 = Screen('MakeTexture', win_ctrl, imgs{img_idx});
                        Screen('DrawTexture', win_ctrl, tex2, [], stim_rect_curr)
                        
                        Screen('DrawDots', win_ctrl, [eyeposx_cur, eyeposy_cur],...
                            dotsz, rgba_cols(:,end), [], 1);
                        
                        if stimulus_pre_dot && ~stimulus_pre_dot_disappear
                            Screen('DrawDots', win_ctrl, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                            Screen('DrawDots', win, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
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
                        Screen('DrawDots', win_ctrl, mvdot_coords,...
                            size_curr, [0 0 0], [], 1);
                        
                        curr_x = mvdot_coords(1);
                        curr_y = mvdot_coords(2);
                        
                        if stimulus_pre_dot
                            Screen('DrawDots', win_ctrl, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                            Screen('DrawDots', win, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                        end
                        mov_dur(end+1) = GetSecs()-t_mov_start;
                    case 'movie'
                        t_mov_start = GetSecs();
                        waitforimage = 0;
                        movtex_new = Screen('GetMovieImage', win_ctrl, moviePtr, waitforimage);
                        if movtex_new>0
                            movtex = movtex_new;
                        end
                        Screen('DrawTexture', win_ctrl, movtex, [],  stim_rect)
                        Screen('DrawTexture', win, movtex, [],  stim_rect)
                        % for movies, re-draw dots
                        Screen('DrawDots', win_ctrl, [eyeposx_cur, eyeposy_cur],...
                            dotsz, rgba_cols(:,end), [], 1);
                        if stimulus_pre_dot
                            Screen('DrawDots', win_ctrl, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                            Screen('DrawDots', win, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                        end
                        mov_dur(end+1) = GetSecs()-t_mov_start;
                    case 'rectangle'
                        Screen('FillRect', win, rect_col, stim_rect);
                        %Screen('DrawDots', win, all_pts(seqidx,:), 5, blackcol, [], 1);
                        Screen('FillRect', win_ctrl, rect_col, stim_rect);
                        %Screen('DrawDots', win_ctrl, all_pts(seqidx,:), 5, blackcol, [], 1)
                        if stimulus_pre_dot
                            Screen('DrawDots', win_ctrl, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                            Screen('DrawDots', win, all_pts(seqidx,:), stim_pre_dot_sz, whitecol, [], 1);
                        end
                end
            elseif seqidx == 0 % draw wakeup trial
                
                if isempty(wake_up_image_displayed{i})
                    wake_up_image_displayed{i} = wu_image_fname_list{wu_img_idx};
                end
                
                % stim_rect_curr = [xcurr-rsx_curr/2 ycurr-rsy_curr/2 xcurr+rsx_curr/2 ycurr+rsy_curr/2];
                tex1 = Screen('MakeTexture', win_ctrl, wake_up_imgs{wu_img_idx});
                Screen('DrawTexture', win, tex1, [], stim_rect)
                tex2 = Screen('MakeTexture', win_ctrl, wake_up_imgs{wu_img_idx});
                Screen('DrawTexture', win_ctrl, tex2, [], stim_rect)
                
                Screen('DrawDots', win_ctrl, [eyeposx_cur, eyeposy_cur],...
                    dotsz, rgba_cols(:,end), [], 1);
                
            elseif inter_rsvp
                %sca
                %keyboard
                Screen('DrawDots', win_ctrl, [eyeposx_cur, eyeposy_cur],...
                    dotsz, rgba_cols(:,end), [], 1);
                inter_rsvp_fr_ctr = inter_rsvp_fr_ctr + 1; 
                %inter_rsvp_fr_ctr 
            end
            
            
            if strcmp(reward_on, 'location')
                %disp(loop_brk_ctr)
                if location_require_quality
                    qual_check = eyetracker_qual(idx_all)<2;
                else
                    qual_check = true;
                end
                
                if rsvp_mode
                    if curr_in_bb && qual_check
                        rsvp_break_ctr = 0; 
                    elseif (~curr_in_bb || ~qual_check)
                        rsvp_break_ctr = rsvp_break_ctr + 1; 
                        %rsvp_break_ctr
                    end
                else
                    if curr_in_bb && qual_check
                        entered_bb = true;
                        loop_brk_ctr = loop_brk_ctr + 1;
                        blink_ctr = 0;
                        if isnan(time_to_bb(i))
                            time_to_bb(i) = GetSecs() - t_start_sec;
                        end
                    elseif (~curr_in_bb || ~qual_check) && loop_brk_ctr>3 && blink_ctr<n_frames_blink
                        blink_ctr = blink_ctr + 1;
                    else
                        blink_ctr = 0;
                        loop_brk_ctr = 0;
                    end
                end
                
            elseif strcmp(reward_on, 'quality')
                if eyetracker_qual(idx_all)<2
                    loop_brk_ctr = loop_brk_ctr + 1;
                else
                    loop_brk_ctr = 0;
                end
            end
            
            vbl = Screen('Flip', win, vbl + halfifi);
            vbl2 = Screen('Flip', win_ctrl, vbl2 + halfifi);
            
            %if strcmp(stim_mode, 'movie')
            %    Screen('Close', movtex);
            %end
            if (strcmp(stim_mode, 'images') || seqidx == 0) && ~inter_rsvp
                Screen('Close', tex1);
                Screen('Close', tex2);
            end
        end
        
        if stfridx == 1
            calib_st_t(i) = GetSecs()-t_start_sec;
        end
        
        if rsvp_mode && rsvpfridx >= stim_frames
            rsvp_ctr = rsvp_ctr + 1;
            inter_rsvp_fr_ctr = 1;
            rsvpfridx = 1;
            % change the image
%             img_idx = img_idx + 1;
%             if img_idx > nr_imgs
%                 img_idx = 1;
%             end
        end
        
        % check if we need to end
        if strcmp(trial_mode, 'trial') || seqidx == 0
            if  seqidx ~= 0 && stfridx >= stim_frames && rsvp_ctr > n_rsvp % finished successfully
                end_stim = 1;
            elseif seqidx == 0 && stfridx >= wake_up_stim_frames
                end_stim = 1;
            end
        elseif strcmp(trial_mode, 'foraging')
            if rsvp_mode
                if rsvp_ctr > n_rsvp
                    end_stim = 1;
                elseif rsvp_break_ctr >= round(time_out_after/1e3/ifi)
                    end_stim = 1; 
                    fprintf('TRIAL BREAK: FIXATION BROKEN, %0.1f s \n', GetSecs() - calib_st_t(i) - t_start_sec);
                end
            else
                if loop_brk_ctr >= round(time_to_reward/1e3/ifi)
                    end_stim = 1;
                elseif ~(curr_in_bb && qual_check) && stfridx >= round(time_out_after/1e3/ifi) && ~entered_bb
                    end_stim = 1;
                    fprintf('TRIAL BREAK: TIMED OUT, %0.1f s \n', GetSecs() - calib_st_t(i) - t_start_sec);
                elseif ~(curr_in_bb && qual_check) && stfridx >= round(time_out_after/1e3/ifi) && entered_bb
                    end_stim = 1;
                    fprintf('TRIAL BREAK: FIXATION BROKEN, %0.1f s \n', GetSecs() - calib_st_t(i) - t_start_sec);
                end
            end
        elseif trial_init_timed_out(i)
            end_stim = 1;
        end
        
        [~,~,kCode] = KbCheck(0);
        if find(kCode) == esc_key
            disp('ESC key recognized, exiting');
            exit_flag = 1;
            calib_end_t(i) =  GetSecs()-t_start_sec;
            trial_aborted(i) = true;
            break
        elseif find(kCode) == rew_key
            manual_reward_t(i) =  GetSecs()-t_start_sec;
            man_reward_ct = man_reward_ct + 1;
            fprintf('manual reward! time = %0.3fs\n', manual_reward_t(i));
            end_stim = 1;
            manual_reward_flag = true;
        end
        
        
        % if ending, note the time:
        if end_stim
            calib_end_t(i) = GetSecs()-t_start_sec;
        end
        
        %Screen('FillRect', win, bg_col_val);
        %Screen('FillRect', win_ctrl, bg_col_val);
        %vbl = Screen('Flip', win, vbl + halfifi);
        %vbl2 = Screen('Flip', win_ctrl, vbl + halfifi);
        
        if end_stim
           % sca
           % keyboard
            Screen('FillRect', win, bg_col_val);
            Screen('FillRect', win_ctrl, bg_col_val);
            
            if seqidx ~= 0 && give_rewards && color_shift_feedback && loop_brk_ctr > 0
                Screen('FillRect', win, col_shift_change(:,loop_brk_ctr), curr_bb);
                Screen('FillRect', win_ctrl, col_shift_change(:,loop_brk_ctr), curr_bb);
            end
            drawInfoText();
            drawBoundingBoxes();
            if draw_retain_bb_pts
                drawGoodEyePts(0);
            end
            vbl = Screen('Flip', win, vbl + halfifi);
            vbl2 = Screen('Flip', win_ctrl, vbl2 + halfifi);
            
            
            if give_rewards
                if strcmp(trial_mode, 'trial') || seqidx == 0 
                    if strcmp(reward_on, 'quality') || seqidx == 0
                        frac_good = sum(eye_data_qual_curr < 2)/numel(eye_data_qual_curr);
                    elseif strcmp(reward_on, 'location')
                        frac_good = sum(eye_in_bb_curr)/numel(eye_in_bb_curr);
                    end
                    
                    if frac_good > reward_thresh
                        reward_this_trial = true;
                    else
                        reward_this_trial = false;
                    end
                    
                elseif strcmp(trial_mode, 'foraging')
                    if rsvp_mode
                        if rsvp_break_ctr >= round(time_out_after/1e3/ifi)
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
                    reward_trial(i) = true;
                    if ~manual_reward_flag
                        reward_ct = reward_ct + 1;
                    end
                    % fprintf('****REWARD TRIAL %d, %0.2f good\n', i, frac_good);
                    fprintf('****REWARD TRIAL %d\n', i);
                    if play_reward_sound
                        sound(aud_y, aud_fs);
                        %PsychPortAudio('Start', pa, 1);
                        %play(aud_playr);
                        %play(aud_playr);
                        disp('sound played');
                        %get(aud_playr)
                    end
                    if ~isempty(reward_pumphand)
                        reward_time(i) = GetSecs()-t_start_sec;
                        if reward_serial
                            writeline(reward_pumphand, 'RUN');
                            WaitSecs(reward_on_dur);
                        else
                            reward_pumphand.digitalWrite(reward_arduino_pin, 1);
                            WaitSecs(reward_on_dur);
                            reward_pumphand.digitalWrite(10, 0);
                        end
                        set(reward_today_hand, 'String', sprintf('%0.3f mL', (reward_ct + man_reward_ct)*reward_vol + start_reward_vol));
                        drawnow;
                        WaitSecs(wait_after_reward);
                        
                    end
                    
                    %WaitSecs(1);
                    
                else
                    fprintf('no reward: trial %d\n', i);
                end
            end
            
            if give_punishments && (give_rewards && ~reward_this_trial)
                punish_trial(i) = true;
                %fprintf('Punish on\n');
                if play_punish_sound
                    sound(aud_pun_y, aud_pun_fs);
                end
                for p = 1:round(punish_length_ms/1e3/ifi)
                    Screen('FillRect', win, blackcol);
                    Screen('FillRect', win_ctrl, blackcol);
                    drawInfoText(1);
                    drawBoundingBoxes();
                    if draw_retain_bb_pts
                        drawGoodEyePts(0);
                    end
                    vbl = Screen('Flip', win, vbl + halfifi);
                    vbl2 = Screen('Flip', win_ctrl, vbl + halfifi);
                end
                %fprintf('Punish off\n');
            end
        end
    end
    
%     if strcmp(stim_mode, 'images') && ~rsvp_mode
%         img_idx = img_idx + 1;
%         if img_idx > nr_imgs
%             img_idx = 1;
%         end
%     end
    
    if seqidx == 0
        wu_img_idx = wu_img_idx +1;
        if wu_img_idx > nr_wu_imgs
            wu_img_idx = 1;
        end
    end
    
    if exit_flag
        break
    end
    
end


%sca
%keyboard

Screen('FillRect', win, bg_col_val);
Screen('FillRect', win_ctrl, bg_col_val);
for j = 1:iti_frames*2
    vbl = Screen('Flip', win);
    vbl2 = Screen('Flip', win_ctrl);
end

t2 = GetSecs;
Screen('CloseAll');
sca

%% gather data for save
calib.n_pts_x = n_pts_x;
calib.n_pts_y = n_pts_y;
calib.repeats_per_stim = repeats_per_stim;
if any(trial_aborted)
    calib.n_completed = find(trial_aborted) - 1;
else
    calib.n_completed = i;
end
calib.pts = all_pts;
calib_st_t(isnan(calib_st_t)) = []; 
calib_end_t(isnan(calib_end_t)) = []; 
calib.start_t = calib_st_t;
calib.end_t = calib_end_t;
calib.time_to_bounding_box = time_to_bb;
calib.sequence = tr_seq;
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

calib_settings.disp_rect = win_rect;
calib_settings.presentation_time = presentation_time;
calib_settings.inter_stim_interval = inter_stim_interval;
calib_settings.iti_random = iti_random;
calib_settings.iti_frames = iti_frames;
calib_settings.position_order = position_order;
calib_settings.bg_col = bg_col;
calib_settings.stim_rect_size_x = stim_rect_size_x;
calib_settings.stim_rect_size_y = stim_rect_size_y;
calib_settings.stim_mode = stim_mode;
calib_settings.img_folder = img_folder;
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
calib_settings.draw_gaze_pt = draw_gaze_pt;
calib_settings.n_gaze_pts_draw = n_gaze_pts_draw;
calib_settings.gaze_pt_dot_col = dot_col;
calib_settings.gaze_pt_dot_sz = gaze_pt_sz;
calib_settings.stimulus_pre_dot = stimulus_pre_dot;
calib_settings.stimulus_pre_time = stimulus_pre_time;
calib_settings.stim_pre_dot_sz = stim_pre_dot_sz;
calib_settings.stimulus_pre_dot_disappear = stimulus_pre_dot_disappear; 

calib_settings.expt_params = expt_params;
calib_settings.wake_up_trials = wake_up_trials; 
calib_settings.wake_up_every = wake_up_every; 
calib_settings.wake_up_img_dir = wake_up_img_fold; 
calib_settings.wake_up_stim_size_x = wake_up_stim_size_x; 
calib_settings.wake_up_stim_size_x = wake_up_stim_size_y; 
calib_settings.wake_up_tr_dur = wake_up_tr_dur; 

% offset in y
calib_settings.gaze_center_adj_y = gaze_center_adj_y;
calib_settings.gaze_center_adj_x = gaze_center_adj_x;
calib_settings.apply_gaze_center_adj = apply_gaze_center_adj;

calib_settings.require_fix_tr_init = require_fix_tr_init;
calib_settings.fixation_to_init = fixation_to_init;
calib_settings.time_out_trial_init_s = time_out_trial_init_s;

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
settings.hostname = hname;
settings.time_save = datestr(now, 'yyyy-dd-mm_HH-MM-SS');
settings.time_start = session_time;
settings.run_time = GetSecs() - t_start_sec;
settings.ifi_monitor = ifi;

reward.give_rewards = true;
reward.nr_rewards_given = sum(reward_trial);
reward.reward_sequence = reward_trial;
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

punish.give_punishments = give_punishments;
punish.nr_puns_given = sum(punish_trial);
punish.pun_sequence = punish_trial;
punish.length_ms = punish_length_ms;
punish.play_sound = play_punish_sound;
punish.sound_file = punish_sound_file;


%% save data
savefname = sprintf('calib_%s_%s.mat', subject, session_time);
save(fullfile(save_data_dir, savefname), 'calib', 'calib_settings', 'eyetrack', 'settings', 'reward', 'punish');
fprintf('session data saved to %s\n', fullfile(save_data_dir, savefname));

% try to save to remote 
try 
    save(fullfile(save_data_dir_remote, savefname), 'calib', 'calib_settings', 'eyetrack', 'settings', 'reward', 'punish');
    fprintf('remote save: session data saved to %s\n', fullfile(save_data_dir_remote, savefname));
catch me
    disp(me); 
    disp('failed to save session data to remote'); 
end


%% execute plots automatically

if contains(calib_settings.expt_params,'center_point')
    get_mean_x_y_pts(fullfile(save_data_dir, savefname))
end

%% add to subject log table automatically


%%% NESTED FUNCTIONS FOR DRAWING ON SCREEN
    function [eyeposx_cur, eyeposy_cur, eye_data_qual] = get_eyetracker_draw_dots()
        %[eyepos_x(idx_all), eyepos_y(idx_all)] = vpx_GetGazePointCorrected(eye);
        %[eyepos_x_tmp, eyepos_y_tmp] = vpx_GetPupilPoint(eye);
        
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
            if ~isempty(cProj)
                cProj
                raw_eye = [eyepos_x_tmp,eyepos_y_tmp];
                T = [cX,cY,cProj];
                transformed_eye = homography_transform(raw_eye,T); 
                eyepos_x(idx_all) = transformed_eye(:,1);
                eyepos_y(idx_all) = transformed_eye(:,2);
            else
                if length(cX) == 3
                    eyepos_x(idx_all) = eyepos_x_tmp*cX(2) + eyepos_y_tmp*cX(3) + cX(1);
                elseif length(cX) == 2
                    eyepos_x(idx_all) = eyepos_x_tmp*cX(2) + cX(1);
                end
                
                if length(cY) == 3
                    eyepos_y(idx_all) = eyepos_y_tmp*cY(2) + eyepos_x_tmp*cY(3) + cY(1);
                elseif length(cY) == 2
                    eyepos_y(idx_all) = eyepos_y_tmp*cY(2) + cY(1);
                end
            end
            
            %eyepos_x(idx_all) = eyepos_x_tmp*cX(1) + eyepos_y_tmp*cX(2) + cX(3);
            %eyepos_y(idx_all) = eyepos_y_tmp*cY(1) + eyepos_x_tmp*cY(2) + cY(3);
            % fprintf('eyepos raw: %0.2f, %0.2f; pix: %0.1f, %0.1f\n', eyepos_x_tmp, eyepos_y_tmp,...
            %    eyepos_x(idx_all),  eyepos_y(idx_all));
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
        
        if draw_gaze_pt
            Screen('DrawDots', win_ctrl, draw_coords,...
                dotsz, rgba_cols(:,(end-min(idx_all, n_gaze_pts_draw)+1):end), [], 1);
            % fprintf('%d, %d\n', eyepos_x_tmp, eyepos_y_tmp);
        end
        eyeposx_cur = draw_coords(1,end);
        eyeposy_cur = draw_coords(2,end);
        %fprintf('eyepos pix: %0.1f %0.1f\n',  eyeposx_cur, eyeposy_cur);
    end

    function drawInfoText(duringPun)
        if nargin == 0
            duringPun = 0;
            textcol = [0 0 0];
        end
        if duringPun
            textcol = [255 255 255];
        else
            textcol = [0 0 0];
        end
        
        Screen('DrawText', win_ctrl, sprintf('mode: %s, %s; tr %d of %d, tr time: %dms', stim_mode, trial_mode, curr_tr, n_trs_tot, round(stfridx*ifi*1e3)), 80, win_rect(4)+30, textcol);
        Screen('DrawText', win_ctrl, sprintf('%d auto reward, %d manual (reward on %s)', reward_ct, man_reward_ct, reward_on), 80, win_rect(4)+55, textcol);
        Screen('DrawText', win_ctrl, sprintf('eyetracking: %s', eye_method), 80, win_rect(4)+80, textcol);
        Screen('FrameRect', win_ctrl, [0 0 0], win_rect, 1);
        Screen('DrawText', win_ctrl, sprintf('display screen'), mean([win_rect(1), win_rect(3)])-70, win_rect(4), textcol);
        
        if strcmp(trial_mode, 'foraging')
            Screen('DrawText', win_ctrl, sprintf('foraging: time out %dms, reward %dms', time_out_after, time_to_reward), ...
                80, win_rect(4)+105, textcol);
            Screen('DrawText', win_ctrl, sprintf('correct time: %dms', round(loop_brk_ctr*ifi *1e3)), ...
                80, win_rect(4) + 130, textcol);
        end
        
        Screen('DrawText', win_ctrl, sprintf('Calibration loaded: %d', apply_calib), 80, win_rect(4)+160, textcol);
        if length(cX) == 3
            Screen('DrawText', win_ctrl, sprintf('calib x: %0.1f, %0.1f, %0.1f', cX(1), cX(2), cX(3)), 80, win_rect(4)+180, textcol);
        elseif length(cX) == 2
            Screen('DrawText', win_ctrl, sprintf('calib x: %0.1f, %0.1f', cX(1), cX(2)), 80, win_rect(4)+180, textcol);
        end
        
        if length(cY) ==3
            Screen('DrawText', win_ctrl, sprintf('calib y: %0.1f, %0.1f, %0.1f', cY(1), cY(2), cY(3)), 80, win_rect(4)+200, textcol);
        elseif length(cY) == 2
            Screen('DrawText', win_ctrl, sprintf('calib y: %0.1f, %0.1f', cY(1), cY(2)), 80, win_rect(4)+200, textcol);
        end
        
        if ~isempty(cProj)
            Screen('DrawText', win_ctrl, sprintf('calib proj: %0.1f, %0.1f', cProj(1), cProj(2)), 80, win_rect(4)+220, textcol);
        end
        
        Screen('DrawText', win_ctrl, sprintf('gaze offset x: %0.3f', gaze_offset_x), 580, win_rect(4)+180, textcol);
        Screen('DrawText', win_ctrl, sprintf('gaze offset y: %0.3f', gaze_offset_y), 580, win_rect(4)+200, textcol);
        
        if give_punishments
            Screen('DrawText', win_ctrl, sprintf('give time-out punish: %d, dur: %dms, nr punish: %d', give_punishments, punish_length_ms,...
                sum(punish_trial)), 800, win_rect(4)+30, textcol);
        else
            Screen('DrawText', win_ctrl, sprintf('time-out punish: %d', give_punishments), 800, win_rect(4)+30, textcol);
        end
        
        Screen('DrawText', win_ctrl, sprintf('Fixation required for tr start: %d', require_fix_tr_init),  800, win_rect(4)+55, textcol);
        Screen('DrawText', win_ctrl, sprintf('Pre-stim time in box: %dms\n', round(fix_pre_fr_ctr*ifi*1e3)), 800, win_rect(4)+80, textcol);
        Screen('DrawText', win_ctrl, sprintf('Time pre-stim total: %dms\n', round(pre_stim_timer*1e3)), 800, win_rect(4)+105, textcol);
        
        
        
        Screen('DrawText', win_ctrl, sprintf('Subject: %s, Session: %s', subject, session_time), 740, win_rect(4)+280);
    end

    function drawBoundingBoxes()
        if show_curr_bounding_box
            if seqidx ~= 0
                Screen('FrameRect', win_ctrl, cols(seqidx+1,:), curr_bb, 5);
            else
                Screen('FrameRect', win_ctrl, whitecol, curr_bb, 5);
            end
        end
        if show_all_bounding_boxes && seqidx ~= 0
            Screen('FrameRect', win_ctrl, cols(2:end,:)', bounding_rects', 1);
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
        if ~isempty(good_pts)
            Screen('DrawDots', win_ctrl, good_pts, dotsz-5, good_pts_cols, [], 1);
        end
    end


end


function draw_coords = ret_dots_in_rect(x,y,disp_rect,eye_rect)
%sca
%keyboard
disp_ht = disp_rect(4)-disp_rect(2);
disp_len = disp_rect(3)-disp_rect(1);

eyewin_xoffset = eye_rect(1) * disp_len;
eyewin_yoffset = eye_rect(2) * disp_ht;
eyewin_ht = (eye_rect(4) - eye_rect(2))*disp_ht;
eyewin_len = (eye_rect(3) - eye_rect(1))*disp_len;

% scale x
%draw_coords(1,:) = round(((1-x)*eyewin_len) + eyewin_xoffset);
draw_coords(1,:) = round((x*eyewin_len) + eyewin_xoffset);

% scale y
draw_coords(2,:) = round(((y)*eyewin_ht) + eyewin_yoffset);
end

