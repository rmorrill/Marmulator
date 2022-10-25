function make_params_TEMPLATE()

%%%%% THIS IS INTENDED AS A COMPLETE TEMPLATE - MAKE A COPY OF THIS FILE TO
%%%%% MAKE YOUR OWN PARAMETER GENERATION FUNCTIONS 

save_params_name = 'grid_3_3_brad_6deg_y_offset_90_x_offset_-20_bbox_200'; % e.g. 'grid_3_3_foraging_yj_wrenches_camels_6deg'; 
save_base = 'C:\MATLAB\Marmulator'; % e.g. 'C:\MATLAB\Marmulator'; 

%% SETTINGS:
% general
window_rect = []; % empty to run fullscreen
skip_sync_tests = 1;

calibration_win_len = 200;
calibration_win_ht = 200; 

n_pts_x = 3; % how many different x grid locations
n_pts_y = 3; % how many different y grid locations 
repeats_per_stim = 10; % along w n_pts_x and n_pts_y, determines number of trials: n_pts_x * n_pts_y * repeats_per_stim 
presentation_time = 500; % ms
inter_stim_interval = [1000 2000];
iti_random = false;
order = 'random'; % or 'sequential'
bg_col = 'gray';
stim_rect_size_x = 90; % now applies to moving dots, rects, images
stim_rect_size_y = 90; %


show_curr_bounding_box = true;
show_all_bounding_boxes = true;

bounding_rect_size_x = 200;
bounding_rect_size_y =200; 

manual_bounding_boxes = []; % empty for automoatic OR e.g. [0, 0, 768, 432; 768,0,1536,432; 0,432,768,864; 768,432,1536,864];

% mode
trial_mode = 'foraging'; % 'trial' or 'foraging' - MUST BE SET TO FORAGING FOR RSVP MODE
reward_on = 'location'; % quality or location
time_to_reward = 500; % ms; only applies on foraging expts
time_out_after = 2000; % ms; ; only applies on foraging expts
location_require_quality = true; 

stim_mode = 'images'; % 'images', 'movie', 'moving_dot', 'rectangle', 'spinning', 'smooth pursuit'

if strcmp(stim_mode,'spinning') || strcmp(stim_mode,'smooth pursuit') 
    diam_circle = 500; % diameter of the circular trajectory % px
    speed = 500; %speed of the stimulus along the circular trajectory % px
end

if strcmp(stim_mode,'smooth pursuit')
    angle_endpts = [45, 135, 225, 315]; 
end

% images
img_folder = {fullfile(save_base, 'stim\yj_camels'), fullfile(save_base, 'stim\yj_wrenches'),}; 
img_folder = {fullfile(save_base,'stim\brad')};
pulse_size = false; % true or false

% movies 
movie_folder = 'C:\Users\Hector\Documents\Nature videos\for_calib_stim'; 
movie_rate = 1; % speed of movie playback 

% moving dot parameters
scale_fact_move = 10; % how quickly should dot move?
scale_fact_size = 3; % how quickly should dot change size? 0 = no change
mvdot_sz = 15;
reverse_rate = 0.2; % how often (per frame) should the stim change direction in x, y, or size
reverse_rate_sz = 0.2;

% rsvp - if n_rsvp>1, this will use 'rsvp mode' and trial_mode must be set
% to 'foraging'
n_rsvp = 1; 
rsvp_iti_t = 200; 
break_after = 300; 

% stimulus pre dot
stimulus_pre_dot = true;
stimulus_pre_dot_disappear = true; 
stimulus_pre_time = 500; %ms, does not apply if require_fix_tr_init
stim_pre_dot_sz = 8;

% require fixation for trial start
require_fix_tr_init = true;
fixation_to_init = 400;
time_out_trial_init_s = 10; 

% offset in y
gaze_center_adj_y = 116;
gaze_center_adj_x = 0;
apply_gaze_center_adj = true;

% gaze point drawing
draw_gaze_pt = true;
n_gaze_pts_draw = 1;
dot_col = distinguishable_colors(1);
gaze_pt_sz = 15;
draw_retain_bb_pts = false; % draw and retain the points that land in the bounding box through the session

% gaze point feedback % SOMEWHAT DEPRECATED
color_shift_feedback = false; 
rew_col_start = [127.5,127.5,127.5]; % gray
rew_col_end = [50,205,50]; % lime green

% reward sound feedback 
play_reward_sound = true; 
reward_sound_file = fullfile(save_base, 'trial_outcome_snds\correct.wav'); 

% punishment params 
give_punishments = true; 
punish_length_ms = 500; 
play_punish_sound = true; 
punish_sound_file = fullfile(save_base, 'trial_outcome_snds\incorrect.wav');  


% eyetracking
eye_method = 'gaze_point'; % USE 'gaze_point' 'pupil', 'pupil-glint', 'gaze_point', 'gaze_point_corrected' 'rand' [for no eyetracker input]
n_frames_blink = 8; 

% reward
give_rewards = true;
reward_thresh = 0.95; % 95% of quality codes should be 0 ('glint is good')
reward_on_dur = 0.1; % this shouldn't matter given the pump program just runs on the rising edge
reward_vol = 0.001; % just for logging purposes, again set by pump program and stored in non-vol mem
ard_pin_rew = 10; % which pin has the reward signal?
wait_after_reward = 4;

% wake-up trials 
wake_up_trials = true;
wake_up_every = [4 8];
wake_up_tr_dur = 6;
wake_up_img_fold = fullfile(save_base, 'stim\naturalistic_scenes');
wake_up_stim_size_x = 592;
wake_up_stim_size_y = 444;

%%% DO NOT MODIFY %%%% 
save_params_here = fullfile(save_base, 'experiment_params'); 
save_full = fullfile(save_params_here, save_params_name); 
save(fullfile(save_params_here, save_params_name)); 
fprintf('saved %s \n', save_params_name); 
disp(save_full);  

