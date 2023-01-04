function logData_bysession(log_dir,session_file_path)

load(session_file_path); 

% check if subject folder exists
if ~exist(fullfile(log_dir,settings.subject))
       mkdir(fullfile(log_dir,settings.subject))
end
% check if table for today's date exists
d_string = strsplit(settings.time_start,'_'); 
date_session= d_string{1};
time_session = d_string{2}; 

% session info 
columns = {'Date','Session','Offset x', 'Offset y', 'Calibration Used?', 'Stimuli', 'Trial mode',...
    'Stimuli Size','Stimuli Location', 'Bounding box','n_rsvp','Response time', 'Hold time', 'Time to initiate','total n_trials','correct n_trials','Hit rate',...
    'cumulative n_trials','reward vol','reward type','cumulative reward vol','n_trials_licked','n_imgs_seen','cumulative n_imgs_seen',...
    'experiment param file'}; 

session_file_dir = fullfile(log_dir,settings.subject,strcat(date_session,'.mat')); 
if ~exist(session_file_dir) 
    T = cell(1,length(columns));
    save(session_file_dir,'T'); 
    curridx= 1;
    cumulative_n_trials = 0; 
    cumulative_reward = 0; 
    cumulative_n_imgs_seen = 0;
else
    load(session_file_dir); 
    curridx = height(T) + 1; 
    cumulative_n_trials = sum(T.('correct n_trials')); 
    cumulative_reward = sum(T.('reward vol')); 
    cumulative_n_imgs_seen = sum(T.('n_imgs_seen')); 
    T = table2cell(T);
end

T{curridx,1} = date_session;
T{curridx,2} = time_session;
    
if isfield(calib,'gaze_offset')
    if ~isempty(calib.gaze_offset)
        T{curridx,3} = calib.gaze_offset(1);
        T{curridx,4} = calib.gaze_offset(2);
    end
end
if eyetrack.calib_applied == 1
    T{curridx,5} = strcat('X: ', strjoin(cellstr(num2str(eyetrack.cX'))),...
        ' Y: ', strjoin(cellstr(num2str(eyetrack.cY'))));
end

T{curridx,6} = calib_settings.stim_mode; 
if strcmp(T{curridx,6}, 'images')
    str_l = ''; 
    if strcmp(class(calib_settings.img_folder),'cell')
        for i = 1:length(calib_settings.img_folder)
            l = strsplit(calib_settings.img_folder{i},'\'); 
            l_bool = ~cellfun(@isempty,l);
            l = l(l_bool); 
            str_l = strcat(str_l,l{length(l)},{' '});
        end
    else
        l = strsplit(calib_settings.img_folder,'\'); 
        l_bool = ~cellfun(@isempty,l);
        l = l(l_bool); 
        str_l = l{length(l)};
    end
    
    if calib_settings.pulsed_img_size == 1
        T{curridx,6} = strcat('pulsating', {' '}, str_l);
    else
        T{curridx,6} = str_l; 
    end
end

if isfield(calib_settings,'expt_params')
    
    str_split = strsplit(calib_settings.expt_params,'\');
    T{curridx,25} = str_split{length(str_split)}; 
    
    if contains(calib_settings.expt_params,'quadrant')
        T{curridx,6} = strcat('quadrant', {' '}, T{curridx,6});
    end
    
    if contains(calib_settings.expt_params,'quadrant')
        T{curridx,10} = repmat([calib_settings.disp_rect(3)/2,calib_settings.disp_rect(4)/2],[length(calib.pts),1]); % quadrant 
    elseif strcmp(calib.reward_on,'quality')
        T{curridx,10} = NaN;
    else
        if isfield(calib_settings,'bounding_rect_size_x')
            T{curridx,10}  = [calib_settings.bounding_rect_size_x,calib_settings.bounding_rect_size_y];
        elseif isfield(calib_settings,'bounding_box_x')
            T{curridx,10} = [calib_settings.bounding_box_x,calib_settings.bounding_box_y]; 
            
        end
    end
end

if isfield(calib,'trial_mode')
    T{curridx,7} = calib.trial_mode;
end

if isfield(calib_settings,'stim_rect_size_x')
    stim_size = [calib_settings.stim_rect_size_x;calib_settings.stim_rect_size_y]; 
    T{curridx,8} = stim_size;
end

T{curridx,9} = calib.pts; 
center_screen = [calib_settings.disp_rect(3)/2, calib_settings.disp_rect(4)/2];

if isfield(calib_settings,'n_rsvp')
    T{curridx,11} = calib_settings.n_rsvp; 
end

if isfield(calib,'time_out_after')
    T{curridx,12} = calib.time_out_after;
end

if isfield(calib_settings,'n_rsvp')
    if calib_settings.n_rsvp >1
        T{curridx,13} = calib_settings.presentation_time * calib_settings.n_rsvp + calib_settings.rsvp_iti_ms * (calib_settings.n_rsvp-1); 
    else
        T{curridx,13} = calib.time_to_reward;
    end
else
    if isfield(calib,'time_to_reward')
        T{curridx,13} = calib.time_to_reward;
    end
end

if isfield(calib_settings,'require_fix_tr_init')

    if calib_settings.require_fix_tr_init == 1
        T{curridx,14} = calib_settings.fixation_to_init;
    end
end

valid_trial_ind = calib.sequence~=0; 

if datetime(date_session) < datetime('2022-10-01')
    T{curridx,15} = calib.n_completed;
else
    finished_sequence = calib.sequence(1:calib.n_completed);
    T{curridx,15} = sum(finished_sequence~=0);
end

if isfield(reward,'correct_trial')
    T{curridx,16} = sum(reward.correct_trial(calib.sequence~=0)); % correct trials
elseif isfield(reward,'n_rewards_given')  % previous Marmulator version
    T{curridx,16} = sum(reward.reward_sequence(calib.sequence~=0)) - sum(reward.reward_sequence == 1 &  calib.trial_init_timed_out ==1);  % remove trials that froze
elseif isfield(reward,'nr_rewards_given')
    T{curridx,16} = reward.nr_rewards_given; 
else
    T{curridx,16} = nan; 
end


T{curridx,17} = T{curridx,16}/T{curridx,15};
T{curridx,18} = cumulative_n_trials + T{curridx,16}; 

if isfield(reward,'nr_rewards_given')
    T{curridx,19} =reward.nr_rewards_given * reward.reward_vol;
elseif isfield(reward,'n_rewards_given')
    T{curridx,19} = reward.n_rewards_given * reward.reward_vol;
end

if isfield(reward,'man_reward_ct')
    T{curridx,19} = T{curridx,19} + reward.man_reward_ct * reward.reward_vol;
end

if isfield(reward,'bonus_reward_ct')
    T{curridx,19} = T{curridx,19} + reward.bonus_reward_ct * reward.reward_vol;
end

if isfield(reward,'reward_type')
    T{curridx,20} = reward.reward_type;
end

T{curridx,21} = cumulative_reward+ T{curridx,19}; 

if isfield(reward,'lick_trial')
    if reward.lickometer == 1
        T{curridx,22} = sum(reward.lick_trial);
    end
end

% number of successful image presentations
% eye duration 
eye_dur = calib.end_t - calib.start_t;

% 
n_imgs_seen = zeros(1,calib.n_completed); 
if isfield(calib_settings,'n_rsvp')
    if calib_settings.n_rsvp > 1
        for n = 1:calib.n_completed
            if calib.sequence(n) == 1
                for nn = 1:calib_settings.n_rsvp
                    if eye_dur(n) >= (round(calib_settings.presentation_time/1e3/(1/calib_settings.iti_frames(n)))-1)* 1/calib_settings.iti_frames(n)* nn +  (round(calib_settings.rsvp_iti_ms/1e3/(1/calib_settings.iti_frames(n)))-1)* 1/calib_settings.iti_frames(n) *(nn-1) 
                        n_imgs_seen(n) = n_imgs_seen(n) + 1; 
                    end
                end
            end
        end
    end
end

n_imgs_seen_session = sum(n_imgs_seen); 

if T{curridx,11} >1
    T{curridx,23} = n_imgs_seen_session;
    T{curridx,24} = cumulative_n_imgs_seen + n_imgs_seen_session;
else
    T{curridx,23} = 0;
    T{curridx,24} = cumulative_n_imgs_seen; 
end

T = cell2table(T); 
T.Properties.VariableNames = columns;
save(session_file_dir,'T'); 
end
