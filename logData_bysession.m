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
    'Stimuli Size','Stimuli Location', 'Bounding box','Response time', 'Hold time', 'Time to initiate','total n_trials','correct n_trials','Hit rate',...
    'cumulative n_trials'}; 

session_file_dir = fullfile(log_dir,settings.subject,strcat(date_session,'.mat')); 
if ~exist(session_file_dir) 
    T = cell(1,length(columns));
    save(session_file_dir,'T'); 
    curridx= 1;
    cumulative_n_trials = 0; 
else
    load(session_file_dir); 
    curridx = height(T) + 1; 
    cumulative_n_trials = sum(T.('correct n_trials')); 
    T = table2cell(T);
end

T{curridx,1} = date_session;
T{curridx,2} = time_session;
    
if ~isempty(calib.gaze_offset)
    T{curridx,3} = calib.gaze_offset(1);
    T{curridx,4} = calib.gaze_offset(2);
end
if eyetrack.calib_applied == 1
    T{curridx,5} = strcat('X: ', strjoin(cellstr(num2str(eyetrack.cX'))),...
        ' Y: ', strjoin(cellstr(num2str(eyetrack.cY'))));
end

T{curridx,6} = calib_settings.stim_mode; 
if strcmp(T{curridx,6}, 'images')
    l = strsplit(calib_settings.img_folder,'\'); 
    l_bool = ~cellfun(@isempty,l);
    l = l(l_bool); 
    if calib_settings.pulsed_img_size == 1
        T{curridx,6} = strcat('pulsating', {' '}, l{length(l)});
    else
        T{curridx,6} = l{length(l)}; 
    end
end

if contains(calib_settings.expt_params,'quadrant')
    T{curridx,6} = strcat('quadrant', {' '}, T{curridx,6});
end

T{curridx,7} = calib.trial_mode;
stim_size = [calib_settings.stim_rect_size_x;calib_settings.stim_rect_size_y]; 
T{curridx,8} = stim_size;
T{curridx,9} = calib.pts; 
center_screen = [calib_settings.disp_rect(3)/2, calib_settings.disp_rect(4)/2];

if contains(calib_settings.expt_params,'quadrant')
    T{curridx,10} = repmat([calib_settings.disp_rect(3)/2,calib_settings.disp_rect(4)/2],[length(calib.pts),1]); % quadrant 
elseif strcmp(calib.reward_on,'quality')
    T{curridx,10} = NaN;
else
    if isfield(calib_settings,'bounding_rect_size_x')
        T{curridx,10}  = [calib_settings.bounding_rect_size_x,calib_settings.bounding_rect_size_y];
    else
        % previous data didn't store bounding box 
        load(calib_settings.expt_params);
        T{curridx,10} = [bounding_rect_size_x,bounding_rect_size_y];
    end
end
T{curridx,11} = calib.time_out_after;
T{curridx,12} = calib.time_to_reward; 
if isfield(calib_settings,'require_fix_tr_init')

    if calib_settings.require_fix_tr_init == 1
        T{curridx,13} = calib_settings.fixation_to_init;
    end
end

T{curridx,14} = calib.n_completed;
T{curridx,15} = reward.nr_rewards_given - sum(reward.reward_sequence == 1 &  calib.trial_init_timed_out ==1);  % remove trials that froze 

T{curridx,16} = T{curridx,15}/calib.n_completed;
T{curridx,17} = cumulative_n_trials + T{curridx,15}; 
T = cell2table(T); 
T.Properties.VariableNames = columns;
save(session_file_dir,'T'); 

end