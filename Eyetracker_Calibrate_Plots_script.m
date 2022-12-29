%%%% DATA LOAD 

close all; 
%session_file = 'calib_West_2022-12-18_14-18-16.mat'; 
session_file = 'calib_West_2022-12-28_13-22-26.mat'; 
session_data_dir = 'C:\Data\West\2022-12-28\calibration\';  

%%
%%%% CALIBARTION SETTINGS 
calc_calib_coeffs_on = 'all_pts'; % 'all pts' OR 'means'
%calc_calib_coeffs_on = 'means'; % 'all pts' OR 'means'
outlier_rejection = 'by_stimulus'; % 'all_pts' or 'by_stimulus'; 
remove_unrewarded_trs = true;
remove_freeze_trs = true; 
SAVE_CALIB = 1;
show_outliers_plot2 = 0; 
filter_by_eye_area = 1; 
%ea_thresh = 0.04; 
ea_thresh = 0.03; 
ar_reject_sds = 2.5; % number of aspect ratio SDs for rejection


%% THREE FIGURES: 
% F1: main 
% F2: manual outlier rejection
% F3: location-specific plotting 

%%
calibration_save_dir = fullfile(session_data_dir, 'calib_coeffs'); 
session_file_full = fullfile(session_data_dir, session_file); 

session_ts = char(regexp(session_file_full, '\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}', 'match'));  
clear sp
set(0, 'defaulttextinterpreter', 'none'); 

load(session_file_full);
if isfield(calib_settings, 'fixation_to_init')
    if ~isnan(calib_settings.fixation_to_init)
        t_win_start = (calib.time_to_reward + calib_settings.fixation_to_init)/1000;
    else
        t_win_start = calib.time_to_reward/1000;
    end
    
else
    t_win_start = calib.time_to_reward/1000;
end

%t_win = [0 1.5]; 
t_win = [-t_win_start 0]; 

t_win_relative_to = 'offset'; % 'onset' or 'offset' 


subject = settings.subject; 

x_data_raw = eyetrack.eyepos_raw(:,1)'; 
y_data_raw = eyetrack.eyepos_raw(:,2)'; 

% SET F1 n_subplots 
f1 = figure('Position', [16    62   948   1200], 'Color', 'w');
n_x = 2;
n_y = 11;

subplot_er(n_y, n_x, 1);
scatter(x_data_raw, y_data_raw, 10, eyetrack.time);
title('raw, color = time')
ylabel('y');
xlabel('x');

disp_x_len = calib_settings.disp_rect(3);
disp_y_len = calib_settings.disp_rect(4);

n_calib_pts = size(calib.pts,1);
cols = distinguishable_colors(n_calib_pts);
calib_x = calib.pts(:,1);
calib_y = calib.pts(:,2);

%% print info about session and performance
fprintf('\nSession info:\n'); 
fprintf('mode: %s, reward on: %s\n', calib.trial_mode, calib.reward_on); 
if strcmp(calib.trial_mode, 'foraging')
    fprintf('response time: %d ms\n', calib.time_out_after); 
    fprintf('hold time for reward: %d ms\n', calib.time_to_reward); 
else
    fprintf('trial duration: %d ms\n', calib_settings.presentation_time); 
end

fprintf('\nAnalysis window: %d ms to %d ms\n', round(t_win(1)*1e3), round(t_win(2)*1e3)); 
fprintf('Include only reward trs: %d\n', remove_unrewarded_trs); 
fprintf('remove trials that were frozen: %d\n', remove_freeze_trs); 
    
calib_seq = calib.sequence(1:calib.n_completed); 
%reward_seq = reward.reward_sequence(1:calib.n_completed); % previous Marmulator version
reward_seq = reward.correct_trial(1:calib.n_completed); % 2022-11-25 YJ. includes correct trials with no licks
fprintf('\nPerformance summary:\n');  
for i = 1:n_calib_pts
    n_trs_curr = sum(calib_seq == i); 
    n_trs_corr_curr = sum(calib_seq == i & reward_seq); 
    fprintf('Position %d: %0.1f%% (%d/%d)correct\n', i, n_trs_corr_curr*100/n_trs_curr, n_trs_corr_curr, n_trs_curr); 
    %fprintf('\tMean time to reward = %0.2f +/- %0.2f ms');
end

%% PLOT FOR EYE AREA PT REMOVAL
eye_area = eyetrack.pupil_size_x .* eyetrack.pupil_size_y * pi; 
ar = eyetrack.pupil_size_x./eyetrack.pupil_size_y; 
qual = eyetrack.quality; 

ar_mu = mean(ar); 
ar_std = std(ar); 
aspect_reject_thresh = ar_mu + ar_std*ar_reject_sds;  

nY3 = 3; 
nX3 = 1; 
f3 = figure('Position', [680   220   371   758], 'color', 'w', 'Name', 'Filter by eye'); 
subplot(nY3,nX3,1)
scatter(x_data_raw, y_data_raw, 10, eye_area, 'filled')
colormap(jet)
cb = colorbar; 
set(get(cb,'Title'),'String','Eye area')
hold on
xlabel('x'); 
ylabel('y'); 
title('Eye area by location'); 

subplot(nY3,nX3,2)
scatter(eye_area, ar, 10, 'k')
ylabel('aspect ratio'); 
xlabel('eye area'); 
hold on
xL = xlim; 
yL = ylim; 
plot(xL, [aspect_reject_thresh, aspect_reject_thresh], 'r--'); 
plot([ea_thresh, ea_thresh], yL, 'r--'); 

subplot(nY3, nX3, 3)
filtidx = eye_area>ea_thresh; 
s1= scatter(x_data_raw(filtidx), y_data_raw(filtidx), 10, 'filled'); 
s1.MarkerFaceAlpha = 0.2; 

hold on
title(sprintf('filter: eye_area>%0.3f\n%0.1f %% removed', ea_thresh, sum(~filtidx)*100/numel(filtidx)));  

%% outlier rejection:
figure(f1)
%outlier_thresh = 3;
outlier_thresh = 3;
pts_keep = [];

subplot_er(n_y, n_x, 2)
cla

if filter_by_eye_area
    reject_on_ar = ar> aspect_reject_thresh | eye_area <= ea_thresh; 
else
    reject_on_ar = ar> aspect_reject_thresh; 
end

if remove_unrewarded_trs
    include_trs = reward_seq; 
else
    include_trs = true(1,calib.n_completed); 
end

if remove_freeze_trs
    ind = find(reward.reward_sequence == 1 &  calib.trial_init_timed_out ==1); 
    include_trs(ind) = 0;
end

switch outlier_rejection
    case 'all_pts'
        % method 1: all data points
        these_pts = [];
        outlier_ct = 0;
        for i = 1:n_calib_pts
            currstim = find(calib_seq == i & include_trs);
            %disp(numel(currstim));
            these_pts_tmp = [];
            tr_idx_tmp = [];
            reject_on_ar_tmp = [];
            for j = 1:numel(currstim)
                if strcmp(t_win_relative_to, 'onset')
                    curr_st = calib.start_t(currstim(j)) + t_win(1);
                    curr_end = calib.start_t(currstim(j)) + t_win(2);
                elseif strcmp(t_win_relative_to, 'offset')
                    curr_st = calib.end_t(currstim(j)) + t_win(1);
                    curr_end = calib.end_t(currstim(j)) + t_win(2);
                end
                idx_pts_tr = find(eyetrack.time > curr_st & eyetrack.time<=curr_end);
                these_pts_tmp = [these_pts_tmp idx_pts_tr];
                tr_idx_tmp = [tr_idx_tmp ones(1,numel(idx_pts_tr))*currstim(j)];
                reject_on_ar_tmp = [reject_on_ar_tmp reject_on_ar(idx_pts_tr)];
                % keyboard
            end
            these_pts{i} = [these_pts_tmp; ones(1,numel(these_pts_tmp))*i; tr_idx_tmp; reject_on_ar_tmp];
        end
        
        these_pts_all = [these_pts{:}];
        reject_ar = these_pts_all(4,:); 
        curr_x = x_data_raw(these_pts_all(1,:));
        curr_y = y_data_raw(these_pts_all(1,:));
        
        outlier_x = isoutlier(curr_x, 'ThresholdFactor', outlier_thresh);
        outlier_y = isoutlier(curr_y, 'ThresholdFactor', outlier_thresh);
        outlier_both = outlier_x | outlier_y | reject_ar;
        outlier_ct = sum(outlier_both) + outlier_ct;
        
        for i = 1:n_calib_pts
            curridx = these_pts_all(2,:)  == i;
            tr_idx_all = these_pts_all(3,:); 
            %disp(numel(unique(tr_idx_all))); 
            tr_idx = tr_idx_all(curridx & ~outlier_both)'; 
            curr_x_in = curr_x(curridx & ~outlier_both)'; 
            curr_x_out = curr_x(curridx & outlier_both)'; 
            curr_y_in = curr_y(curridx & ~outlier_both)'; 
            curr_y_out = curr_y(curridx & outlier_both)'; 
            pts_keep{i} = [curr_x_in curr_y_in tr_idx];
            disp(numel(unique(tr_idx))); 
            s_in = scatter(curr_x_in, curr_y_in, 10, cols(i,:), 'filled');
            s_in.MarkerFaceAlpha = 0.3;
            
            hold on
            
            if show_outliers_plot2
                s_out = scatter(curr_x_out, curr_y_out, 10, cols(i,:), 'filled');
                s_out.MarkerFaceAlpha = 0.2;
                s_out = scatter(curr_x_out, curr_y_out, 30, [0.5 0.5 0.5], 'Marker', 'x', 'LineWidth', 0.1);
            end
            
            hold on
        end
        
    case 'by_stimulus'
        % method 2: by stimulus
        
        outlier_ct = 0;
        for i = 1:n_calib_pts
            
            currstim = find(calib_seq == i & include_trs);
            these_pts = [];
            tr_idx = [];
            for j = 1:numel(currstim)
                %disp(numel(currstim));
                if strcmp(t_win_relative_to, 'onset')
                    curr_st = calib.start_t(currstim(j)) + t_win(1);
                    curr_end = calib.start_t(currstim(j)) + t_win(2);
                elseif strcmp(t_win_relative_to, 'offset')
                    curr_st = calib.end_t(currstim(j)) + t_win(1);
                    curr_end = calib.end_t(currstim(j)) + t_win(2);
                end
                these_ts = find(eyetrack.time > curr_st & eyetrack.time<=curr_end);
                these_pts = [these_pts these_ts];
                tr_idx = [tr_idx ones(1,numel(these_ts))*currstim(j)];
                % keyboard
            end

            curr_x = x_data_raw(these_pts);
            curr_y = y_data_raw(these_pts);
            
            outlier_x = isoutlier(curr_x, 'ThresholdFactor', outlier_thresh);
            outlier_y = isoutlier(curr_y, 'ThresholdFactor', outlier_thresh);
            %outlier_x = false(1,numel(curr_x)); 
            %outlier_y = false(1,numel(curr_y)); 
            
            outlier_both = outlier_x | outlier_y;
            outlier_ct = sum(outlier_both) + outlier_ct;
            scatter(curr_x(~outlier_both), curr_y(~outlier_both), 10, cols(i,:));
            hold on
            if show_outliers_plot2
                s_out = scatter(curr_x(outlier_both), curr_y(outlier_both), 10, cols(i,:));
                %s_out.MarkerEdgeAlpha = 0.5;
                s_out = scatter(curr_x(outlier_both), curr_y(outlier_both), 30, [0.5 0.5 0.5], 'Marker', 'x', 'LineWidth', 0.1);
                hold on
            end
            pts_keep{i} = [curr_x(~outlier_both)' curr_y(~outlier_both)' tr_idx(~outlier_both)'];
            disp(numel(unique(tr_idx(~outlier_both))));
            %disp(unique(tr_idx));
            
        end
end
axis ij
ylabel('y');
xlabel('x');
box off
title(sprintf('outliers: %s, thresh = %d MAD\n%0.2f %% (%d) outliers AUTO removed', outlier_rejection, outlier_thresh, outlier_ct*100/numel(x_data_raw), outlier_ct));

%% manual outlier removal
f3 = figure('Color', 'w', 'Position', [670   496   857   502]);

clear p_mu p_sdx p_sdy
p_mu = []; 
p_sdx = []; 
p_sdy = []; 
for i = 1:n_calib_pts
    if ~isempty(pts_keep{i})
        s = scatter(pts_keep{i}(:,1), pts_keep{i}(:,2), 10, cols(i,:), 'filled');
        hold on
        s.MarkerFaceAlpha = 0.3;
        curr_sdx = std(pts_keep{i}(:,1));
        curr_mux = median(pts_keep{i}(:,1));
        curr_sdy = std(pts_keep{i}(:,2));
        curr_muy = median(pts_keep{i}(:,2));
        p_mu(end+1) = plot(curr_mux, curr_muy, 'Marker', '.', 'Color', cols(i,:));
        hold on
        p_sdx(end+1) = plot([curr_mux - curr_sdx curr_mux + curr_sdx], [curr_muy curr_muy], 'Color', cols(i,:), 'LineWidth', 2);
        p_sdy(end+1) = plot([curr_mux curr_mux], [curr_muy - curr_sdy curr_muy + curr_sdy], 'Color', cols(i,:), 'LineWidth', 2);
        hold on
    end
end
uistack([p_mu p_sdx p_sdy], 'top'); 
axis ij
xlabel('x');
ylabel('y'); 
title('Draw polygon around points for inclusion, double click to bypass');
roi = drawpolygon(gca);

if ~isempty(roi.Position)
    
    xQ = roi.Position(:,1);
    yQ = roi.Position(:,2);
    
    f4 = figure('Color', 'w', 'Position', [670   496   857   502]);
    
    p_mu = [];
    p_sdx = [];
    p_sdy = [];
    for i = 1:n_calib_pts
        if ~isempty(pts_keep{i})
            x_curr = pts_keep{i}(:,1);
            y_curr = pts_keep{i}(:,2);
            in_idx = inpolygon(x_curr, y_curr, xQ, yQ);
            fprintf('Position: %d: %0.2f %% included\n',i,  sum(in_idx)*100/numel(in_idx));
            
            pts_keep{i} = [x_curr(in_idx) y_curr(in_idx) pts_keep{i}(in_idx,3)];
            s = scatter(pts_keep{i}(:,1), pts_keep{i}(:,2), 10, cols(i,:), 'filled');
            hold on
            s.MarkerFaceAlpha = 0.3;
            curr_sdx = std(pts_keep{i}(:,1));
            curr_mux = mean(pts_keep{i}(:,1));
            curr_sdy = std(pts_keep{i}(:,2));
            curr_muy = mean(pts_keep{i}(:,2));
            p_mu(end+1) = plot(curr_mux, curr_muy, 'Marker', '.', 'Color', cols(i,:));
            hold on
            p_sdx(end+1) = plot([curr_mux - curr_sdx curr_mux + curr_sdx], [curr_muy curr_muy], 'Color', cols(i,:), 'LineWidth', 2);
            p_sdy(end+1) = plot([curr_mux curr_mux], [curr_muy - curr_sdy curr_muy + curr_sdy], 'Color', cols(i,:), 'LineWidth', 2);
            hold on
        end
        
    end
    axis ij
    
    title('Included data points with means');
    
end

%% add a second figure for location-specific outlier removal
f2 = figure('Position', [1490  397 1531 889], 'Color', 'w');
%keyboard

n_y2 = calib.n_pts_y;
n_x2 = calib.n_pts_x;
clear sp 
for i = 1:n_calib_pts
    pts_curr = pts_keep{i}; 
    if ~isempty(pts_curr)  
        %sp(i) = subplot_er(n_y2, n_x2, i);
        sp(i) = subplot(n_y2, n_x2, i);
        trs_un = unique(pts_curr(:,3));
        for j = 1:numel(trs_un)
            idx_curr_tr = find(pts_curr(:,3) == trs_un(j));
            pts_curr_x = pts_curr(idx_curr_tr,1);
            pts_curr_y = pts_curr(idx_curr_tr,2);
            
            s1 = scatter(pts_curr_x, pts_curr_y, 10, cols(i,:), 'filled');
            hold on
            s1.MarkerFaceAlpha = 0.2;
            
            mu_x = median(pts_curr_x);
            mu_y = median(pts_curr_y);
            std_x = std(pts_curr_x);
            std_y = std(pts_curr_y);
            plot(mu_x, mu_y, 'Marker', '.', 'Color', cols(i,:));
            hold on
            plot([mu_x - std_x mu_x + std_x], [mu_y mu_y], 'Color', cols(i,:));
            plot([mu_x mu_x], [mu_y - std_y mu_y + std_y], 'Color', cols(i,:));
            text(mu_x, mu_y, sprintf('%d', trs_un(j)), 'Color', cols(i,:), ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
            
            set(gca, 'XGrid', 'on', 'YGrid', 'on')
        end
    end
    title(sprintf('%d', i));
    axis tight
    axis square
    axis ij
    box off
end
if ~isnan(calib_settings.fixation_to_init)
    t_win_start = (calib.time_to_reward + calib_settings.fixation_to_init)/1000;
else
    t_win_strat = calib.time_to_reward/1000;
end
sp = sp(isgraphics(sp)); 
xL_all = get(sp, 'XLim');
xL_set = [min([xL_all{:}]) max([xL_all{:}])];
yL_all = get(sp, 'YLim');
yL_set = [min([yL_all{:}]) max([yL_all{:}])];


set(sp, 'YLim', yL_set);
set(sp, 'XLim', xL_set);

sgtitle('Points included by location')
%
xL_all = get(sp, 'XLim');
xL_set = [min([xL_all{:}]) max([xL_all{:}])];
yL_all = get(sp, 'YLim');
yL_set = [min([yL_all{:}]) max([yL_all{:}])];


set(sp, 'YLim', yL_set);
set(sp, 'XLim', xL_set);

%% OR MANUALLY SET THE X AND Y LIMS
%xL_set = [0.5 0.6];
%yL_set = [0.35  0.5];

%set(sp, 'YLim', yL_set);
%set(sp, 'XLim', xL_set);

%%
% remove outliers

figure(f1);

subplot_er(n_y, n_x, 3);
mu_x_all = zeros(1,n_calib_pts); 
mu_y_all = zeros(1,n_calib_pts);
mu_xy_all = zeros(1,n_calib_pts); 
for i = 1:n_calib_pts
    if ~isempty(pts_keep{i})
    pts_curr_x = pts_keep{i}(:,1);
    pts_curr_y = pts_keep{i}(:,2);
    %mu_x = median(pts_curr_x);
    %mu_y = median(pts_curr_y);
    mu_x = mean(pts_curr_x);
    mu_y = mean(pts_curr_y);
    std_x = std(pts_curr_x);
    std_y = std(pts_curr_y);
    plot(mu_x, mu_y, 'Marker', '.', 'Color', cols(i,:));
    hold on
    plot([mu_x - std_x mu_x + std_x], [mu_y mu_y], 'Color', cols(i,:));
    plot([mu_x mu_x], [mu_y - std_y mu_y + std_y], 'Color', cols(i,:));
    text(mu_x, mu_y, sprintf('%d', i), 'Color', cols(i,:), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    mu_x_all(i) = mu_x;
    mu_y_all(i) = mu_y;
    mu_xy_all(i) = mu_x * mu_y; 
    end
end
box off
axis ij

title('Means +/- SDs'); 
%%
subplot_er(n_y, n_x, 4)

for i = 1:n_calib_pts
    plot(calib.pts(i, 1), calib.pts(i,2), 'Color', cols(i,:), 'Marker', 's', 'MarkerSize', 10, 'LineWidth', 1.5)
    hold on
    text(calib.pts(i, 1), calib.pts(i,2), sprintf('%d', i), 'Color', cols(i,:), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
axis ij
xlim([0 disp_x_len])
ylim([0 disp_y_len])
box off
ylabel('pixel y');
xlabel('pixel x');

%% X regress plot

subplot_er(n_y, n_x, 5);
scatter(mu_x_all(mu_x_all~=0), calib_x(mu_x_all~=0), 20, cols(mu_x_all~=0, :), 'filled')
title('x regress');
xvals_eval = linspace(min(xlim), max(xlim), 10);
for i = 1:n_calib_pts
    if mu_x_all(i) ~= 0
        text(mu_x_all(i), calib.pts(i,1), sprintf('%d', i), 'Color', cols(i,:), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end
ylabel('stim location on x');
xlabel('eyetracker measured x');

hold on

%% Y regress plot
subplot_er(n_y, n_x, 6);
scatter(mu_y_all(mu_y_all~=0), calib_y(mu_y_all~=0), 20, cols(mu_y_all~=0,:), 'filled')
title('y regress');
yvals_eval = linspace(min(xlim), max(xlim), 10);
for i = 1:n_calib_pts
    if mu_y_all(i) ~= 0 
    text(mu_y_all(i), calib.pts(i,2), sprintf('%d', i), 'Color', cols(i,:), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end
hold on
ylabel('stim location on y');
xlabel('eyetracker measured y'); 
%axis square

%% PLOT DIFFERENT TRANSFORMS 
tts = 0.3; % transform text spacing

 %% PROJECTIVE TRANSFORM
 lin_reg_mode = 'projective';
 exclude_vec = [];

 p = [];
 q = [];
 for i = 1:n_calib_pts
     if ~isempty(pts_keep{i})
         p = [p;mu_x_all(i),mu_y_all(i)];
         q = [q;calib.pts(i,:)];
     end
 end

if size(q,1)>=4
    v = homography_solve(p,q); 
    % all_pts 
    % p = [all_x_pts, all_y_pts]; 
    % q = [calib_x_all, calib_y_all]; 
    coeff_X = v(:,1); 
    coeff_Y = v(:,2);
    coeff_Proj = v(:,3);

    subplot_er(n_y, n_x, 15)
    plotMuEvalPts(lin_reg_mode, calib, calib_settings, coeff_X, coeff_Y, mu_x_all, mu_y_all, mu_xy_all, cols, v)

    %%
    subplot_er(n_y, n_x, 16)
    writeCoeffsTextIntoBox(coeff_X, coeff_Y, coeff_Proj)

    subplot_er(n_y, n_x, 17)
    plotAllEvalPts(calib_settings, pts_keep, lin_reg_mode, coeff_X, coeff_Y, cols, mu_x_all, mu_y_all, mu_xy_all, v)
    
%     for i = 1:n_calib_pts
%          if ~isempty(pts_keep{i})
% 
%             %text(calib.pts(i, 1), calib.pts(i,2), sprintf('%d', i), 'Color', cols(i,:), ...
%             %    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
% 
%             eval_curr = homography_transform(pts_keep{i}(:,1:2),v);
%             xeval_curr = eval_curr(:,1);
%             yeval_curr = eval_curr(:,2); 
%                 s2 = scatter(xeval_curr, yeval_curr, 10, cols(i,:), 'filled'); 
%             s2.MarkerFaceAlpha = 0.2;
%     %         if any(ismember(exclude_vec, i))
%     %             plot(xeval_curr, yeval_curr, 'x', 'Color', [1 0 0], 'MarkerSize', 10);
%     %         end
%          end
%         hold on
% 
%     end
% 
% 
%     axis ij
%     plot([0, disp_x_len], [0, 0], 'k-'); 
%     plot([0, 0], [0, disp_y_len], 'k-'); 
%     plot([0, disp_x_len], [disp_y_len, disp_y_len], 'k-'); 
%     plot([disp_x_len, disp_x_len], [0, disp_y_len], 'k-'); 
% 
%     ylabel('screen y'); 
%     xlabel('screen x');
%     title('all (included) points: calibration applied');
%     
    %%
    if SAVE_CALIB
        saveCalib(coeff_X, coeff_Y, coeff_Proj, pts_keep, calib, calib_settings, settings,...
            lin_reg_mode, calibration_save_dir, session_ts, calc_calib_coeffs_on,...
            outlier_rejection, outlier_thresh, filter_by_eye_area, ea_thresh, aspect_reject_thresh)
    end
end

%% Gather the points for regression 
% this uses all the points
switch calc_calib_coeffs_on
    case 'means'
        inc_vec = true(1,size(calib.pts,1));
        inc_vec(exclude_vec) = false;
        
        calib_x_locs = calib.pts(inc_vec,1);
        calib_y_locs = calib.pts(inc_vec,2);
        
        x_in = mu_x_all(inc_vec);
        y_in = mu_y_all(inc_vec);
        xy_in = x_in.*y_in; 
    case 'all_pts'
        calib_x_pts = calib.pts(:,1);
        calib_y_pts = calib.pts(:,2);
        
        calib_x_locs = [];
        calib_y_locs = [];
        x_in = [];
        y_in = [];
        xy_in = [];
        for i = 1:n_calib_pts
            if ~isempty(pts_keep{i})
                n_pts_curr = length(pts_keep{i}(:,1));
                x_in = [x_in; pts_keep{i}(:,1)];
                y_in = [y_in; pts_keep{i}(:,2)];
                xy_in = [xy_in; pts_keep{i}(:,1).*pts_keep{i}(:,2)];
                calib_x_locs = [calib_x_locs; ones(n_pts_curr, 1)*calib_x_pts(i)];
                calib_y_locs = [calib_y_locs; ones(n_pts_curr, 1)*calib_y_pts(i)];
            end
        end
end


%% CROSS TERM
lin_reg_mode = 'cross term';
coeff_Proj = []; 
[coeff_X, coeff_Y] = calcTransform(lin_reg_mode, calib_x_locs, calib_y_locs, x_in, y_in, xy_in); 

figure(f1)
subplot_er(n_y, n_x, 7); 
plotMuEvalPts(lin_reg_mode, calib, calib_settings, coeff_X, coeff_Y, mu_x_all, mu_y_all, mu_xy_all, cols)

subplot_er(n_y, n_x, 8)
writeCoeffsTextIntoBox(coeff_X, coeff_Y); 

subplot_er(n_y, n_x, 9)
plotAllEvalPts(calib_settings, pts_keep, lin_reg_mode, coeff_X, coeff_Y, cols, mu_x_all, mu_y_all, mu_xy_all)

if SAVE_CALIB
    saveCalib(coeff_X, coeff_Y, coeff_Proj, pts_keep, calib, calib_settings, settings,...
        lin_reg_mode, calibration_save_dir, session_ts, calc_calib_coeffs_on,...
        outlier_rejection, outlier_thresh, filter_by_eye_area, ea_thresh, aspect_reject_thresh)
end


%% NO CROSS TERM
lin_reg_mode = 'no cross term';
coeff_Proj = []; 
[coeff_X, coeff_Y] = calcTransform(lin_reg_mode, calib_x_locs, calib_y_locs, x_in, y_in, xy_in);  

subplot_er(n_y, n_x, 11)
plotMuEvalPts(lin_reg_mode, calib, calib_settings, coeff_X, coeff_Y, mu_x_all, mu_y_all, mu_xy_all, cols);

subplot_er(n_y, n_x, 12); 
writeCoeffsTextIntoBox(coeff_X, coeff_Y); 

subplot_er(n_y, n_x, 13)
plotAllEvalPts(calib_settings, pts_keep, lin_reg_mode, coeff_X, coeff_Y, cols, mu_x_all, mu_y_all, mu_xy_all)


if SAVE_CALIB
    saveCalib(coeff_X, coeff_Y, coeff_Proj, pts_keep, calib, calib_settings, settings,...
        lin_reg_mode, calibration_save_dir, session_ts, calc_calib_coeffs_on,...
        outlier_rejection, outlier_thresh, filter_by_eye_area, ea_thresh, aspect_reject_thresh)
end


%% CROSS TERM AND XY TERM
lin_reg_mode = 'cross term and xy term';
coeff_Proj = []; 

[coeff_X, coeff_Y] = calcTransform(lin_reg_mode, calib_x_locs, calib_y_locs, x_in, y_in, xy_in); 

subplot_er(n_y, n_x, 19)
plotMuEvalPts(lin_reg_mode, calib, calib_settings, coeff_X, coeff_Y, mu_x_all, mu_y_all, mu_xy_all, cols);

subplot_er(n_y, n_x, 20)
writeCoeffsTextIntoBox(coeff_X, coeff_Y)

subplot_er(n_y, n_x, 21)
plotAllEvalPts(calib_settings, pts_keep, lin_reg_mode, coeff_X, coeff_Y, cols, mu_x_all, mu_y_all, mu_xy_all)


if SAVE_CALIB
    saveCalib(coeff_X, coeff_Y, coeff_Proj, pts_keep, calib, calib_settings, settings,...
        lin_reg_mode, calibration_save_dir, session_ts, calc_calib_coeffs_on,...
        outlier_rejection, outlier_thresh, filter_by_eye_area, ea_thresh, aspect_reject_thresh)
end

%%

sgtitle(session_file)


%% plot means
function plotMuEvalPts(lin_reg_mode, calib, calib_settings, coeff_X, coeff_Y, mu_x_all, mu_y_all, mu_xy_all, cols, v)

if nargin == 9 
    v = []; 
end

n_calib_pts = size(calib.pts, 1);
disp_x_len = calib_settings.disp_rect(3);
disp_y_len = calib_settings.disp_rect(4); 
for i = 1:n_calib_pts
    
    plot(calib.pts(i, 1), calib.pts(i,2), 'Color', cols(i,:), 'Marker', 's', 'MarkerSize', 20, 'LineWidth', 1.5)
    hold on
    if strcmp(lin_reg_mode, 'cross term')
        xeval_curr = mu_x_all(i)*coeff_X(2) + mu_y_all(i)*coeff_X(3) + coeff_X(1);
        yeval_curr = mu_y_all(i)*coeff_Y(2) + mu_x_all(i)*coeff_Y(3) + coeff_Y(1);
    elseif strcmp(lin_reg_mode, 'no cross term')
        xeval_curr = mu_x_all(i)*coeff_X(2) + coeff_X(1);
        yeval_curr = mu_y_all(i)*coeff_Y(2)  + coeff_Y(1);
    elseif strcmp(lin_reg_mode,'cross term and xy term')
        xeval_curr = mu_x_all(i)*coeff_X(2) + mu_y_all(i)*coeff_X(3) + mu_xy_all(i) * coeff_X(4) + coeff_X(1);
        yeval_curr = mu_y_all(i)*coeff_Y(2) + mu_x_all(i)*coeff_Y(3) + mu_xy_all(i) * coeff_Y(4) + coeff_Y(1);
    elseif strcmp(lin_reg_mode, 'projective')
        eval_curr = homography_transform([mu_x_all(i),mu_y_all(i)],v);
        xeval_curr = eval_curr(:,1);
        yeval_curr = eval_curr(:,2);
        
    end
    
    plot(xeval_curr, yeval_curr, '.', 'Color', cols(i,:), 'MarkerSize', 20)
    
    %     if any(ismember(exclude_vec, i))
    %         plot(xeval_curr, yeval_curr, 'x', 'Color', [1 0 0], 'MarkerSize', 10);
    %     end
end

axis ij
xlim([0 disp_x_len])
ylim([0 disp_y_len])
box off
ylabel('pixel y');
%xlabel('pixel x');

%title(sprintf('calibration applied. %s . exlude points: %s', lin_reg_mode,sprintf('%d, ', exclude_vec)));
title(sprintf('calibration applied. %s ', lin_reg_mode));

end

%% plot all eval points
function plotAllEvalPts(calib_settings, pts_keep, lin_reg_mode, coeff_X, coeff_Y, cols, mu_x_all, mu_y_all, mu_xy_all, v)

if nargin == 9 
    v = []; 
end

disp_x_len = calib_settings.disp_rect(3);
disp_y_len = calib_settings.disp_rect(4);

for i = 1:numel(pts_keep)
    if ~isempty(pts_keep{i})
        x_raw = pts_keep{i}(:,1);
        y_raw = pts_keep{i}(:,2);
        xy_raw = pts_keep{i}(:,1) .* pts_keep{i}(:,2);
        if strcmp(lin_reg_mode,'cross term')
            xeval = x_raw*coeff_X(2) + y_raw*coeff_X(3) + coeff_X(1);
            yeval = y_raw*coeff_Y(2) + x_raw*coeff_Y(3) + coeff_Y(1);
        elseif strcmp(lin_reg_mode,'no cross term')
            xeval = x_raw*coeff_X(2)  + coeff_X(1);
            yeval = y_raw*coeff_Y(2)  + coeff_Y(1);
        elseif strcmp(lin_reg_mode,'cross term and xy term')
            xeval = x_raw*coeff_X(2) + y_raw*coeff_X(3) + xy_raw * coeff_X(4) + coeff_X(1);
            yeval = y_raw*coeff_Y(2) + x_raw*coeff_Y(3) + xy_raw * coeff_Y(4) + coeff_Y(1);
        elseif strcmp(lin_reg_mode, 'projective')
            eval_curr = homography_transform(pts_keep{i}(:,1:2),v);
            xeval = eval_curr(:,1);
            yeval = eval_curr(:,2);
        end
        
        s2 = scatter(xeval, yeval, 10, cols(i,:), 'filled');
        s2.MarkerFaceAlpha = 0.2;
        hold on
    end
end

if strcmp(lin_reg_mode,'cross term')
    xeval_means = mu_x_all*coeff_X(2) + mu_y_all*coeff_X(3) + coeff_X(1);
    yeval_means = mu_y_all*coeff_Y(2) + mu_x_all*coeff_Y(3) + coeff_Y(1);
elseif strcmp(lin_reg_mode,'no cross term')
    xeval_means = mu_x_all*coeff_X(2) + coeff_X(1);
    yeval_means = mu_y_all*coeff_Y(2) +  coeff_Y(1);
elseif strcmp(lin_reg_mode, 'cross term and xy term')
    xeval_means = mu_x_all*coeff_X(2) + mu_y_all*coeff_X(3) + mu_xy_all * coeff_X(4) + coeff_X(1);
    yeval_means = mu_y_all*coeff_Y(2) + mu_x_all*coeff_Y(3) + mu_xy_all * coeff_Y(4) + coeff_Y(1);
elseif strcmp(lin_reg_mode, 'projective')
    eval_curr = homography_transform([mu_x_all' mu_y_all'],v);
    xeval_means = eval_curr(:,1);
    yeval_means = eval_curr(:,2);
end

hold on
scatter(xeval_means(mu_x_all~=0), yeval_means(mu_x_all~=0), 20, cols((mu_x_all~=0),:))
axis ij

% plot screen coords

plot([0, disp_x_len], [0, 0], 'k-');
plot([0, 0], [0, disp_y_len], 'k-'); 
plot([0, disp_x_len], [disp_y_len, disp_y_len], 'k-'); 
plot([disp_x_len, disp_x_len], [0, disp_y_len], 'k-'); 

ylabel('screen y'); 
%xlabel('screen x'); 
title('all (included) points: calibration applied'); 
end

%% write coeffs into box
function writeCoeffsTextIntoBox(coeff_X, coeff_Y, coeff_Proj)
tts = 0.3; 
fs = 13; 
for i = 1:length(coeff_X)
    text(0.1, tts*(4-i), sprintf('cX(%d) = %0.3f', i, coeff_X(i)));
    text(0.4, tts*(4-i), sprintf('cY(%d) = %0.3f', i, coeff_Y(i)));
    if nargin == 3
        text(0.8, tts*(4-i), sprintf('cProj(%d) = %0.3f', i, coeff_Proj(i)));
    end
    axis tight
end
axis tight
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
box off;
title('coeffs', 'FontSize', 14)
end

%% calculate transform
function [coeff_X, coeff_Y] = calcTransform(lin_reg_mode, calib_x_locs, calib_y_locs, x_in, y_in, xy_in)

if strcmp(lin_reg_mode,'cross term')
    X = [ones(length(x_in),1) x_in y_in];
    Y = [ones(length(x_in),1) y_in x_in];
elseif strcmp(lin_reg_mode, 'no cross term')
    X = [ones(length(x_in),1) x_in ];
    Y = [ones(length(x_in),1) y_in ];
elseif strcmp(lin_reg_mode,'cross term and xy term')
    X = [ones(length(x_in),1), x_in, y_in, xy_in];
    Y = [ones(length(y_in),1), y_in, x_in, xy_in];
end

coeff_X = X\calib_x_locs;
coeff_Y = Y\calib_y_locs;

end

%% save calibration
function saveCalib(coeff_X, coeff_Y, coeff_Proj, pts_keep, calib, calib_settings, settings,...
    lin_reg_mode, calibration_save_dir, session_ts, calc_calib_coeffs_on,...
    outlier_rejection, outlier_thresh, filter_by_eye_area, ea_thresh, aspect_reject_thresh)

subject = settings.subject;

if ~isempty(coeff_Proj)
    coeffs.coeff_Proj = coeff_Proj;
end
coeffs.x = coeff_X;
coeffs.y = coeff_Y;

useidx = ~cellfun(@isempty, pts_keep);
calib_pts = calib.pts(useidx,:);
disp_rect = calib_settings.disp_rect;
calib_savename = sprintf('calibCoeffs_%s_%s_%s.mat', subject, session_ts,lin_reg_mode);
calib_save_full = fullfile(calibration_save_dir, calib_savename);
save_timestamp = datestr(now);

if ~isdir(calibration_save_dir); mkdir(calibration_save_dir); end

save(calib_save_full, 'coeff_X', 'coeff_Y', 'calc_calib_coeffs_on',...
    'lin_reg_mode', 'outlier_rejection', 'outlier_thresh', 'filter_by_eye_area', ...
    'ea_thresh', 'aspect_reject_thresh', 'save_timestamp');

fprintf('calibration (%s) saved to %s\n', lin_reg_mode, calib_save_full);

end

