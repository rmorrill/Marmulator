function [offset_x, offset_y] = get_mean_x_y_pts(filepath)
    %eye_data_dir = 'C:\MATLAB\eyetracker_calibration_071222\calib_data';
    %eye_data_file = 'calib_TEST_2022-09-13_14-24-52.mat'; 
    %eye_data_file = 'calib_Blizzard_2022-09-13_15-22-14.mat'; 
    %eye_data_file = 'calib_Blizzard_2022-09-13_16-35-52.mat'; 
    %eye_data_file = 'calib_Blizzard_2022-09-19_14-11-06.mat'; 
    %eye_data_file = 'calib_Blizzard_2022-09-20_13-33-07.mat'; 
    %eye_data_file = 'calib_Blizzard_2022-09-21_14-13-05.mat';
    %eye_data_file = 'calib_Blizzard_2022-09-22_13-41-57.mat';
    %eye_data_file = 'calib_Blizzard_2022-10-01_15-23-39.mat';
    %eye_data_file = 'calib_Blizzard_2022-10-03_15-02-08.mat'; 
    %eye_data_file = 'calib_Blizzard_2022-10-06_14-03-59.mat'; 

    %E = load(fullfile(eye_data_dir, eye_data_file)); 
    E = load(filepath); 
    x_vals_all = E.eyetrack.eyepos_raw(:,1); 
    y_vals_all = E.eyetrack.eyepos_raw(:,2); 

    t_win = [-0.6 0]; 
    t_win_relative_to = 'offset'; % 'onset' or 'offset' 
    remove_unrewarded_trs = true; 

    n_trs_complete = E.calib.n_completed;
    stim_start_t = E.calib.start_t;
    stim_end_t = E.calib.end_t;
    eye_time = E.eyetrack.time; 

    for i = 1:n_trs_complete
        if strcmp(t_win_relative_to, 'onset')
            curr_st = stim_start_t(i) + t_win(1);
            curr_end = stim_start_t(i) + t_win(2);
        elseif strcmp(t_win_relative_to, 'offset')
            curr_st = stim_end_t(i) + t_win(1);
            curr_end = stim_end_t(i) + t_win(2);
        end
        idx_pts_tr = find(eye_time > curr_st & eye_time<=curr_end);
        these_pts_tmp{i} = [x_vals_all(idx_pts_tr) y_vals_all(idx_pts_tr)]'; 
    end

    if remove_unrewarded_trs 
        these_pts_tmp = these_pts_tmp(E.reward.reward_sequence(1:n_trs_complete) & E.calib.sequence(1:n_trs_complete)>0); 
    else % then only remove catch wake-up trials
        these_pts_tmp = these_pts_tmp(E.calib.sequence(1:n_trs_complete)>0)
    end

    use_pts = [these_pts_tmp{:}]'; 
    fprintf('%d pts prior to outlier removal\n', length(use_pts)); 

    %% manual outlier removal
    f3 = figure('Color', 'w', 'Position', [670   496   857   502]);

    clear p_mu p_sdx p_sdy
    p_mu = [];
    p_sdx = [];
    p_sdy = [];
    currcol = distinguishable_colors(1); 

    s = scatter(use_pts(:,1), use_pts(:,2), 10, currcol, 'filled');
    hold on
    s.MarkerFaceAlpha = 0.3;
    curr_sdx = std(use_pts(:,1));
    curr_mux = mean(use_pts(:,1));
    curr_sdy = std(use_pts(:,2));
    curr_muy = mean(use_pts(:,2));
    p_mu(end+1) = plot(curr_mux, curr_muy, 'Marker', '.', 'Color', currcol);
    hold on
    p_sdx(end+1) = plot([curr_mux - curr_sdx curr_mux + curr_sdx], [curr_muy curr_muy], 'Color', currcol, 'LineWidth', 2);
    p_sdy(end+1) = plot([curr_mux curr_mux], [curr_muy - curr_sdy curr_muy + curr_sdy], 'Color', currcol, 'LineWidth', 2);
    hold on

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
        x_curr = use_pts(:,1);
        y_curr = use_pts(:,2);
        in_idx = inpolygon(x_curr, y_curr, xQ, yQ);
        fprintf('Center positioin: %0.2f%% included\n', sum(in_idx)*100/numel(in_idx));

        use_pts = [x_curr(in_idx) y_curr(in_idx)];
        s = scatter(use_pts(:,1), use_pts(:,2), 10, currcol, 'filled');
        hold on
        s.MarkerFaceAlpha = 0.3;
        curr_sdx = std(use_pts(:,1));
        curr_mux = mean(use_pts(:,1));
        curr_sdy = std(use_pts(:,2));
        curr_muy = mean(use_pts(:,2));
        p_mu(end+1) = plot(curr_mux, curr_muy, 'Marker', '.', 'Color', currcol);
        hold on
        p_sdx(end+1) = plot([curr_mux - curr_sdx curr_mux + curr_sdx], [curr_muy curr_muy], 'Color', currcol, 'LineWidth', 2);
        p_sdy(end+1) = plot([curr_mux curr_mux], [curr_muy - curr_sdy curr_muy + curr_sdy], 'Color', currcol, 'LineWidth', 2);
        hold on


        axis ij

        title('Included data points with means');

    end

    %%

    figure('Color', 'w', 'Position', [438   408   492   767]);
    
    
    x_mu = mean(use_pts(:,1));
    y_mu = mean(use_pts(:,2));
    x_mu = median(use_pts(:,1));
    y_mu = median(use_pts(:,2));
    
    
    bins = linspace(0,1, 101);
    bincents = diff(bins)/2 + bins(1:end-1);
    x_n = histcounts(use_pts(:,1), bins);
    y_n = histcounts(use_pts(:,2), bins);
    
    subplot(2,1,1);
    bar(bincents, x_n);
    yL = get(gca, 'YLim');
    hold on
    ylabel('counts');
    xlabel('x location');
    plot([x_mu, x_mu], yL, 'r-');
    title(sprintf('mean x: %0.3f', x_mu));
    set(gca, 'FontSize', 14);
    box off
    
    subplot(2,1,2);
    bar(bincents, y_n);
    hold on
    yL = get(gca, 'YLim');
    plot([y_mu, y_mu], yL, 'r-');
    ylabel('counts');
    xlabel('y location');
    title(sprintf('mean y: %0.3f', y_mu));
    set(gca, 'FontSize', 14);
    box off
    
    center_pt_x_pix = (E.calib_settings.gaze_center_adj_x + E.calib_settings.disp_rect(3)/2);
    center_pt_y_pix = (E.calib_settings.gaze_center_adj_y + E.calib_settings.disp_rect(4)/2);
    center_pt_x = center_pt_x_pix/E.calib_settings.disp_rect(3);
    center_pt_y = center_pt_y_pix/E.calib_settings.disp_rect(4);
    
    fprintf('gaze center pt x: %d pix, %0.3f on screen\n', center_pt_x_pix, center_pt_x);
    fprintf('gaze center pt y: %d pix, %0.3f on screen\n\n', center_pt_y_pix, center_pt_y);
    offset_x = center_pt_x-x_mu;
    offset_y = center_pt_y-y_mu;
    fprintf('offset x: %0.3f\n', offset_x);
    fprintf('offset y: %0.3f\n', offset_y);
    

end
