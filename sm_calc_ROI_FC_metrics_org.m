function sm_calc_ROI_FC_metrics(varargin)
global h
% seed_idx = h.inv_soln(h.current_inv_soln).plv_seed_idx;     % seed dipole indices of current lead fields
% seed_idx1 = h.inv_soln(h.current_inv_soln).classifier_metrics.Hits; % hit seed dipole indices in lead field
% seed_idx2 = h.inv_soln(h.current_inv_soln).plv_seed_idx;
% seed_idx = nan(size(seed_idx1)); % putting plv_seed_idx in order of hits; needed when misses occur
% for v=1:length(seed_idx1) 
%     if ~isnan(seed_idx1(v))
%         v_idx = find_nearest_voxel(h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(seed_idx1(v),:),h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(seed_idx2,:));
%         seed_idx(v) = seed_idx2(v_idx); 
%     end
% end
seed_idx = h.inv_soln(h.current_inv_soln).plv_seed_idx;

comp_idx  = h.inv_soln(h.current_inv_soln).plv_comp_idx;    % comparison dipole indices of current lead fields
plv_idx = h.inv_soln(h.current_inv_soln).plv_contrast_idx;  % comparison indices for contrasts of "seed_swf"  and "comp_swf" below; clmn 1 = seed_idx;  clmn 2 = comp_idx
plv_contrasts = h.inv_soln(h.current_inv_soln).plv_contrasts;

%% calculating ROI metrics
if isfield(h.cfg.source,'TFR_results')
    if ~isempty(h.cfg.source.TFR_results)        
        %% %%%%% TFR Total Power
        if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'avg_seed_wt')
            if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt)
                %% prepare TFR data for ROI metrics
                if h.radio_3D_PLV_PLI_ROI_stepwise_stats.Value == 1 || h.radio_3D_PLV_PLI_ROI_surg_std_stats.Value == 1 % calculate stepwise ROC stats for ROI
                    %% getting seedPLV idx only = seed_plv_idx
                    true_idx = find(~isnan(h.inv_soln(h.current_inv_soln).plv_seed_idx));     % matching up with FC_true indices so the true PLV source1v2 matches with hit seed source1v2
                    
                    FC_true = permute(h.cfg.source.TFR_results.avg_true_wt,[1 3 2]);
                    FC_seed = zeros(size(FC_true));  % setting all missed source PLV to zeros
%                     FC_seed(:,true_idx,:) = permute(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt(:,:,true_idx),[1 3 2]);
                    FC_seed(:,true_idx,:) = permute(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt(:,:,true_idx),[1 3 2]);
                    FC_comp = permute(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt,[1 3 2]);
                    tfr_comp_indices = comp_idx; % leadfield indices ordered by tfr_comp_idx;
                    
                    tfr_true_locs = h.cfg.source.vx_locs;
                    tfr_seed_locs = nan(size(FC_seed,2),3);
                    tfr_seed_locs(~isnan(seed_idx),:) = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(seed_idx(~isnan(seed_idx)),:);
                    tfr_comp_locs = nan(size(FC_comp,2),3);
                    tfr_comp_locs(~isnan(comp_idx),:) = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(comp_idx(~isnan(comp_idx)),:);
                    
                    vol =[];
                    ROI_freq_int = str2num(h.edit_3D_PLV_PLI_ROI_freq_int.String);
                    ROI_time_int = str2num(h.edit_3D_PLV_PLI_ROI_time_int.String);
                    freq_int = ROI_freq_int;
                    time_int = ROI_time_int;
                    freqs = h.inv_soln(h.current_inv_soln).TFR_results.TFR_freqs;
                    lat = h.cfg.study.lat_sim;
                    plot_flag = 0;
                    fc_scale = [];
                    fc_map_scale = [];
                    ROI_only_flag=1; % only calculating for ROI
                end
                %% Step-wise thresholding for TFR ROI metrics
                if h.radio_3D_PLV_PLI_ROI_stepwise_stats.Value == 1 % calculate stepwise ROC stats for ROI
                    thresh_vals = [1e-10:str2num(h.edit_3D_PLV_PLI_ROI_step_size.String):1]*10; %10 dB
                    
                    [h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics] = calc_MCC_TFR_thresh(FC_true, FC_seed, FC_comp,...
                        freq_int, time_int, ROI_freq_int, ROI_time_int, thresh_vals, freqs, lat, plot_flag, ROI_only_flag);
                    % Fitting FPR & TPR data to create fitted ROC curve
                    [xData, yData] = prepareCurveData( squeeze([h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.ROI_perf.FPR]),...
                        squeeze([h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.ROI_perf.TPR]) );
                    if length(xData)>1 % data exists
                        ft = fittype( 'smoothingspline' ); opts = fitoptions( 'Method', 'SmoothingSpline' ); opts.Normalize = 'on'; opts.SmoothingParam = 0.9999;
                        % Fit model to data.
                        [fitresult, gof] = fit( xData, yData, ft, opts );
                        x_plv= 0:.01:1; y_plv = feval(fitresult,x_plv); y_plv(y_plv>1)=1;
                        h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR = x_plv;
                        h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.fitted_ROI_TPR = y_plv;
                        h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.fitted_ROI_AUC = trapz(x_plv,y_plv);
                    else
                        h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR = [];
                        h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.fitted_ROI_TPR = [];
                        h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.fitted_ROI_AUC = [];
                    end
                end
            end
        end
        %% %%%%% TFR Evoked Power
        if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'avg_seed_wt_evk')
            if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_evk)
                %% prepare TFR data for ROI metrics
                if h.radio_3D_PLV_PLI_ROI_stepwise_stats.Value == 1 || h.radio_3D_PLV_PLI_ROI_surg_std_stats.Value == 1 % calculate stepwise ROC stats for ROI
                    %% getting seedPLV idx only = seed_plv_idx
                     true_idx = find(~isnan(h.inv_soln(h.current_inv_soln).plv_seed_idx));     % matching up with FC_true indices so the true PLV source1v2 matches with hit seed source1v2
                    
                    FC_true = permute(h.cfg.source.TFR_results.avg_true_wt_evk,[1 3 2]);
                    FC_seed = zeros(size(FC_true));  % setting all missed source PLV to zeros
%                     FC_seed(:,true_idx,:) = permute(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt(:,:,true_idx),[1 3 2]);
                    FC_seed(:,true_idx,:) = permute(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_evk(:,:,true_idx),[1 3 2]);
                    FC_comp = permute(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_evk,[1 3 2]);
                    tfr_comp_indices = comp_idx; % leadfield indices ordered by tfr_comp_idx;

                    
                    tfr_true_locs = h.cfg.source.vx_locs;
                    tfr_seed_locs = nan(size(FC_seed,2),3);
                    tfr_seed_locs(~isnan(seed_idx),:) = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(seed_idx(~isnan(seed_idx)),:);
                    tfr_comp_locs = nan(size(FC_comp,2),3);
                    tfr_comp_locs(~isnan(comp_idx),:) = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(comp_idx(~isnan(comp_idx)),:);
                    
                    vol =[];
                    ROI_freq_int = str2num(h.edit_3D_PLV_PLI_ROI_freq_int.String);
                    ROI_time_int = str2num(h.edit_3D_PLV_PLI_ROI_time_int.String);
                    freq_int = ROI_freq_int;
                    time_int = ROI_time_int;
                    freqs = h.inv_soln(h.current_inv_soln).TFR_results.TFR_freqs;
                    lat = h.cfg.study.lat_sim;
                    plot_flag = 0;
                    fc_scale = [];
                    fc_map_scale = [];
                    ROI_only_flag=1; % only calculating for ROI
                    end
                %% Step-wise thresholding for TFR ROI metrics
                if h.radio_3D_PLV_PLI_ROI_stepwise_stats.Value == 1 % calculate stepwise ROC stats for ROI
                    thresh_vals = [1e-10:str2num(h.edit_3D_PLV_PLI_ROI_step_size.String):1]*10; %10 dB
                    
                    [h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics] = calc_MCC_TFR_thresh(FC_true, FC_seed, FC_comp,...
                        freq_int, time_int, ROI_freq_int, ROI_time_int, thresh_vals, freqs, lat, plot_flag, ROI_only_flag);
                    % Fitting FPR & TPR data to create fitted ROC curve
                    [xData, yData] = prepareCurveData( squeeze([h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.ROI_perf.FPR]),...
                        squeeze([h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.ROI_perf.TPR]) );
                    if length(xData)>1 % data exists
                        ft = fittype( 'smoothingspline' ); opts = fitoptions( 'Method', 'SmoothingSpline' ); opts.Normalize = 'on'; opts.SmoothingParam = 0.9999;
                        % Fit model to data.
                        [fitresult, gof] = fit( xData, yData, ft, opts );
                        x_plv= 0:.01:1; y_plv = feval(fitresult,x_plv); y_plv(y_plv>1)=1;
                        h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR = x_plv;
                        h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.fitted_ROI_TPR = y_plv;
                        h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.fitted_ROI_AUC = trapz(x_plv,y_plv);
                    else
                        h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR = [];
                        h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.fitted_ROI_TPR = [];
                        h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.fitted_ROI_AUC = [];
                    end
                end
            end
        end
        %% %%%%% TFR Induced Power
        if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'avg_seed_wt_ind')
            if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_ind)
                %% prepare TFR data for ROI metrics
                if h.radio_3D_PLV_PLI_ROI_stepwise_stats.Value == 1 || h.radio_3D_PLV_PLI_ROI_surg_std_stats.Value == 1 % calculate stepwise ROC stats for ROI
                    %% getting seedPLV idx only = seed_plv_idx
                    true_idx = find(~isnan(h.inv_soln(h.current_inv_soln).plv_seed_idx));     % matching up with FC_true indices so the true PLV source1v2 matches with hit seed source1v2
                    
                    FC_true = permute(h.cfg.source.TFR_results.avg_true_wt_ind,[1 3 2]);
                    FC_seed = zeros(size(FC_true));  % setting all missed source PLV to zeros
%                     FC_seed(:,true_idx,:) = permute(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_ind(:,:,true_idx),[1 3 2]);
                    FC_seed(:,true_idx,:) = permute(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_ind(:,:,true_idx),[1 3 2]);
                    FC_comp = permute(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_ind,[1 3 2]);
                    tfr_comp_indices = comp_idx; % leadfield indices ordered by tfr_comp_idx;
                    
                    tfr_true_locs = h.cfg.source.vx_locs;
                    tfr_seed_locs = nan(size(FC_seed,2),3);
                    tfr_seed_locs(~isnan(seed_idx),:) = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(seed_idx(~isnan(seed_idx)),:);
                    tfr_seed_locs = nan(size(FC_comp,2),3);
                    tfr_comp_locs(~isnan(comp_idx),:) = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(comp_idx(~isnan(comp_idx)),:);
                    
                    vol =[];
                    ROI_freq_int = str2num(h.edit_3D_PLV_PLI_ROI_freq_int.String);
                    ROI_time_int = str2num(h.edit_3D_PLV_PLI_ROI_time_int.String);
                    freq_int = ROI_freq_int;
                    time_int = ROI_time_int;
                    freqs = h.inv_soln(h.current_inv_soln).TFR_results.TFR_freqs;
                    lat = h.cfg.study.lat_sim;
                    plot_flag = 0;
                    fc_scale = [];
                    fc_map_scale = [];
                    ROI_only_flag=1; % only calculating for ROI
                end
                %% Step-wise thresholding for TFR ROI metrics
                if h.radio_3D_PLV_PLI_ROI_stepwise_stats.Value == 1 % calculate stepwise ROC stats for ROI
                    thresh_vals = [1e-10:str2num(h.edit_3D_PLV_PLI_ROI_step_size.String):1]*10; %10 dB
                    
                    [h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics] = calc_MCC_TFR_thresh(FC_true, FC_seed, FC_comp,...
                        freq_int, time_int, ROI_freq_int, ROI_time_int, thresh_vals, freqs, lat, plot_flag, ROI_only_flag);
                    % Fitting FPR & TPR data to create fitted ROC curve
                    [xData, yData] = prepareCurveData( squeeze([h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.ROI_perf.FPR]),...
                        squeeze([h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.ROI_perf.TPR]) );
                    if length(xData)>1 % data exists
                        ft = fittype( 'smoothingspline' ); opts = fitoptions( 'Method', 'SmoothingSpline' ); opts.Normalize = 'on'; opts.SmoothingParam = 0.9999;
                        % Fit model to data.
                        [fitresult, gof] = fit( xData, yData, ft, opts );
                        x_plv= 0:.01:1; y_plv = feval(fitresult,x_plv); y_plv(y_plv>1)=1;
                        h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR = x_plv;
                        h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.fitted_ROI_TPR = y_plv;
                        h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.fitted_ROI_AUC = trapz(x_plv,y_plv);
                    else
                        h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR = [];
                        h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.fitted_ROI_TPR = [];
                        h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.fitted_ROI_AUC = [];
                    end
                end
            end
        end
        
        %% %%%%% PLV
        if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'plv_based')
            if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.plv_based)
                
                %% prepare PLV data for ROI metrics
                if h.radio_3D_PLV_PLI_ROI_stepwise_stats.Value == 1 || h.radio_3D_PLV_PLI_ROI_surg_std_stats.Value == 1 % calculate stepwise ROC stats for ROI
                    %% getting seedPLV idx only = seed_plv_idx
                    clear xv xt;
                    tidx = nchoosek(seed_idx,2);
                    xidx = perms(seed_idx); xidx = xidx(:,1:2);
                    for v=1:size(tidx,1)
                        if ~isempty(find(ismember(xidx,tidx(v,:),"rows")))
                            xt(v) = find(ismember(xidx,tidx(v,:),"rows"));
                        else
                            xt(v) = nan;
                        end
                    end % find true plv idx
                    ridx = setdiff(1:length(xidx),xt);   % PLV to be removed from calcuations because they are repeats "1 2" and "2 1"
                    for v=1:length(xidx)
                        if ~isempty(find(ismember(plv_contrasts,xidx(v,:),"rows")))
                            xv(v) = find(ismember(plv_contrasts,xidx(v,:),"rows"));
                        else
                            xv(v) = nan;
                        end
                    end
                    xv = xv(~isnan(xv));
                    
                    seed_plv_idx = find(ismember(plv_contrasts,tidx,"rows"));
                    true_idx = sum(~isnan(tidx),2)==2;     % matching up with FC_true indices so the true PLV source1v2 matches with hit seed source1v2
                    
                    FC_true = h.cfg.source.TFR_results.plv_based;
                    FC_seed = zeros(size(FC_true));  % setting all missed source PLV to zeros
                    FC_seed(:,true_idx,:) = h.inv_soln(h.current_inv_soln).TFR_results.plv_based(:,seed_plv_idx,:);
                    FC_comp = h.inv_soln(h.current_inv_soln).TFR_results.plv_based; FC_comp(:,xv,:) = nan;
                    plv_comp_indices = h.inv_soln(h.current_inv_soln).plv_contrasts; % leadfield indices ordered by plv_comp_idx;
                    plv_comp_idx = plv_idx; % nchoosek indices
                    
                    plv_true_locs = h.cfg.source.vx_locs;
                    plv_seed_locs = nan(size(FC_seed,2),3);
                    plv_seed_locs(~isnan(seed_idx),:) = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(seed_idx(~isnan(seed_idx)),:);
                     plv_comp_locs = nan(size(FC_comp,2),3);
                   plv_comp_locs = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(comp_idx(~isnan(comp_idx)),:);
                    
                    vol =[];
                    ROI_freq_int = str2num(h.edit_3D_PLV_PLI_ROI_freq_int.String);
                    ROI_time_int = str2num(h.edit_3D_PLV_PLI_ROI_time_int.String);
                    freq_int = ROI_freq_int;
                    time_int = ROI_time_int;
                    freqs = h.inv_soln(h.current_inv_soln).TFR_results.TFR_freqs;
                    lat = h.cfg.study.lat_sim;
                    plot_flag = 0;
                    fc_scale = [];
                    fc_map_scale = [];
                    ROI_only_flag=1; % only calculating for ROI
                end
                %% Step-wise thresholding for ROI metrics
                if h.radio_3D_PLV_PLI_ROI_stepwise_stats.Value == 1 % calculate stepwise ROC stats for ROI
                    thresh_vals = 1e-10:str2num(h.edit_3D_PLV_PLI_ROI_step_size.String):1;
                    
                    [h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics] = calc_MCC_TFR_PLV_thresh(FC_true, FC_seed, FC_comp,...
                        plv_comp_indices, plv_comp_idx, plv_true_locs, plv_seed_locs, plv_comp_locs,...
                        vol, freq_int, time_int, ROI_freq_int, ROI_time_int, thresh_vals, freqs, lat,...
                        plot_flag, fc_scale, fc_map_scale, ROI_only_flag);
                    % Fitting FPR & TPR data to create fitted ROC curve
                    [xData, yData] = prepareCurveData( squeeze([h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.ROI_perf.FPR]),...
                        squeeze([h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.ROI_perf.TPR]) );
                    if length(xData)>1 % data exists
                        ft = fittype( 'smoothingspline' ); opts = fitoptions( 'Method', 'SmoothingSpline' ); opts.Normalize = 'on'; opts.SmoothingParam = 0.9999;
                        % Fit model to data.
                        [fitresult, gof] = fit( xData, yData, ft, opts );
                        x_plv= 0:.01:1; y_plv = feval(fitresult,x_plv); y_plv(y_plv>1)=1;
                        h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.fitted_ROI_FPR = x_plv;
                        h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.fitted_ROI_TPR = y_plv;
                        h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.fitted_ROI_AUC = trapz(x_plv,y_plv);
                    else
                        h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.fitted_ROI_FPR = [];
                        h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.fitted_ROI_TPR = [];
                        h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.fitted_ROI_AUC = [];
                    end
                end
                %% Surrogate StdDev thresholding for ROI metrics
                if h.radio_3D_PLV_PLI_ROI_surg_std_stats.Value == 1
                    if isfield(h.cfg.source.TFR_results,'Noise') && isfield(h.inv_soln(h.current_inv_soln).TFR_results,'Noise')% calculate stepwise ROC stats for ROI
                        thresh_vals = 1e-10:str2num(h.edit_3D_PLV_PLI_ROI_step_size.String):1;
                        FC_true_noise_surg = h.cfg.source.TFR_results.Noise.plv_surg_based_mean;
                        FC_true_noise_surg_std = h.cfg.source.TFR_results.Noise.plv_surg_based_std;
                        FC_comp_noise_surg = h.inv_soln(h.current_inv_soln).TFR_results.Noise.plv_surg_based_mean;
                        FC_comp_noise_surg_std = h.inv_soln(h.current_inv_soln).TFR_results.Noise.plv_surg_based_std;
                        
                        std_step_size = str2num(h.edit_3D_PLV_PLI_ROI_std_step_size.String);
                        std_stop = str2num(h.edit_3D_PLV_PLI_ROI_std_stop_size.String);
                        [h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics] = calc_MCC_TFR_PLV_surg_noise(...
                            FC_true, FC_true_noise_surg, FC_true_noise_surg_std, FC_seed, FC_comp, FC_comp_noise_surg, FC_comp_noise_surg_std,...
                            std_step_size, std_stop, plv_comp_indices, plv_comp_idx, plv_true_locs, plv_seed_locs, plv_comp_locs,...
                            vol, freq_int, time_int, ROI_freq_int, ROI_time_int, freqs, lat, plot_flag, fc_scale, fc_map_scale, ROI_only_flag);
                        
                        % Fitting FPR & TPR data to create fitted ROC curve
                         % Fitting FPR & TPR data to create fitted ROC curve
                        x = squeeze([h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.ROI_perf.FPR]); 
                        y = squeeze([h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.ROI_perf.TPR]); 
                        if max(x)<1 
                            x1=x; x1(2:end+1)=x; x1(1) = 1; x=x1; 
                            y1=y; y1(2:end+1)=y; y1(1) = y(1); y=y1;
                        end
                        [xData, yData] = prepareCurveData(x,y);
                        if length(xData)>1 % data exists
                            ft = fittype( 'smoothingspline' ); opts = fitoptions( 'Method', 'SmoothingSpline' ); opts.Normalize = 'on'; opts.SmoothingParam = 0.9999;
                            % Fit model to data.
                            [fitresult, gof] = fit( xData, yData, ft, opts );
                            x_plv= 0:.01:1; y_plv = feval(fitresult,x_plv); y_plv(y_plv>1)=1;
                            h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.fitted_ROI_FPR = x_plv;
                            h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.fitted_ROI_TPR = y_plv;
                            h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.fitted_ROI_AUC = trapz(x_plv,y_plv);
                        else
                            h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.fitted_ROI_FPR = [];
                            h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.fitted_ROI_TPR = [];
                            h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.fitted_ROI_AUC = [];
                        end
                    else
                        hm = warndlg('Please check box "include projected Noise FC"','Need to include Noise FC');
                        hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
                    end
                end
                
                
            else
                hm = warndlg('Please Calculate PLV Connectivity','PLV Connectivity Not Calculated');
                hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
            end
        else
            hm = warndlg('Please Calculate PLV Connectivity','PLV Connectivity Not Calculated');
            hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
        end
        %% %%%%% PLI
        if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'pli_based')
            if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.pli_based)
                %% prepare PLI data for ROI metrics
                if h.radio_3D_PLV_PLI_ROI_stepwise_stats.Value == 1 || h.radio_3D_PLV_PLI_ROI_surg_std_stats.Value == 1 % calculate stepwise ROC stats for ROI
                    %% getting seedPLV idx only = seed_plv_idx
                    clear xv xt;
                    tidx = nchoosek(seed_idx,2);
                    xidx = perms(seed_idx); xidx = xidx(:,1:2);
                    for v=1:size(tidx,1);
                        if ~isempty(find(ismember(xidx,tidx(v,:),"rows")))
                            xt(v) = find(ismember(xidx,tidx(v,:),"rows"));
                        else
                            xt(v) = nan;
                        end
                    end % find true plv idx
                    ridx = setdiff(1:length(xidx),xt);   % PLV to be removed from calcuations because they are repeats "1 2" and "2 1"
                    for v=1:length(xidx)
                        if ~isempty(find(ismember(plv_contrasts,xidx(v,:),"rows")))
                            xv(v) = find(ismember(plv_contrasts,xidx(v,:),"rows"));
                        else
                            xv(v) = nan;
                        end
                    end
                    xv = xv(~isnan(xv));
                    
                    seed_plv_idx = find(ismember(plv_contrasts,tidx,"rows"));
                    true_idx = sum(~isnan(tidx),2)==2;     % matching up with FC_true indices so the true PLV source1v2 matches with hit seed source1v2
                    
                    FC_true = h.cfg.source.TFR_results.pli_based;
                    FC_seed = zeros(size(FC_true));  % setting all missed source PLV to zeros
                    FC_seed(:,true_idx,:) = h.inv_soln(h.current_inv_soln).TFR_results.pli_based(:,seed_plv_idx,:);
                    FC_comp = h.inv_soln(h.current_inv_soln).TFR_results.pli_based; FC_comp(:,xv,:) = nan;
                    plv_comp_indices = h.inv_soln(h.current_inv_soln).plv_contrasts; % leadfield indices ordered by plv_comp_idx;
                    plv_comp_idx = plv_idx; % nchoosek indices
                    
                    plv_true_locs = h.cfg.source.vx_locs;
                    plv_seed_locs = nan(size(FC_seed,2),3);
                    plv_seed_locs(~isnan(seed_idx),:) = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(seed_idx(~isnan(seed_idx)),:);
                    plv_comp_locs = nan(size(FC_comp,2),3);
                    plv_comp_locs = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(comp_idx(~isnan(comp_idx)),:);
                    
                    vol =[];
                    ROI_freq_int = str2num(h.edit_3D_PLV_PLI_ROI_freq_int.String);
                    ROI_time_int = str2num(h.edit_3D_PLV_PLI_ROI_time_int.String);
                    freq_int = ROI_freq_int;
                    time_int = ROI_time_int;
                    freqs = h.inv_soln(h.current_inv_soln).TFR_results.TFR_freqs;
                    lat = h.inv_soln(h.current_inv_soln).TFR_results.pli_lat; % h.cfg.study.lat_sim;
                    plot_flag = 0;
                    fc_scale = [];
                    fc_map_scale = [];
                    ROI_only_flag=1; % only calculating for ROI
                end
                
                %% Step-wise thresholding for ROI metrics
                if h.radio_3D_PLV_PLI_ROI_stepwise_stats.Value == 1 % calculate stepwise ROC stats for ROI
                    thresh_vals = 1e-10:str2num(h.edit_3D_PLV_PLI_ROI_step_size.String):1;
                    
                    [h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics] = calc_MCC_TFR_PLV_thresh(FC_true, FC_seed, FC_comp,...
                        plv_comp_indices, plv_comp_idx, plv_true_locs, plv_seed_locs, plv_comp_locs,...
                        vol, freq_int, time_int, ROI_freq_int, ROI_time_int, thresh_vals, freqs, lat,...
                        plot_flag, fc_scale, fc_map_scale, ROI_only_flag);
                    % Fitting FPR & TPR data to create fitted ROC curve
                    [xData, yData] = prepareCurveData( squeeze([h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.ROI_perf.FPR]),...
                        squeeze([h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.ROI_perf.TPR]) );
                    if length(xData)>1 % data exists
                        ft = fittype( 'smoothingspline' ); opts = fitoptions( 'Method', 'SmoothingSpline' ); opts.Normalize = 'on'; opts.SmoothingParam = 0.9999;
                        % Fit model to data.
                        [fitresult, gof] = fit( xData, yData, ft, opts );
                        x_plv= 0:.01:1; y_plv = feval(fitresult,x_plv); y_plv(y_plv>1)=1;
                        h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.fitted_ROI_FPR = x_plv;
                        h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.fitted_ROI_TPR = y_plv;
                        h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.fitted_ROI_AUC = trapz(x_plv,y_plv);
                    else
                        h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.fitted_ROI_FPR = [];
                        h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.fitted_ROI_TPR = [];
                        h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.fitted_ROI_AUC = [];
                    end
                end
                
                %% Surrogate StdDev thresholding for ROI metrics
                if h.radio_3D_PLV_PLI_ROI_surg_std_stats.Value == 1
                    if isfield(h.cfg.source.TFR_results,'Noise') && isfield(h.inv_soln(h.current_inv_soln).TFR_results,'Noise')% calculate stepwise ROC stats for ROI
                        thresh_vals = 1e-10:str2num(h.edit_3D_PLV_PLI_ROI_step_size.String):1;
                        FC_true_noise_surg = h.cfg.source.TFR_results.Noise.pli_surg_based_mean;
                        FC_true_noise_surg_std = h.cfg.source.TFR_results.Noise.pli_surg_based_std;
                        FC_comp_noise_surg = h.inv_soln(h.current_inv_soln).TFR_results.Noise.pli_surg_based_mean;
                        FC_comp_noise_surg_std = h.inv_soln(h.current_inv_soln).TFR_results.Noise.pli_surg_based_std;
                        
                        std_step_size = str2num(h.edit_3D_PLV_PLI_ROI_std_step_size.String);
                        std_stop = str2num(h.edit_3D_PLV_PLI_ROI_std_stop_size.String);
                        [h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics] = calc_MCC_TFR_PLV_surg_noise(...
                            FC_true, FC_true_noise_surg, FC_true_noise_surg_std, FC_seed, FC_comp, FC_comp_noise_surg, FC_comp_noise_surg_std,...
                            std_step_size, std_stop, plv_comp_indices, plv_comp_idx, plv_true_locs, plv_seed_locs, plv_comp_locs,...
                            vol, freq_int, time_int, ROI_freq_int, ROI_time_int, freqs, lat, plot_flag, fc_scale, fc_map_scale, ROI_only_flag);
                        
                        % Fitting FPR & TPR data to create fitted ROC curve
                        x = squeeze([h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.ROI_perf.FPR]); 
                        y = squeeze([h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.ROI_perf.TPR]); 
                        if max(x)<1 
                            x1=x; x1(2:end+1)=x; x1(1) = 1; x=x1; 
                            y1=y; y1(2:end+1)=y; y1(1) = y(1); y=y1;
                        end
                        [xData, yData] = prepareCurveData( x,y );
                        if length(xData)>1 % data exists
                            ft = fittype( 'smoothingspline' ); opts = fitoptions( 'Method', 'SmoothingSpline' ); opts.Normalize = 'on'; opts.SmoothingParam = 0.9999;
                            % Fit model to data.
                            [fitresult, gof] = fit( xData, yData, ft, opts );
                            x_plv= 0:.01:1; y_plv = feval(fitresult,x_plv); y_plv(y_plv>1)=1;
                            h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.fitted_ROI_FPR = x_plv;
                            h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.fitted_ROI_TPR = y_plv;
                            h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.fitted_ROI_AUC = trapz(x_plv,y_plv);
                        else
                            h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.fitted_ROI_FPR = [];
                            h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.fitted_ROI_TPR = [];
                            h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.fitted_ROI_AUC = [];;
                        end
                    else
                        hm = warndlg('Please check box "include projected Noise FC"','Need to include Noise FC');
                        hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
                    end
                end
            else
                hm = warndlg('Please Calculate PLI Connectivity','PLI Connectivity Not Calculated');
                hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
            end
        else
            hm = warndlg('Please Calculate PLI Connectivity','PLI Connectivity Not Calculated');
            hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
        end
    else
        hm = warndlg('Please Calculate True Connectivity','True Connectivity Not Calculated');
        hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
    end
else
    hm = warndlg('Please Calculate True Connectivity','True Connectivity Not Calculated');
    hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
end
%%
try; 
    sm_plot_ROI_FC_metrics;
catch 
end


