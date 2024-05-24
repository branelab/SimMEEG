function sm_calc_results_loc_performance(varargin)
global h

h.calc_results_flag = 1; % do not update plots while scanning through calculations of results 
fprintf('Calculating Results: Step-wise thresholding for MCC, TPR, and FPR\n'); 

%% setting slider to 0 and max(P.img) to get full range for calulating 
h.inv_soln(h.current_inv_soln).soln.plot_min_max = [0 max(abs(h.inv_soln(h.current_inv_soln).soln.P.img))];
update_image_thresh_txt; 


%%
if ~isempty(h.inv_soln(h.current_inv_soln).peak_voxels)
    img_vals = h.inv_soln(h.current_inv_soln).peak_voxels(:,4); % peak values
    if h.radio_3D_stepwise_thresh.Value ==1
        num_thresh_vals = 50;  % >100 takes time to compute.
    else
        num_thresh_vals = 4;  % >100 takes time to compute.
    end
    thresh_vals = linspace(0,h.inv_soln(h.current_inv_soln).soln.plot_min_max(2)*1.1,num_thresh_vals);
    
    
    %% Step-wise (absolute) thresholding of image from 0 to 1.1*max(img)
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_MCC = nan(size(thresh_vals));
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_TPR = nan(size(thresh_vals));
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_FRP = nan(size(thresh_vals));
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals = thresh_vals;
    % hw = waitbar(0,'Calculating Performance Results: Step-wise thresholding');
    
    
    for t=1:length(thresh_vals)
        %     waitbar(t/length(thresh_vals),hw,'Calculating Performance Results: Step-wise thresholding')
        h.current_inv_soln_show_peak_idx = find(h.inv_soln(h.current_inv_soln).peak_voxels(:,4)>thresh_vals(t));
        h.current_inv_soln_hide_peak_idx = find(h.inv_soln(h.current_inv_soln).peak_voxels(:,4)<=thresh_vals(t));  % peaks to be hidden when below threshold
        h.slider_3D_image_thresh.Value = thresh_vals(t);
        h.slider_3D_image_thresh.Max = nanmax(img_vals)*1.1;
        h.slider_3D_image_thresh.Min = 0;
        
        sm_search_for_hits('slider thresh');
        sm_calc_localizer_performance();
        bs_calc_errors_inv_soln();
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_MCC(t) = h.inv_soln(h.current_inv_soln).classifier_performance.MCC;
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_TPR(t) = h.inv_soln(h.current_inv_soln).classifier_performance.TPR;
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_FPR(t) = h.inv_soln(h.current_inv_soln).classifier_performance.FPR;
        
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_hits(t).idx = h.inv_soln(h.current_inv_soln).classifier_metrics.Hits;
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_miss(t).idx = h.inv_soln(h.current_inv_soln).classifier_metrics.Hits;
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_fa(t).idx = h.inv_soln(h.current_inv_soln).classifier_metrics.FA;
        
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_hits(t) = sum(~isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.Hits));
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_miss(t) = sum(isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.Hits));
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_fa(t) = sum(~isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.FA));
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_CR(t) = h.inv_soln(h.current_inv_soln).classifier_metrics.num_CR;
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_loc_error(t,:) = h.inv_soln(h.current_inv_soln).classifier_metrics.loc_error;
        if isempty(h.inv_soln(h.current_inv_soln).ori_error_mse)   % no sources found
            h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_ori_error_mse(t,:) = nan(size(h.sim_data.cfg.source.vx_idx));
            h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_act(t,:) = nan(size(h.sim_data.cfg.source.vx_idx));
            h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_ctrl(t,:) = nan(size(h.sim_data.cfg.source.vx_idx));
        else
            h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_ori_error_mse(t,:) = h.inv_soln(h.current_inv_soln).ori_error_mse;
            h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_act(t,:) = h.inv_soln(h.current_inv_soln).wave_error_mse_norm_act;
            h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_ctrl(t,:) = h.inv_soln(h.current_inv_soln).wave_error_mse_norm_ctrl;
        end
        
    end
    % close(hw);
    idx = find(h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_hits==0);
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_act(idx,:) = 100;    % setting all "miss" to 100% wave error
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_ctrl(idx,:) = 100;   % setting all "miss" to 100% wave error
    
    %% Standard Deviation (std) thresholding of image from 0 to 3 stdev
    % std_vals = linspace(0,10,num_thresh_vals);
    % h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_MCC = nan(size(std_vals));
    % h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_TPR = nan(size(std_vals));
    % h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_FRP = nan(size(std_vals));
    % h.inv_soln(h.current_inv_soln).classifier_results.std_thresh_vals = std_vals;
    % h.inv_soln(h.current_inv_soln).classifier_results.stddev_img = std(h.inv_soln(h.current_inv_soln).soln.P.img);
    %
    % for t=1:length(std_vals)
    %     null_thresh = h.inv_soln(h.current_inv_soln).classifier_results.stddev_img*std_vals(t);
    %     h.current_inv_soln_show_peak_idx = find( h.inv_soln(h.current_inv_soln).peak_voxels(:,4) > null_thresh );
    %     h.current_inv_soln_hide_peak_idx = find( h.inv_soln(h.current_inv_soln).peak_voxels(:,4) <= null_thresh );  % peaks to be hidden when below threshold
    %     sm_search_for_hits('slider thresh');
    %     sm_calc_localizer_performance();
    %     h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_MCC(t) = h.inv_soln(h.current_inv_soln).classifier_performance.MCC;
    %     h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_TPR(t) = h.inv_soln(h.current_inv_soln).classifier_performance.TPR;
    %     h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_FPR(t) = h.inv_soln(h.current_inv_soln).classifier_performance.FPR;
    %
    %     h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_num_hits(t) = sum(~isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.Hits));
    %     h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_num_miss(t) = sum(isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.Hits));
    %     h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_num_fa(t) = sum(~isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.FA));
    %     h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_num_CR(t) = h.inv_soln(h.current_inv_soln).classifier_metrics.num_CR;
    %
    % end
    
    
    
    %% Plotting results
    sm_plot_results_loc_performance();
    
    %% update plots allowed again
    h.calc_results_flag = 0;
    
    %% resetting image and slider values back to before this functiona call was made
    set_3D_image_thresh();
    
end
