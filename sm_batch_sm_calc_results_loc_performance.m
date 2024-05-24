function h = sm_batch_sm_calc_results_loc_performance(h)

h.calc_results_flag = 1; % do not update plots while scanning through calculations of results 
fprintf('Calculating Results: Step-wise thresholding for MCC, TPR, and FPR\n'); 

%% setting slider to 0 and max(P.img) to get full range for calulating 
h.inv_soln(h.current_inv_soln).soln.plot_min_max = [0 max(abs(h.inv_soln(h.current_inv_soln).soln.P.img))];
% update_image_thresh_txt; 

%%
num_thresh_vals = 50;  % >100 takes time to compute.
%% Step-wise (absolute) thresholding of image from 0 to 1.1*max(img)
thresh_vals = linspace(0,h.inv_soln(h.current_inv_soln).soln.plot_min_max(2)*1.1,num_thresh_vals);
nan_vals = nan(size(thresh_vals));
h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_MCC = nan_vals;
h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_TPR = nan_vals;
h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_FPR = nan_vals;

h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_hits = nan_vals;
h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_miss = nan_vals;
h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_fa = nan_vals;
h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_CR = nan_vals;
nan_vals = nan(num_thresh_vals,3);
h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_loc_error = nan_vals;
h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_ori_error_mse = nan_vals;
h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_act = nan_vals;
h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_ctrl = nan_vals;

h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals = thresh_vals;

if ~isempty(h.inv_soln(h.current_inv_soln).peak_voxels)
    img_vals = h.inv_soln(h.current_inv_soln).peak_voxels(:,4); % peak values
    
    fprintf('Percent completed: ');
    for t=1:length(thresh_vals)
        fprintf('%.f ',100*t/length(thresh_vals))
        h.current_inv_soln_show_peak_idx = find(h.inv_soln(h.current_inv_soln).peak_voxels(:,4)>thresh_vals(t));
        h.current_inv_soln_hide_peak_idx = find(h.inv_soln(h.current_inv_soln).peak_voxels(:,4)<=thresh_vals(t));  % peaks to be hidden when below threshold
        h.slider_3D_image_thresh.Value = thresh_vals(t);
        h.slider_3D_image_thresh.Max = nanmax(img_vals)*1.1;
        h.slider_3D_image_thresh.Min = 0;
        
        h = sm_batch_sm_search_for_hits(h,'slider thresh');
        h = sm_batch_bs_calc_errors_inv_soln(h);
        h = sm_batch_sm_calc_localizer_performance(h);
        
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_MCC(t) = h.inv_soln(h.current_inv_soln).classifier_performance.MCC;
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_TPR(t) = h.inv_soln(h.current_inv_soln).classifier_performance.TPR;
        h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_FPR(t) = h.inv_soln(h.current_inv_soln).classifier_performance.FPR;
        
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
    fprintf('\n');
    idx = find(h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_hits==0);
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_act(idx,:) = 100;    % setting all "miss" to 100% wave error
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_ctrl(idx,:) = 100;   % setting all "miss" to 100% wave error
else
    % updating values for Not Hits and all misses
    zero_vals = zeros(size(thresh_vals));
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_MCC = zero_vals;
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_TPR = zero_vals;
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_FPR = zero_vals;
    
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_hits = zero_vals;
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_miss = zero_vals+length(h.cfg.source.vx_idx);
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_fa = zero_vals;
    h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_CR = zero_vals+size(h.anatomy.leadfield.H,3)-length(h.cfg.source.vx_idx);
end

%% update plots allowed again 
h.calc_results_flag = 0; 
