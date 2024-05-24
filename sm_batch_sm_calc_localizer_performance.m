function h = sm_batch_sm_calc_localizer_performance(h)
%function [metrics]=Sim_study_metrics_v2(eeg_sim,ssp,true_locs,true_idx,inv_soln,HeadModelMat,inv_type,act_samps,inv_data,ssp_data);
% This program will calculate the SimStudy metrics to evaluate the LCMV & MCMV performances.
%
metrics.true_idx = h.sim_data.cfg.source.vx_idx;
metrics.dist_thresh = (sum(3*(str2num(h.edit_inv_dist_search_thresh.String)).^2)).^.5; % scalar distances = radius

%% creating an empty array if no inv_soln sources exist
true_locs = h.sim_data.cfg.source.vx_locs;  % true source locations

%     h.current_inv_soln_show_peak_idx

if ~isempty(h.current_inv_soln_show_peak_idx) % no voxels found
    
    
    %% localizer performance  = hits / (FA+Miss)
    peak_locs = nan(size(true_locs));
    lf_idx = h.current_peak_hit_lf_idx; %h.current_peak_hit_lf_idx(~isnan(h.current_peak_hit_lf_idx));
    hit_idx = find(~isnan(lf_idx)==1);
    miss_idx = find(isnan(lf_idx)==1);
    if ~isempty(lf_idx)
        peak_locs(hit_idx,:) = h.current_3D_peak_voxels(hit_idx,1:3); %h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(lf_idx(hit_idx),:);
    end
    %% first criteria find peaks within dist_thresh surrounding true_locs
    metrics.loc_error = nan(size(peak_locs,1),size(true_locs,1));
    xd = true_locs - peak_locs;
    metrics.loc_error = (nanmean(xd.^2,2)).^.5;
    metrics.loc_error(miss_idx) = nan;
    
else
    peak_voxels =[]; peak_idx =[];
    peak_locs = [];
    metrics.loc_error = nan(size(metrics.true_idx));
    h.current_peak_hit_lf_idx = nan(size(metrics.true_idx));
    h.current_peak_miss_lf_idx = metrics.true_idx(isnan(h.current_peak_hit_lf_idx));    % Missed true locations
    h.current_peak_fa_lf_idx = [];
end

metrics.Hits = h.current_peak_hit_lf_idx;
metrics.FA = h.current_peak_fa_lf_idx;
metrics.Miss = h.current_peak_miss_lf_idx;
metrics.num_CR = size(h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,1) - ( sum(~isnan(metrics.Hits)) + sum(~isnan(metrics.FA)) + sum(~isnan(metrics.Miss)) );

h.inv_soln(h.current_inv_soln).classifier_metrics = metrics;

% inv_soln's classifier performance
[perf]=sm_calc_classifier_performance(sum(~isnan(metrics.Hits)),sum(~isnan(metrics.Miss)),metrics.num_CR,sum(~isnan(metrics.FA)));
if isnan(perf.MCC); perf.MCC=0; end
h.inv_soln(h.current_inv_soln).classifier_performance = perf;



