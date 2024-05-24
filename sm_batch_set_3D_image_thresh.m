function h = sm_batch_set_3D_image_thresh(h)

min_max = str2num(h.edit_3D_min_max.String);

h.current_3D_thresh = h.slider_3D_image_thresh.Value;
h.axes_3D_images.CLim = min_max;

% updating inv_soln min_max
try
    h.inv_soln(h.current_inv_soln).soln.plot_min_max = min_max;
    h.inv_soln(h.current_inv_soln).soln.plot_thresh = h.slider_3D_image_thresh.Value;
catch
end

%% Update --> "h.current_3D_peak_idx" and "h.current_3D_peak_voxels" based on slider threshold
% resetting to original peak_voxels when threshold=0; during run_source_modeling.m
h.current_3D_peak_voxels = h.inv_soln(h.current_inv_soln).peak_voxels;
h.current_3D_peak_idx = h.inv_soln(h.current_inv_soln).peak_idx;

if ~isempty(h.current_3D_peak_voxels)
% Thresholding base on slider threshold
vox_num = h.current_3D_peak_voxels( (h.current_3D_peak_voxels(:,4)>h.current_3D_thresh),5); % leadfield voxel number of source found
h.current_inv_soln_show_peak_idx = []; h.current_inv_soln_plot_map3Dpeaks_idx = [];
for v=1:length(vox_num)
    h.current_inv_soln_show_peak_idx(v) = find(h.inv_soln(h.current_inv_soln).peak_voxels(:,5)==vox_num(v));
    h.current_inv_soln_plot_map3Dpeaks_idx(v) = find(h.inv_soln(h.current_inv_soln).peak_voxels(:,5)==vox_num(v));   %
end
h.current_inv_soln_hide_peak_idx = setdiff(1:size(h.inv_soln(h.current_inv_soln).peak_voxels),h.current_inv_soln_show_peak_idx);  % map3d_locs peaks to be hidden when below threshold
h.current_3D_peak_idx(h.current_inv_soln_hide_peak_idx) = nan;

h = sm_batch_sm_search_for_hits(h,'slider thresh');

if size(h.current_3D_peak_voxels,1) < length(h.sim_data.cfg.source.vx_idx)
    peak_voxels = nan(length(h.sim_data.cfg.source.vx_idx),5);
    peak_voxels(~isnan(h.current_3D_peak_idx),:) = h.current_3D_peak_voxels(~isnan(h.current_3D_peak_idx),:);
    h.current_3D_peak_voxels = peak_voxels;
    h.current_3D_peak_idx = peak_voxels(:,5)';
end

h = sm_batch_bs_calc_errors_inv_soln(h);
h = sm_batch_sm_calc_localizer_performance(h);
else
    % no peak voxels found - no hits 
    h.current_peak_hit_lf_idx = [];
    h.current_peak_miss_lf_idx = h.cfg.source.vx_idx; % missed all sources
    h.current_peak_fa_lf_idx = [];
    h = sm_batch_bs_calc_errors_inv_soln(h);
    h = sm_batch_sm_calc_localizer_performance(h);
    
end



