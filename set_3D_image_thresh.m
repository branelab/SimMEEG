function set_3D_image_thresh(varargin)
global h

update_image_thresh_txt;
min_max = str2num(h.edit_3D_min_max.String);
if min_max(1)<min_max(2)
    h.current_3D_thresh = h.slider_3D_image_thresh.Value;
    h.axes_3D_images.CLim = min_max;
else
    h.current_3D_thresh = 0;
    h.axes_3D_images.CLim = [0 2];
    h.edit_3D_min_max.String = '0 2.0';
end

% updating inv_soln min_max
try
    h.inv_soln(h.current_inv_soln).soln.plot_min_max = min_max;
    h.inv_soln(h.current_inv_soln).soln.plot_thresh = h.slider_3D_image_thresh.Value;
catch
end

set_3D_transparency; 


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

sm_search_for_hits('slider thresh');
if size(h.current_3D_peak_voxels,1) < length(h.sim_data.cfg.source.vx_idx)
    peak_voxels = nan(length(h.sim_data.cfg.source.vx_idx),5);
    peak_voxels(~isnan(h.current_3D_peak_idx),:) = h.current_3D_peak_voxels(~isnan(h.current_3D_peak_idx),:);
    h.current_3D_peak_voxels = peak_voxels;
    h.current_3D_peak_idx = peak_voxels(:,5)';
end




%% replot map3D_locs
sm_plot_3Dmap_locs(); 
update_listbox_peaks_found(); 

toggle_true_locs;
%% Show/Hide peaks below threshold
% if ~isempty(h.current_inv_soln_plot_map3Dpeaks_idx)   % Hiding peaks and slices - need to do first because some peaks share slices
%     sm_show_peaks(h.current_inv_soln_plot_map3Dpeaks_idx,'off');  % shows ('on') or hides ('off') selected peaks in 3D image
% end
if ~isempty(h.current_inv_soln_plot_map3Dpeaks_idx)      % showing peaks and slices
    sm_show_peaks(h.current_inv_soln_plot_map3Dpeaks_idx,'on');  % shows ('on') or hides ('off') selected peaks in 3D image
end

%% threshold maps
for s=1:size(h.func3D_midline,1)
    if isempty(h.func3D_midline(s).UserData) % stor original slice data;
%         h.func3D_midline(s).CData(isnan(h.func3D_midline(s).CData))=0;
    h.func3D_midline(s).UserData = h.func3D_midline(s).CData; 
    end
    h.func3D_midline(s).CData = h.func3D_midline(s).UserData; 
    h.func3D_midline(s).CData(h.func3D_midline(s).UserData<h.slider_3D_image_thresh.Value) = nan;
end
for s=1:size(h.func3D,1)
    if isempty(h.func3D(s).UserData) % stor original slice data;
%         h.func3D(s).CData(isnan(h.func3D(s).CData))=0;
    h.func3D(s).UserData = h.func3D(s).CData; 
    end
    h.func3D(s).CData = h.func3D(s).UserData; 
    h.func3D(s).CData(h.func3D(s).UserData<h.slider_3D_image_thresh.Value) = nan;
end



%% visualize
sm_view_map_slices(); 
bs_calc_errors_inv_soln();
sm_calc_localizer_performance(); 
update_listbox_peaks_found(); 
listbox_peaks_found_Callback();
toggle_true_locs(); toggle_peak_locs(); 


end
