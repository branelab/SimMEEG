function update_listbox_peaks_found(varargin)
global h

peak_clr = h.ln_clr;

if ~isempty(h.current_3D_peak_voxels)
    peak_voxels = nan(size(h.sim_data.cfg.source.vx_locs,1),5);
    peak_voxels = h.current_3D_peak_voxels; %
else
    peak_voxels = nan(size(h.sim_data.cfg.source.vx_locs,1),5);
end

%% hiding False Alarms that are below slider threshold

% nan_idx = find(isnan(peak_voxels(:,1))); nan_idx = find(ismember(1:length(h.sim_data.cfg.source.vx_idx),nan_idx));
% good_idx = setdiff(1:size(peak_voxels,1),nan_idx);
% peak_voxels = peak_voxels(good_idx,:);


%% Updating listbox
h.listbox_peaks_found.UserData.clr_str=''; h.listbox_peaks_found.String=''; 
if h.listbox_peaks_found.Value > length(h.listbox_peaks_found.String)
h.listbox_peaks_found.Value = 1;
end


pre = '<HTML><FONT color="'; post = '</FONT></HTML>';
for v=1:size(peak_voxels,1)
    if isnan(peak_voxels(v,4))
        peak_name = 'Miss';
    else
        if peak_voxels(v,4)<.001    
            val = compose("%.1e", peak_voxels(v,4));    % scientific notation
        else
             val = compose("%.3f", peak_voxels(v,4));    % reg notation
       end
        peak_name = sprintf('%.f @ %.1f %.1f %.1f mm (value=%s)',peak_voxels(v,5),peak_voxels(v,1:3),val);
        
    end
    clr_str = reshape( dec2hex( round(peak_clr(v,:)*255), 2 )',1, 6);
    h.listbox_peaks_found.UserData.clr_str{v} = [pre clr_str '">' peak_name post];
end
h.listbox_peaks_found.String = h.listbox_peaks_found.UserData.clr_str;

