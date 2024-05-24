function sm_invert_wave(varargin)
global h

for v=1:length(h.listbox_peaks_found.Value)
    h.current_inv_swf_plots(h.listbox_peaks_found.Value(v)).YData = h.current_inv_swf_plots(h.listbox_peaks_found.Value(v)).YData*-1;
    h.current_norm_peak_swf(:,h.listbox_peaks_found.Value(v)) = h.current_norm_peak_swf(:,h.listbox_peaks_found.Value(v))*-1; 
    vx_idx = h.current_3D_peak_voxels(h.listbox_peaks_found.Value(v),5);
    h.inv_soln(h.current_inv_soln).soln.ori(vx_idx,:) = h.inv_soln(h.current_inv_soln).soln.ori(vx_idx,:)*-1;
end

sm_calc_localizer_performance();
bs_calc_errors_inv_soln();

