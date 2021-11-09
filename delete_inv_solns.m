
function delete_inv_solns(varargin)
global h

h.current_inv_soln = h.listbox_inv_solns.Value;
keep_idx = setdiff(1:length(h.inv_soln),h.current_inv_soln);
h.inv_soln = h.inv_soln(keep_idx);
bs_update_3D_listbox;
h.listbox_inv_solns.Value=length(h.inv_soln);
h.current_3D_peak_idx = []; h.current_3D_peak_voxels = [];
try
    update_listbox_peaks_found
    h.axes_3D_images.clo;
    h.axes_source_waves.clo;
    h.axes_source_fft.clo;
    h.axes_invSoln_halfmax_width.clo;
    h.axes_invSoln_errors_locs.clo; h.axes_invSoln_errors_ori.clo; h.axes_invSoln_errors_waves.clo;
catch
end

