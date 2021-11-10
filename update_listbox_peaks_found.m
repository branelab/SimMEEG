function update_listbox_peaks_found(varargin)
global h


peak_clr = h.map3D_peak_locs(1).CData;

peak_voxels = h.current_3D_peak_voxels;
peak_idx = h.current_3D_peak_idx;

h.listbox_peaks_found.UserData.clr_str=''; h.listbox_peaks_found.String=''; h.listbox_peaks_found.Value = 1;
pre = '<HTML><FONT color="'; post = '</FONT></HTML>';
for v=1:length(peak_idx)
    peak_name = sprintf('%.f @ %.1f %.1f %.1f mm (value=%.3f)',peak_idx(v),peak_voxels(v,1:4) );
    clr_str = reshape( dec2hex( round(peak_clr(v,:)*255), 2 )',1, 6);
    h.listbox_peaks_found.UserData.clr_str{v} = [pre clr_str '">' peak_name post];
end
h.listbox_peaks_found.String = h.listbox_peaks_found.UserData.clr_str; 

