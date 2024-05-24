function sm_view_map_slices(varargin) 
global h        

if h.radio_3D_peak_flag.Value == 0  % show peaks only
    sm_show_peaks(h.current_inv_soln_plot_map3Dpeaks_idx,'on');  % shows ('on') or hides ('off') selected peaks in 3D image
    set(h.func3D,'Visible','on');
    set(h.func3D_midline,'Visible','off');
    
elseif h.radio_3D_peak_flag.Value == 1 && ~isempty(h.func3D_midline) % plot midline only
    set(h.func3D,'Visible','off');
    set(h.func3D_midline,'Visible','on');
    
end