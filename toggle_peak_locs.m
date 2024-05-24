function toggle_peak_locs(varargin)
global h

h.map3D_peak_ori = handle(h.map3D_peak_ori);
if h.radio_3D_peak_locs.Value == 1
    try set(h.map3D_peak_locs,'Visible','on');  catch; end 
    try set(h.map3D_peak_ori,'Visible','on'); catch; end
    try set(h.current_peak_selected,'Visible','on'); catch; end
elseif h.radio_3D_peak_locs.Value == 0
    try set(h.map3D_peak_locs,'Visible','off');catch; end 
    try set(h.map3D_peak_ori,'Visible','off'); catch; end
    try set(h.current_peak_selected,'Visible','off'); catch; end
end
bs_plot_peak_waves;
