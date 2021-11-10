function toggle_peak_locs(varargin)
global h

h.map3D_peak_ori = handle(h.map3D_peak_ori);
if h.radio_3D_peak_locs.Value == 1
    h.map3D_peak_locs(1).Visible = 'on'; %h.map3D_peak_locs(2).Visible = 'on';
    h.map3D_peak_locs(2).Visible = 'on'; %h.map3D_peak_locs(2).Visible = 'on';
    for v=1:length(h.map3D_peak_ori); h.map3D_peak_ori(v).Visible='on'; end
elseif h.radio_3D_peak_locs.Value == 0
    h.map3D_peak_locs(1).Visible = 'off'; %h.map3D_peak_locs(2).Visible = 'off';
    h.map3D_peak_locs(2).Visible = 'off'; %h.map3D_peak_locs(2).Visible = 'off';
    for v=1:length(h.map3D_peak_ori); h.map3D_peak_ori(v).Visible='off'; end
end
bs_plot_peak_waves;
