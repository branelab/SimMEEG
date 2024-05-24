function toggle_true_locs(varargin)
global h
for v=1:length(h.cfg.source.vx_idx)
    if h.radio_3D_true_locs.Value == 1
        set(handle(h.map3D_true_locs(1,v)),'Visible','on'); set(handle(h.map3D_true_locs(2,v)),'Visible','on');
    elseif h.radio_3D_true_locs.Value == 0
        set(handle(h.map3D_true_locs(1,v)),'Visible','off'); set(handle(h.map3D_true_locs(2,v)),'Visible','off');
    end
end
bs_plot_peak_waves;