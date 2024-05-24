function toggle_sens_waves(varargin)
global h

lg.String = h.axes_source_waves.Legend.String; 
if h.radio_3D_plot_erp_waves.Value == 1
    sdata = squeeze(nanmean(h.sim_data.sens_final,3));
    sdata = sdata/max(max(abs(sdata))) * max(abs(h.axes_source_waves.YLim));
    
    if isfield(h,'current_3D_sens_waves') % updating
       delete(h.current_3D_sens_waves)
        h.current_3D_sens_waves = plot(h.axes_source_waves, h.sim_data.cfg.study.lat_sim, sdata,'color',[1 1 1]*.8);
    else
        h.current_3D_sens_waves = plot(h.axes_source_waves, h.sim_data.cfg.study.lat_sim, sdata,'color',[1 1 1]*.8);
    end
else
    if isvalid(h.current_3D_sens_waves)
    set(h.current_3D_sens_waves,'Visible','off');
    end
end
h.axes_source_waves.Children = flipud(h.axes_source_waves.Children); 
% h.axes_source_waves.Legend.String = lg.String; h.axes_source_waves.Legend.String = lg.String;
legend(h.axes_source_waves, lg.String); 
