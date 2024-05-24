function sm_set_topo_caxis(varargin)
global h
y_scale = str2num(h.edit_yscale.String);

s_pos = h.slider_topo_scale.Position;
h.slider_topo_scale_text_val.Position(1) = s_pos(1) - h.slider_topo_scale_text_val.Position(3) - .01;
h.slider_topo_scale_text_val.Position(2) = s_pos(2)+(s_pos(4)*h.slider_topo_scale.Value*.7)+(s_pos(4)*.1);

if h.menu_sens_type.Value==1 % MEG
    h.slider_topo_scale_text_val.String = [sprintf('%.f ',h.slider_topo_scale.Value*y_scale(2)) 'fT'];
elseif h.menu_sens_type.Value==2 % EEG
    h.slider_topo_scale_text_val.String = [sprintf('%.f ',h.slider_topo_scale.Value*y_scale(2)) char(181) 'V'];
end


h.topo_scale = y_scale*h.slider_topo_scale.Value;
h.axes_anatomy.CLim = h.topo_scale;