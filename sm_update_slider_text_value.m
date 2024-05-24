function sm_update_slider_text_value(~,~,h_slider,~)
% this function will update the slider text value associated with the selected slider 
switch h_slider.UserData
    case 'slider_3D_image_thresh'
        update_image_thresh_txt;
    case 'slider_topo_scale'
        sm_set_topo_caxis;
end
