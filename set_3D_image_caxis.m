function set_3D_image_caxis(varargin)
global h
min_max = str2num(h.edit_3D_min_max.String);
h.slider_3D_image_thresh.Min = min_max(1); h.slider_3D_image_thresh.Max = min_max(2);
h.slider_3D_image_thresh_text_max.String = num2str(min_max(2));

if h.slider_3D_image_thresh.Value<min_max(1)
    h.slider_3D_image_thresh.Value = min_max(1);
elseif h.slider_3D_image_thresh.Value>min_max(2)
    h.slider_3D_image_thresh.Value = min_max(2);
end

h.slider_3D_image_thresh_text_max.String = num2str(min_max(2));
h.axes_3D_images.CLim = min_max;
h.inv_soln(h.current_inv_soln).soln.plot_min_max = min_max;
h.inv_soln(h.current_inv_soln).soln.plot_thresh = h.slider_3D_image_thresh.Value;

drawnow;
% set_3D_image_thresh; 
