function set_3D_image_thresh(varargin)
global h

update_image_thresh_txt;
min_max = str2num(h.edit_3D_min_max.String);

h.current_3D_thresh = h.slider_3D_image_thresh.Value;
h.axes_3D_images.CLim = min_max;

% updating inv_soln min_max
try
    h.inv_soln(h.current_inv_soln).soln.plot_min_max = min_max;
    h.inv_soln(h.current_inv_soln).soln.plot_thresh = h.slider_3D_image_thresh.Value;
catch
end

set_3D_transparency; 
bs_plot_inv_soln;

toggle_true_locs;

