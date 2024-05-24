function h = sm_batch_update_image_thresh_txt(h)

pos_ratio = ( h.slider_3D_image_thresh.Value-h.slider_3D_image_thresh.Min) / range([h.slider_3D_image_thresh.Min h.slider_3D_image_thresh.Max]);

if h.slider_3D_image_thresh.Value<.001
    txt = compose("%.1e", h.slider_3D_image_thresh.Value);  % scientific notation
else
    txt = compose("%.3f", h.slider_3D_image_thresh.Value);  % reg notation
end

min_max = [h.slider_3D_image_thresh.Min h.slider_3D_image_thresh.Max]; %str2num(h.edit_3D_min_max.String);

if min_max(1)<.001; txt_min = compose("%.1e",min_max(1)); else txt_min = compose("%.3f",min_max(1)); end
if min_max(2)<.001; txt_max = compose("%.1e",min_max(2)); else txt_max= compose("%.3f",min_max(2)); end
if min_max(1)==0; txt_min='0'; end
txt_minmax = [char(txt_min) ' ' char(txt_max)];

h.slider_3D_image_thresh_text_val.String = txt;
h.slider_3D_image_thresh_text_max.String = txt_max;
h.edit_3D_min_max.String = sprintf('%.3f %.3f', txt_minmax);

