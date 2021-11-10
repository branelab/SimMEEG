function update_image_thresh_txt(varargin)
global h

min_max = str2num(h.edit_3D_min_max.String);
pos_ratio = ( h.slider_3D_image_thresh.Value-h.slider_3D_image_thresh.Min) / range([h.slider_3D_image_thresh.Min h.slider_3D_image_thresh.Max]);
s_pos = h.slider_3D_image_thresh.Position;
h.slider_3D_image_thresh_text_val.Position(1) = s_pos(1) - h.slider_3D_image_thresh_text_val.Position(3) - .01;
h.slider_3D_image_thresh_text_val.Position(2) = (range(h.slider_3D_image_thresh.Position([2 4]))*pos_ratio*.625) + h.slider_3D_image_thresh.Position(2);
h.slider_3D_image_thresh_text_val.String = num2str(h.slider_3D_image_thresh.Value);
h.slider_3D_image_thresh_text_max.String = num2str(min_max(2));