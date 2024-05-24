function toggle_light_OnOff_Callback(varargin)
global h

% axes(h.axes_3D_images); 

if h.toggle_light_OnOff.Value == 0    % off
    % lighting off;
    delete(findall(h.axes_3D_images,'Type','light'));
    h.toggle_light_OnOff.ForegroundColor = [1 0 0];  h.toggle_light_OnOff.BackgroundColor = [1 1 1]*.9;
    h.toggle_light_OnOff.String = 'Light Off';
    
elseif h.toggle_light_OnOff.Value == 1
      delete(findall(h.axes_3D_images,'Type','light'))
        b_val = .95; lgt_type = 'infinite';
        h.camlight1 = light(h.axes_3D_images, 'Color',[1 1 1]*b_val, 'Style',lgt_type,'Position',[1 0 0]); 
        h.camlight2 = light(h.axes_3D_images, 'Color',[1 1 1]*b_val, 'Style',lgt_type,'Position',[-1 0 0]); 
        h.camlight3 = light(h.axes_3D_images, 'Color',[1 1 1]*b_val, 'Style',lgt_type,'Position',[0 1 0]); 
        h.camlight4 = light(h.axes_3D_images, 'Color',[1 1 1]*b_val, 'Style',lgt_type,'Position',[0 -1 0]); 
        h.camlight5 = light(h.axes_3D_images, 'Color',[1 1 1]*b_val, 'Style',lgt_type,'Position',[0 0 1]); 
        h.camlight6 = light(h.axes_3D_images, 'Color',[1 1 1]*b_val, 'Style',lgt_type,'Position',[0 0 -1]); 

        
        lighting(h.axes_3D_images, 'phong'); material(h.axes_3D_images, 'dull');
    h.toggle_light_OnOff.ForegroundColor = [0 .6 0];  h.toggle_light_OnOff.BackgroundColor = [1 1 1];
    h.toggle_light_OnOff.String = 'Light On';
end
