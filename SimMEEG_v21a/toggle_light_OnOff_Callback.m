function toggle_light_OnOff_Callback(varargin)
global h

axes(h.axes_3D_images); 

if h.toggle_light_OnOff.Value == 0    % off
    % lighting off;
    delete(findall(gca,'Type','light'));
    h.toggle_light_OnOff.ForegroundColor = [1 0 0];  h.toggle_light_OnOff.BackgroundColor = [1 1 1]*.9;
    h.toggle_light_OnOff.String = 'Light Off';
    
elseif h.toggle_light_OnOff.Value == 1
      delete(findall(gca,'Type','light'))
%     if ~isvalid(h.camlight)
        h.camlight1 = camlight;
        h.camlight2 = camlight; % double camlight to make it brighter
        lighting gouraud; material dull;
%     end
    h.toggle_light_OnOff.ForegroundColor = [0 .6 0];  h.toggle_light_OnOff.BackgroundColor = [1 1 1];
    h.toggle_light_OnOff.String = 'Light On';
end
