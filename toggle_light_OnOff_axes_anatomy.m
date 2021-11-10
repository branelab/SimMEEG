function toggle_light_OnOff_axes_anatomy(varargin)
global h

axes(h.axes_anatomy); 

if h.toggle_light_OnOff_anat.Value == 0    % off
    % lighting off;
    delete(findall(gca,'Type','light'));
    h.toggle_light_OnOff_anat.ForegroundColor = [1 0 0];  h.toggle_light_OnOff_anat.BackgroundColor = [1 1 1]*.9;
    h.toggle_light_OnOff_anat.String = 'Light Off';
    
elseif h.toggle_light_OnOff_anat.Value == 1
      delete(findall(gca,'Type','light'))
%     if ~isvalid(h.camlight)
        h.camlight1 = camlight;
        h.camlight2 = camlight; % double camlight to make it brighter
        lighting gouraud; material dull;
%     end
    h.toggle_light_OnOff_anat.ForegroundColor = [0 .6 0];  h.toggle_light_OnOff_anat.BackgroundColor = [1 1 1];
    h.toggle_light_OnOff_anat.String = 'Light On';
end
