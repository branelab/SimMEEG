function sm_open_ARM_params(varargin)
global h

if h.btn_ARM_open_panel.Value == 1
    h.panel_ARM_params.Visible = 'on';
else
    h.panel_ARM_params.Visible = 'off';
    if h.btn_ARM_change_source_params.Value == 1 % turn off the source params panel as well
    h.btn_ARM_change_source_params.Value = 0; 
    h.panel_ARM_source_params.Visible = 'off'; 
    end
    
end

