function sm_change_sensors(varargin)
global h

switch varargin{end}
    case 'on'
        h.panel_sensor_change.Visible = 'on';
        uistack(h.panel_sensor_change,'top')
end

%% UITABLE of sensors (editable)
sens_pos = h.anatomy.sens.elecpos;
sens_label = h.anatomy.sens.label';
sens_type = h.anatomy.sens.chantype;
if ~isempty(sens_label)
for v=1:length(sens_label)
    data(v,:) = {sens_label{v} sens_pos(v,1) sens_pos(v,2) sens_pos(v,3) sens_type{v}};
end
%     h.sens_table = uitable(h.panel_sensor_change,'Data', data, 'ColumnName', {'Label' 'X' 'Y' 'Z'},'Units','normalized','Position',[.05 .05 .45 .9],'ColumnEditable',true(1,length(sens_pos)));
h.sens_table.Data = data;
end

end







