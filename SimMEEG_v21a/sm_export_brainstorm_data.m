function sm_export_brainstorm_data(varargin)
global h

% Francois to write code here for exporting source data and sensor data to brainstorm

if isappdata(0, 'BrainstormRunning')
    fprintf('Exported data to Brainstorm\n');
    
    
else
    w=warndlg(sprintf('\n\nPlease run Brainstorm first and\nthen call SimMEEG within Brainstorm.\n'),'Brainstorm is NOT Running!');
    w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = 11; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
end