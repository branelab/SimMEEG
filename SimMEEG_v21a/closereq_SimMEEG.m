function closereq_SimMEEG(varargin)
% Closes SimMEEG and checks if BrainStorm was called
global h
try
if isvalid(h.main_fig)
    if h.bst_called_flag == 1   % Ask User to define what to export to Brainstorm before closing
      
        %% Francois to place code in here on what to export to Brainstorm
        % check to see what exists within h.sim_data

      delete(h.main_fig)
      
    else
        answ = questdlg('Are you sure you want to exit SimMEEG?','Close SimMEEG Study?','Yes','No','No');
        switch answ
            case 'Yes'
                delete(h.main_fig)
            case 'No'
        end
    end
    
else
    delete(h.main_fig);
end
catch
    delete(h.main_fig); %close all; 
end

