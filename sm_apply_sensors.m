function [x]=sm_apply_sensors(varargin)
global h
lbl = h.sens_table.Data(:,1);
x = cell2mat(h.sens_table.Data(:,2:4));
sens_type = h.sens_table.Data(:,5);
if h.menu_sens_type.Value == 1      %MEG
    msgbox('Changing MEG sensors is not available yet.');
    %     h.anatomy.sens_meg.elecpos = x;
    
elseif h.menu_sens_type.Value == 2      %EEG
    h.anatomy.sens_eeg.elecpos = x;
    h.anatomy.sens_eeg.chanpos = x;
    h.anatomy.sens_eeg.label = lbl;
    h.anatomy.sens_eeg.tra = eye(size(x,1)); % h.anatomy.sens_eeg.tra(1:size(x,1),1:size(x,1));
    h.anatomy.sens_eeg.good_sensors = 1:size(x,1);
    h.anatomy.sens_eeg.type = 'eeg';
    h.anatomy.sens_eeg.cfg=[];
    h.anatomy.sens_eeg.chantype = sens_type;
    h.anatomy.sens_eeg.chanunit = repmat({'V'},size(x,1),1);

    h.anatomy.sens = h.anatomy.sens_eeg;
    h.anatomy.leadfield = []; h.anatomy.leadfield_eeg_vol =[]; h.anatomy.leadfield_eeg_cortex = [];
    h.anatomy.headmodel = []; h.anatomy.headmodel_eeg_cortex = []; h.anatomy.headmodel_eeg_vol = []; h.anatomy.headmodel_meg_cortex = []; h.anatomy.headmodel_meg_vol = [];
    h.fcn_handle.menu_head_model_CallBack('');
    h.fcn_handle.plot_3D_mri('');
end


