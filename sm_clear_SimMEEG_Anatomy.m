function sm_clear_SimMEEG_Anatomy(varargin)
global h

btn = questdlg(sprintf('Are you sure you want to CLEAR all ANATOMY?'),'Clear Anatomy?','Yes','No','No');

switch btn
    case 'Yes'
        
        h.anatomy = struct('mri','','sens_eeg','','sens_meg','','mesh_volumes',[],'headmodel_eeg_vol','','headmodel_meg_vol','','headmodel_eeg_cortex','','headmodel_meg_cortex','','leadfield_eeg_vol','','leadfield_eeg_cortex','','leadfield_meg_vol','','leadfield_meg_cortex','','sens','','headmodel','','leadfield','','mesh_cortex','');
        h.anat_path = ' ';
        h.anat_file = ' ';
        update_anatomy_fields;
    case 'No'
end


