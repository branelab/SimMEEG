function sm_load_FieldTrip_Anatomy(varargin)
global h;

h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Loading Field Trip''s Anatomy\n\n%s',h.anat_file); drawnow;

[fname,fpath]=uigetfile({'*.mat','FieldTrip Anatomy (*.mat)'; '*.*','All files'});
h.anat_path = fpath;
h.anat_file = fname;
[~,~,fext]=fileparts(h.anat_file);

switch fext
    case '.mat' % if matlab file then load in mri mat file that should also contain the segmentedmri and mesh_volumes
       h.anatomy = struct('mri','','sens_eeg','','sens_meg','','mesh_volumes',[],'headmodel_eeg_vol','','headmodel_meg_vol','','headmodel_eeg_cortex','','headmodel_meg_cortex','','leadfield_eeg_vol','','leadfield_eeg_cortex','','leadfield_meg_vol','','leadfield_meg_cortex','','sens','','headmodel','','leadfield','','mesh_cortex','');

        x = load(fullfile(fpath, fname));
        h.anatomy.mri = x.mri;
        %          = segmentedmri;
        h.anatomy.mesh_volumes = x.mesh_volumes;
        h.anatomy.mesh_volumes_org = x.mesh_volumes;
        msg_string = {'mri' 'mesh_volumes'};
        if isfield(x,'sens_meg')
            h.anatomy.sens_meg = x.sens_meg;
            h.anatomy.sens = x.sens_meg; h.menu_sens_type.Value=1;
            msg_string = [msg_string {'MEG sensors'}];
        end
        if isfield(x,'sens_eeg')
            h.anatomy.sens_eeg = x.sens_eeg;
            h.anatomy.sens = x.sens_eeg; h.menu_sens_type.Value=2;
            msg_string = [msg_string 'EEG sensors'];
        end
        
        %% Volume
        if isfield(x,'headmodel_meg_vol')
            h.anatomy.headmodel_meg_vol = x.headmodel_meg_vol;
            h.anatomy.headmodel = x.headmodel_meg_vol; h.menu_sens_type.Value = 1;
            msg_string = [msg_string 'Head Model MEG Volume'];
        end
        if isfield(x,'headmodel_eeg_vol')
            h.anatomy.headmodel_eeg_vol = x.headmodel_eeg_vol;
            h.anatomy.headmodel = x.headmodel_eeg_vol;  h.menu_sens_type.Value = 2;
            msg_string = [msg_string 'Head Model EEG Volume'];
        end
        
        if isfield(x,'leadfield_meg_vol')
            h.anatomy.leadfield_meg_vol = x.leadfield_meg_vol;
            h.anatomy.leadfield_vol = x.leadfield_meg_vol; h.menu_sens_type.Value = 1;
            msg_string = [msg_string 'Lead Field MEG Volume'];
        end
        if isfield(x,'leadfield_eeg_vol')
            h.anatomy.leadfield_eeg_vol = x.leadfield_eeg_vol;
            h.anatomy.leadfield = x.leadfield_eeg_vol;  h.menu_sens_type.Value = 2;
            msg_string = [msg_string 'Lead Field EEG Volume'];
        end
        
        %% Cortex
        if isfield(x,'headmodel_meg_cortex')
            h.anatomy.headmodel_meg_cortex = x.headmodel_meg_cortex;
            h.anatomy.headmodel = x.headmodel_meg_cortex; h.menu_sens_type.Value = 1;
            msg_string = [msg_string 'Head Model MEG Cortex'];
        end
        if isfield(x,'headmodel_eeg_cortex')
            h.anatomy.headmodel_eeg_cortex = x.headmodel_eeg_cortex;
            h.anatomy.headmodel = x.headmodel_eeg_cortex;  h.menu_sens_type.Value = 2;
            msg_string = [msg_string 'Head Model EEG Cortex'];
        end
        
        if isfield(x,'leadfield_meg_cortex')
            h.anatomy.leadfield_meg_cortex = x.leadfield_meg_cortex;
            h.anatomy.leadfield_cortex = x.leadfield_meg_cortex; h.menu_sens_type.Value = 1;
            msg_string = [msg_string 'Lead Field MEG Cortex'];
        end
        if isfield(x,'leadfield_eeg_cortex')
            h.anatomy.leadfield_eeg_cortex = x.leadfield_eeg_cortex;
            h.anatomy.leadfield = x.leadfield_eeg_cortex;  h.menu_sens_type.Value = 2;
            msg_string = [msg_string 'Lead Field EEG Cortex'];
        end
        msgbox([sprintf('Loaded pre-existing:\n') sprintf('     %s\n',msg_string{:})]);
end

update_anatomy_fields();

h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');

