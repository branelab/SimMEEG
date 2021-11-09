function sm_convert_bst2ft(varargin)
global bst
global h

h.anatomy = struct('mri','','sens_eeg','','sens_meg','','mesh_volumes',[],'headmodel_eeg_vol','','headmodel_meg_vol','','headmodel_eeg_cortex','','headmodel_meg_cortex','','leadfield_eeg_vol','','leadfield_eeg_cortex','','leadfield_meg_vol','','leadfield_meg_cortex','','sens','','headmodel','','leadfield','','mesh_cortex','');

include_megref=0;
%% %%%%% Converting ANATOMY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MRI
bst.subj_MriFile = fullfile(bst.subj_anat_dir,bst.anat_mri_file.String{bst.anat_mri_file.Value});
[fpath,fname,fext] = fileparts(bst.subj_MriFile);
h.anat_file = sprintf('%s%s',fname,fext);

if ~exist(bst.subj_MriFile,'file')
    warndlg(sprintf('MRI file does not exist\n\n%s\n',bst.subj_MriFile),'Warning!');
else
    try
        fprintf('Loading "MRI" file and converting to FieldTrip format: %s\n',bst.subj_MriFile);
        h.anatomy.mri = out_fieldtrip_mri(bst.subj_MriFile, 'anatomy');
    catch
        warndlg('You must have BrainStorm open to load in brainstorm anatomy','Warning');
        return
    end
end
%% Scalp
bst.subj_scalpFile = fullfile(bst.subj_anat_dir,bst.anat_scalp_file.String{bst.anat_scalp_file.Value});
if ~exist(bst.subj_scalpFile,'file'); warndlg(sprintf('MRI file does not exist\n\n%s\n',bst.subj_scalpFile),'Warning!');
else
    fprintf('Loading "Scalp" file and converting to FieldTrip format: %s\n',bst.subj_scalpFile);
    h.anatomy.mesh_volumes = out_fieldtrip_tess(bst.subj_scalpFile);
end
%% Skull
bst.subj_skullFile = fullfile(bst.subj_anat_dir,bst.anat_skull_file.String{bst.anat_skull_file.Value});
if ~exist(bst.subj_skullFile,'file'); warndlg(sprintf('Skull file does not exist\n\n%s\n',bst.subj_skullFile),'Warning!');
else
    fprintf('Loading "Skull" file and converting to FieldTrip format: %s\n',bst.subj_skullFile);
    h.anatomy.mesh_volumes(2) = out_fieldtrip_tess(bst.subj_skullFile);
end
%% Brain Hull
bst.subj_brainFile = fullfile(bst.subj_anat_dir,bst.anat_brain_hull_file.String{bst.anat_brain_hull_file.Value});
if ~exist(bst.subj_brainFile,'file'); warndlg(sprintf('Brain Hull file does not exist\n\n%s\n',bst.subj_brainFile),'Warning!');
else
    fprintf('Loading "Brain Hull" file and converting to FieldTrip format: %s\n',bst.subj_brainFile);
    h.anatomy.mesh_volumes(3) = out_fieldtrip_tess(bst.subj_brainFile);
end
%% Brain Cortex
bst.subj_cortexFile = fullfile(bst.subj_anat_dir,bst.anat_brain_cortex_file.String{bst.anat_brain_cortex_file.Value});
if ~exist(bst.subj_cortexFile,'file'); warndlg(sprintf('Brain Cortex file does not exist\n\n%s\n',bst.subj_cortexFile),'Warning!');
else
    fprintf('Loading "Brain Cortex" file and converting to FieldTrip format: %s\n',bst.subj_cortexFile);
    h.anatomy.mesh_volumes(4) = out_fieldtrip_tess(bst.subj_cortexFile);
end
%% Brain White Matter
bst.subj_wmFile = fullfile(bst.subj_anat_dir,bst.anat_brain_wm_file.String{bst.anat_brain_wm_file.Value});
if ~exist(bst.subj_wmFile,'file'); warndlg(sprintf('White Matter file does not exist\n\n%s\n',bst.subj_wmFile),'Warning!');
else
    fprintf('Loading "White Matter" file and converting to FieldTrip format: %s\n',bst.subj_wmFile);
    h.anatomy.mesh_volumes(5) = out_fieldtrip_tess(bst.subj_wmFile);
end
%% Brain Cortex High-res
bst.subj_pialFile = fullfile(bst.subj_anat_dir,bst.anat_brain_pial_file.String{bst.anat_brain_pial_file.Value});
if ~exist(bst.subj_pialFile,'file'); warndlg(sprintf('Brain Pial file does not exist\n\n%s\n',bst.subj_pialFile),'Warning!');
else
    fprintf('Loading "Pial Surface" file and converting to FieldTrip format: %s\n',bst.subj_pialFile);
    h.anatomy.mesh_cortex = out_fieldtrip_tess(bst.subj_pialFile);
end

%% %%%%% Converting Sensors, HeadModels & LeadFields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sensor MEG
bst.subj_sens_meg_file = fullfile(bst.subj_data_dir,bst.anat_sens_meg_file.String{bst.anat_sens_meg_file.Value});
meg_idx = [];
if ~exist(bst.subj_sens_meg_file,'file')
    warning(sprintf('MEG Sensor file does not exist:   %s',bst.subj_sens_meg_file));
else
    fprintf('Loading "MEG Channel" file and converting to FieldTrip format: %s\n',bst.subj_sens_meg_file);
    chan_mat=load(bst.subj_sens_meg_file);
    meg_idx = find(strcmpi({chan_mat.Channel(:).Type},'MEG')==1);
    if ~isempty(meg_idx)
        include_megref=0; % (0)=do not include MEG ref sensors in leadfield, (1)=include MEG ref sensors in leadfield,
        [~, h.anatomy.sens_meg] = out_fieldtrip_channel(chan_mat, include_megref);
        if ~isempty(h.anatomy.sens_meg); h.anatomy.sens_meg.type = 'meg'; end
    else
        bst.anat_sens_meg_file.Value=1; h.anatomy.sens_meg =[];
        warning( sprintf('MEG Sensors do no exist in file: %s',bst.subj_sens_meg_file) )
    end
end
%% MEG HeadModel Volume
bst.subj_hdm_meg_vol_file = fullfile(bst.subj_data_dir,bst.anat_hdm_meg_vol_file.String{bst.anat_hdm_meg_vol_file.Value});
if exist(bst.subj_hdm_meg_vol_file,'file') && ~isempty(meg_idx)
    fprintf('Loading "MEG HeadModel & LeadField: Volume Surface" file and converting to FieldTrip format: %s\n',bst.subj_hdm_meg_vol_file);
    [h.anatomy.headmodel_meg_vol, h.anatomy.leadfield_meg_vol] = out_fieldtrip_headmodel(bst.subj_hdm_meg_vol_file,chan_mat,meg_idx,include_megref);
else
    warning(sprintf('No MEG sensors found or  HeadModel MEG Volume file does not exist:   %s',bst.subj_hdm_meg_vol_file));
    h.anatomy.headmodel_meg_vol=[]; h.anatomy.leadfield_meg_vol=[];
end
%% MEG HeadModel Cortex
bst.subj_hdm_meg_cortex_file = fullfile(bst.subj_data_dir,bst.anat_hdm_meg_cortex_file.String{bst.anat_hdm_meg_cortex_file.Value});
if exist(bst.subj_hdm_meg_cortex_file,'file') && ~isempty(meg_idx)
    fprintf('Loading "MEG HeadModel & LeadField: Cortical Surface" file and converting to FieldTrip format: %s\n',bst.subj_hdm_meg_cortex_file);
    [h.anatomy.headmodel_meg_cortex, h.anatomy.leadfield_meg_cortex] = out_fieldtrip_headmodel(bst.subj_hdm_meg_cortex_file,chan_mat,meg_idx,include_megref);
else
    warning(sprintf('No MEG sensors found or HeadModel MEG Cortex file does not exist:   %s',bst.subj_hdm_meg_cortex_file));
    h.anatomy.headmodel_meg_cortex=[]; h.anatomy.leadfield_meg_cortex=[];
end

%% Sensor EEG
bst.subj_sens_eeg_file = fullfile(bst.subj_data_dir,bst.anat_sens_eeg_file.String{bst.anat_sens_eeg_file.Value});
eeg_idx = [];
if ~exist(bst.subj_sens_eeg_file,'file')
    warning(sprintf('EEG Sensor file does not exist:   %s',bst.subj_sens_eeg_file));
else
    fprintf('Loading "EEG Channel" file and converting to FieldTrip format: %s\n',bst.subj_sens_eeg_file);
    chan_mat=load(bst.subj_sens_eeg_file);
    eeg_idx = find(strcmpi({chan_mat.Channel(:).Type},'EEG')==1);
    if ~isempty(eeg_idx)
        include_megref=0; % (0)=do not include MEG ref sensors in leadfield, (1)=include MEG ref sensors in leadfield,
        [h.anatomy.sens_eeg,~] = out_fieldtrip_channel(chan_mat, include_megref);
        if ~isempty(h.anatomy.sens_eeg); h.anatomy.sens_eeg.type = 'eeg'; end
    else
        bst.anat_sens_eeg_file.Value=1; h.anatomy.sens_eeg =[];
        warning( sprintf('MEG Sensors do no exist in file: %s\n',bst.subj_sens_eeg_file) )
    end
end
%% EEG HeadModel Volume
bst.subj_hdm_eeg_vol_file = fullfile(bst.subj_data_dir,bst.anat_hdm_eeg_vol_file.String{bst.anat_hdm_eeg_vol_file.Value});
if exist(bst.subj_hdm_eeg_vol_file,'file') && ~isempty(eeg_idx)
    fprintf('Loading "EEG HeadModel & LeadField: Volume Surface" file and converting to FieldTrip format: %s\n',bst.subj_hdm_eeg_vol_file);
    [h.anatomy.headmodel_eeg_vol, h.anatomy.leadfield_eeg_vol] = out_fieldtrip_headmodel(bst.subj_hdm_eeg_vol_file,chan_mat,eeg_idx,include_megref);
else
    warning(sprintf('No EEG sensors found or  HeadModel EEG Volume file does not exist:   %s',bst.subj_hdm_eeg_vol_file));
    h.anatomy.headmodel_eeg_vol=[]; h.anatomy.leadfield_eeg_vol=[];
end
%% EEG HeadModel Cortex
bst.subj_hdm_eeg_cortex_file = fullfile(bst.subj_data_dir,bst.anat_hdm_eeg_cortex_file.String{bst.anat_hdm_eeg_cortex_file.Value});
if exist(bst.subj_hdm_eeg_cortex_file,'file') && ~isempty(eeg_idx)
    fprintf('Loading "EEG HeadModel & LeadField: Cortical Surface" file and converting to FieldTrip format: %s\n',bst.subj_hdm_eeg_cortex_file);
    [h.anatomy.headmodel_eeg_cortex, h.anatomy.leadfield_eeg_cortex] = out_fieldtrip_headmodel(bst.subj_hdm_eeg_cortex_file,chan_mat,eeg_idx,include_megref);
else
    warning(sprintf('No EEG sensors found or HeadModel EEG Cortex file does not exist:   %s',bst.subj_hdm_eeg_cortex_file));
    h.anatomy.headmodel_eeg_cortex=[]; h.anatomy.leadfield_eeg_cortex=[];
end


%% converted flag
w = msgbox(sprintf('\nSuccessfully converted Anatomy Files\n\nPlease click "Close" to load "Anatomy" structure into SimMEEG'));
w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
bst.converted_flag=1;
