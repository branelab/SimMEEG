function sm_bst2ft_anatomy_from_bst_files(bst)
% This program is linked with SimMEEG_v20+ in order to load in BrainStorm Anatomy when SimMEEG is directly called from Brainstorm.
% It will load and convert the files to FieldTrip and SimMEEG formats to work in SimMEEG.
%
% 'bst' is a structure of file names for anatomy and sensors. 
%   If the type of file doesn't exist for the subject then put it as [] --> (e.g., bst.subj_cortexFile = [];)  
%
%   bst.FieldTrip_dir = FieldTrip's directory version "fieldtrip-20200911" tested, may work with newer versions but not tested (this are not used in this program but within SimMEEG main program)
%   bst.subj_anat_dir = Brainstorm's subject anatomy directory file (this are not used in this program but within SimMEEG main program)
%   bst.subj_data_dir = Brainstorm's subject data directory file with sensor files (this are not used in this program but within SimMEEG main program)
%
%   bst.subj_MriFile = Brainstorm's subject mri .mat file (needs to be fullfile with directory included) 
%   bst.subj_scalpFile = Brainstorm's subject scalp mesh .mat file
%   bst.subj_skullFile = Brainstorm's subject skull mesh .mat file
%   bst.subj_brainFile = Brainstorm's subject brain hull mesh .mat file
%   bst.subj_cortexFile = Brainstorm's subject brain cortex mesh .mat file (the one used for generating surface headmodel & leadfields)
%   bst.subj_wmFile = Brainstorm's subject white matter mesh .mat file
%   bst.subj_pialFile = Brainstorm's subject pial matter mesh .mat file
% 
%   bst.subj_sens_meg_file = Brainstorm's subject MEG sensor .mat file
%   bst.subj_hdm_meg_vol_file = Brainstorm's subject MEG HeadModel Volume .mat file
%   bst.subj_hdm_meg_cortex_file = Brainstorm's subject MEG HeadModel Surface (Cortex) .mat file
% 
%   bst.subj_sens_eeg_file = Brainstorm's subject EEG sensor .mat file
%   bst.subj_hdm_eeg_vol_file = Brainstorm's subject EEG HeadModel Volume .mat file
%   bst.subj_hdm_eeg_cortex_file = Brainstorm's subject EEG HeadModel Surface (Cortex) .mat file


% global bst
global h

if ~isappdata(0, 'BrainstormRunning')   % Brainstorm not running
    w=warndlg(sprintf('\n\nPlease run Brainstorm first and\nthen call SimMEEG within Brainstorm.\n'),'Brainstorm is NOT Running!');
    w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = 11; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
    return
else

h.anatomy = struct('mri','','sens_eeg','','sens_meg','','mesh_volumes',[],'headmodel_eeg_vol','','headmodel_meg_vol','','headmodel_eeg_cortex','','headmodel_meg_cortex','','leadfield_eeg_vol','','leadfield_eeg_cortex','','leadfield_meg_vol','','leadfield_meg_cortex','','sens','','headmodel','','leadfield','','mesh_cortex','');

include_megref=0;
%% %%%%% Converting ANATOMY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MRI
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
if ~exist(bst.subj_scalpFile,'file'); warndlg(sprintf('MRI file does not exist\n\n%s\n',bst.subj_scalpFile),'Warning!');
else
    fprintf('Loading "Scalp" file and converting to FieldTrip format: %s\n',bst.subj_scalpFile);
    h.anatomy.mesh_volumes = out_fieldtrip_tess(bst.subj_scalpFile);
end
%% Skull
if ~exist(bst.subj_skullFile,'file'); warndlg(sprintf('Skull file does not exist\n\n%s\n',bst.subj_skullFile),'Warning!');
else
    fprintf('Loading "Skull" file and converting to FieldTrip format: %s\n',bst.subj_skullFile);
    h.anatomy.mesh_volumes(2) = out_fieldtrip_tess(bst.subj_skullFile);
end

%% Brain Hull
if ~exist(bst.subj_brainFile,'file'); warndlg(sprintf('Brain Hull file does not exist\n\n%s\n',bst.subj_brainFile),'Warning!');
else
    fprintf('Loading "Brain Hull" file and converting to FieldTrip format: %s\n',bst.subj_brainFile);
    h.anatomy.mesh_volumes(3) = out_fieldtrip_tess(bst.subj_brainFile);
end
%% Brain Cortex
if ~exist(bst.subj_cortexFile,'file'); warndlg(sprintf('Brain Cortex file does not exist\n\n%s\n',bst.subj_cortexFile),'Warning!');
else
    fprintf('Loading "Brain Cortex" file and converting to FieldTrip format: %s\n',bst.subj_cortexFile);
    h.anatomy.mesh_volumes(4) = out_fieldtrip_tess(bst.subj_cortexFile);
end
%% Brain White Matter
if ~exist(bst.subj_wmFile,'file'); warndlg(sprintf('White Matter file does not exist\n\n%s\n',bst.subj_wmFile),'Warning!');
else
    fprintf('Loading "White Matter" file and converting to FieldTrip format: %s\n',bst.subj_wmFile);
    h.anatomy.mesh_volumes(5) = out_fieldtrip_tess(bst.subj_wmFile);
end
%% Brain Cortex High-res
if ~exist(bst.subj_pialFile,'file'); warndlg(sprintf('Brain Pial file does not exist\n\n%s\n',bst.subj_pialFile),'Warning!');
else
    fprintf('Loading "Pial Surface" file and converting to FieldTrip format: %s\n',bst.subj_pialFile);
    h.anatomy.mesh_cortex = out_fieldtrip_tess(bst.subj_pialFile);
end

%% %%%%% Converting Sensors, HeadModels & LeadFields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sensor MEG
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
        h.menu_sens_type.Value=2; h.anatomy.sens_meg =[];
        warning( sprintf('MEG Sensors do no exist in file: %s',bst.subj_sens_meg_file) )
    end
end
%% MEG HeadModel Volume
if exist(bst.subj_hdm_meg_vol_file,'file') && ~isempty(meg_idx)
    fprintf('Loading "MEG HeadModel & LeadField: Volume Surface" file and converting to FieldTrip format: %s\n',bst.subj_hdm_meg_vol_file);
    [h.anatomy.headmodel_meg_vol, h.anatomy.leadfield_meg_vol] = out_fieldtrip_headmodel(bst.subj_hdm_meg_vol_file,chan_mat,meg_idx,include_megref);
else
    warning(sprintf('No MEG sensors found or  HeadModel MEG Volume file does not exist:   %s',bst.subj_hdm_meg_vol_file));
    h.anatomy.headmodel_meg_vol=[]; h.anatomy.leadfield_meg_vol=[];
end
%% MEG HeadModel Cortex
if exist(bst.subj_hdm_meg_cortex_file,'file') && ~isempty(meg_idx)
    fprintf('Loading "MEG HeadModel & LeadField: Cortical Surface" file and converting to FieldTrip format: %s\n',bst.subj_hdm_meg_cortex_file);
    [h.anatomy.headmodel_meg_cortex, h.anatomy.leadfield_meg_cortex] = out_fieldtrip_headmodel(bst.subj_hdm_meg_cortex_file,chan_mat,meg_idx,include_megref);
else
    warning(sprintf('No MEG sensors found or HeadModel MEG Cortex file does not exist:   %s',bst.subj_hdm_meg_cortex_file));
    h.anatomy.headmodel_meg_cortex=[]; h.anatomy.leadfield_meg_cortex=[];
end

%% Sensor EEG
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
        h.menu_sens_type.Value=1; h.anatomy.sens_eeg =[];
        warning( sprintf('MEG Sensors do no exist in file: %s\n',bst.subj_sens_eeg_file) )
    end
end
%% EEG HeadModel Volume
if exist(bst.subj_hdm_eeg_vol_file,'file') && ~isempty(eeg_idx)
    fprintf('Loading "EEG HeadModel & LeadField: Volume Surface" file and converting to FieldTrip format: %s\n',bst.subj_hdm_eeg_vol_file);
    [h.anatomy.headmodel_eeg_vol, h.anatomy.leadfield_eeg_vol] = out_fieldtrip_headmodel(bst.subj_hdm_eeg_vol_file,chan_mat,eeg_idx,include_megref);
else
    warning(sprintf('No EEG sensors found or  HeadModel EEG Volume file does not exist:   %s',bst.subj_hdm_eeg_vol_file));
    h.anatomy.headmodel_eeg_vol=[]; h.anatomy.leadfield_eeg_vol=[];
end
%% EEG HeadModel Cortex
if exist(bst.subj_hdm_eeg_cortex_file,'file') && ~isempty(eeg_idx)
    fprintf('Loading "EEG HeadModel & LeadField: Cortical Surface" file and converting to FieldTrip format: %s\n',bst.subj_hdm_eeg_cortex_file);
    [h.anatomy.headmodel_eeg_cortex, h.anatomy.leadfield_eeg_cortex] = out_fieldtrip_headmodel(bst.subj_hdm_eeg_cortex_file,chan_mat,eeg_idx,include_megref);
else
    warning(sprintf('No EEG sensors found or HeadModel EEG Cortex file does not exist:   %s',bst.subj_hdm_eeg_cortex_file));
    h.anatomy.headmodel_eeg_cortex=[]; h.anatomy.leadfield_eeg_cortex=[];
end


%% converted flag
% w = msgbox(sprintf('\nSuccessfully converted Anatomy Files\n\nPlease click "Close" to load "Anatomy" structure into SimMEEG'));
% w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
bst.converted_flag=1;


%% making anatomy compatible with SimMEEG

if strcmp(h.anat_file,'ANATOMY_DEFAULT_adult_SimMEEG_v16.mat')  % reordering meshc_volumes because leadfiels in v16 were based on white matter (mesh_volumes(5)) and not pial surface (mesh_volumes(4))
    h.anatomy.mesh_cortex = h.anatomy.mesh_volumes(4);
    h.anatomy.mesh_volumes(4) = h.anatomy.mesh_volumes(5);
end

for a=1:length(h.anatomy.mesh_volumes); h.anatomy.mesh_volumes(a) = ft_convert_units(h.anatomy.mesh_volumes(a),'mm'); end
%% convert all units to 'mm'
if ~isempty(h.anatomy.sens_eeg);  h.anatomy.sens_eeg = ft_convert_units(h.anatomy.sens_eeg,'mm'); end
if ~isempty(h.anatomy.sens_meg);  h.anatomy.sens_meg = ft_convert_units(h.anatomy.sens_meg,'mm'); end
if ~isempty(h.anatomy.mri);       h.anatomy.mri = ft_convert_units(h.anatomy.mri,'mm'); end


if ~isempty(h.anatomy.headmodel_eeg_cortex);    h.anatomy.headmodel_eeg_cortex = ft_convert_units(h.anatomy.headmodel_eeg_cortex,'mm'); end
if ~isempty(h.anatomy.headmodel_eeg_vol);       h.anatomy.headmodel_eeg_vol = ft_convert_units(h.anatomy.headmodel_eeg_vol,'mm'); end
if ~isempty(h.anatomy.headmodel_meg_cortex);    h.anatomy.headmodel_meg_cortex = ft_convert_units(h.anatomy.headmodel_meg_cortex,'mm'); end
if ~isempty(h.anatomy.headmodel_meg_vol);       h.anatomy.headmodel_meg_vol = ft_convert_units(h.anatomy.headmodel_meg_vol,'mm'); end


if ~isempty(h.anatomy.leadfield_eeg_cortex); h.anatomy.leadfield_eeg_cortex = ft_convert_units(h.anatomy.leadfield_eeg_cortex,'mm'); end
if ~isempty(h.anatomy.leadfield_eeg_vol); h.anatomy.leadfield_eeg_vol = ft_convert_units(h.anatomy.leadfield_eeg_vol,'mm'); end
if ~isempty(h.anatomy.leadfield_meg_cortex); h.anatomy.leadfield_meg_cortex = ft_convert_units(h.anatomy.leadfield_meg_cortex,'mm'); end
if ~isempty(h.anatomy.leadfield_meg_vol); h.anatomy.leadfield_meg_vol = ft_convert_units(h.anatomy.leadfield_meg_vol,'mm'); end

%% make sens compatbile with SimMEEG
if ~isempty(h.anatomy.sens_meg)
    h.anatomy.sens_meg = sm_make_sens_compatible(h.anatomy.sens_meg,'meg');
    cfg.sens_type = 'meg'; cfg.sens_idx = 1:length(h.anatomy.sens_meg.label); h.anatomy.sens_meg = bs_select_meeg_sensors(cfg,h.anatomy.sens_meg);
end
if ~isempty(h.anatomy.sens_eeg)
    h.anatomy.sens_eeg = sm_make_sens_compatible(h.anatomy.sens_eeg,'eeg');
    cfg.sens_type = 'eeg'; cfg.sens_idx = 1:length(h.anatomy.sens_eeg.label); h.anatomy.sens_eeg = bs_select_meeg_sensors(cfg,h.anatomy.sens_eeg);
end

%% add BRANELab leadfield format
if ~isempty(h.anatomy.leadfield_eeg_cortex)
    x=cell2mat(h.anatomy.leadfield_eeg_cortex.leadfield(h.anatomy.leadfield_eeg_cortex.inside==1));
    h.anatomy.leadfield_eeg_cortex.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
    h.anatomy.leadfield_eeg_cortex.voxel_pos=h.anatomy.leadfield_eeg_cortex.pos(h.anatomy.leadfield_eeg_cortex.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
end
if ~isempty(h.anatomy.leadfield_meg_cortex)
    x=cell2mat(h.anatomy.leadfield_meg_cortex.leadfield(h.anatomy.leadfield_meg_cortex.inside==1));
    h.anatomy.leadfield_meg_cortex.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
    h.anatomy.leadfield_meg_cortex.voxel_pos=h.anatomy.leadfield_meg_cortex.pos(h.anatomy.leadfield_meg_cortex.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
end
if ~isempty(h.anatomy.leadfield_eeg_vol)
    x=cell2mat(h.anatomy.leadfield_eeg_vol.leadfield(h.anatomy.leadfield_eeg_vol.inside==1));
    h.anatomy.leadfield_eeg_vol.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
    h.anatomy.leadfield_eeg_vol.voxel_pos=h.anatomy.leadfield_eeg_vol.pos(h.anatomy.leadfield_eeg_vol.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
end
if ~isempty(h.anatomy.leadfield_meg_vol)
    x=cell2mat(h.anatomy.leadfield_meg_vol.leadfield(h.anatomy.leadfield_meg_vol.inside==1));
    h.anatomy.leadfield_meg_vol.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
    h.anatomy.leadfield_meg_vol.voxel_pos=h.anatomy.leadfield_meg_vol.pos(h.anatomy.leadfield_meg_vol.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
end

%% default starting with MEG
if ~isempty(h.anatomy.sens_meg)
    h.anatomy.sens = h.anatomy.sens_meg; h.menu_sens_type.Value = 1;
    if ~isempty(h.anatomy.headmodel_meg_vol)
        h.anatomy.headmodel = h.anatomy.headmodel_meg_vol;
        h.anatomy.leadfield = h.anatomy.leadfield_meg_vol;
    elseif ~isempty(h.anatomy.headmodel_meg_cortex)
        h.anatomy.headmodel = h.anatomy.headmodel_meg_cortex;
        h.anatomy.leadfield = h.anatomy.leadfield_meg_cortex;
    else
        h.anatomy.headmodel = [];
        h.anatomy.leadfield =[];
    end
elseif ~isempty(h.anatomy.sens_eeg) && isempty(h.anatomy.sens_meg)
    h.anatomy.sens = h.anatomy.sens_eeg; h.menu_sens_type.Value = 1;
    if ~isempty(h.anatomy.headmodel_eeg_vol)
        h.anatomy.headmodel = h.anatomy.headmodel_eeg_vol;
        h.anatomy.leadfield = h.anatomy.leadfield_eeg_vol;
    elseif ~isempty(h.anatomy.headmodel_eeg_cortex)
        h.anatomy.headmodel = h.anatomy.headmodel_eeg_cortex;
        h.anatomy.leadfield = h.anatomy.leadfield_eeg_cortex;
    else
        h.anatomy.headmodel = [];
        h.anatomy.leadfield =[];
    end
end


if length(h.anatomy.mesh_volumes)>3
    h.anatomy.mesh_cortex = h.anatomy.mesh_volumes(5);
else
    h.anatomy.mesh_cortex = h.anatomy.mesh_volumes(3);
end
update_anatomy_fields;


end
