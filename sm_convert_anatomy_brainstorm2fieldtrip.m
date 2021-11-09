%% converting to field trip format

%% ********** Brainstorm must be open to run the following scripts
cd C:\BRANELab\matlab_progs\general_progs\brainstorm3; brainstorm; 
addpath C:\Users\ath\OneDrive\SimMEEG\SimMEEG_v19_under_construction\; 
addpath C:\BRANELab\matlab_progs\general_progs\spm12;
cd C:\Users\ath\OneDrive\SimMEEG\fieldtrip-20200911\; ft_defaults;

%% folders and files
anatdir = 'F:\brainstorm_db\SimMEEG_anat\anat\@default_subject';
datadir = 'F:\brainstorm_db\SimMEEG_anat\data\S001_MEG_EEG\@default_study';
chan_file = fullfile(datadir,'channel_ctf_acc1.mat'); 
MriFile = fullfile(anatdir,'subjectimage_T1.mat'); 
scalpFile = fullfile(anatdir,'tess_head.mat');      % mesh_volumes(1)
skullFile = fullfile(anatdir,'tess_outerskull.mat'); % mesh_volumes(2)
brainFile = fullfile(anatdir,'tess_innerskull.mat'); % mesh_volumes(3)
cortexFile = fullfile(anatdir,'tess_cortex_pial_low.mat'); % mesh_volumes(4)
wmFile = fullfile(anatdir,'tess_cortex_white_low.mat'); % mesh_volumes(5)
pialFile = fullfile(anatdir,'tess_cortex_pial_cereb_low.mat'); % mesh_cortex
%% combined MEG+EEG Headmodel and Leadfields for volume MRI

% renaming MEG and EEG to 'MEEG' to include as single leadfield for combined MEG+EEG
chan_mat=load(chan_file);
meg_idx = find(strcmp({chan_mat.Channel(:).Type},'MEG')==1); 
eeg_idx = find(strcmp({chan_mat.Channel(:).Type},'EEG')==1); 
% creating headmodel and leadfield
include_megref=0; % (0)=do not include MEG ref sensors in leadfield, (1)=include MEG ref sensors in leadfield,

%% %%%% MRI Volume(5mm grid) Leadfields %%%%%%%%%%%%%%%%%%
%% MEG leadfield from MRI volume
[headmodel_meg_vol, leadfield_meg_vol] = out_fieldtrip_headmodel(fullfile(datadir,'headmodel_vol_openmeeg_os_meg.mat'),chan_file,meg_idx,include_megref);
%% EEG leadfield from MRI volume
[headmodel_eeg_vol, leadfield_eeg_vol] = out_fieldtrip_headmodel(fullfile(datadir,'headmodel_vol_openmeeg_os_meg.mat'),chan_mat,eeg_idx,include_megref);
%% combine MEG+EEG leadfield from MRI volume
headmodel_meeg_vol = headmodel_eeg_vol; 
leadfield_meeg_vol = leadfield_meg_vol; 
leadfield_meeg_vol.label = [leadfield_meg_vol.label leadfield_eeg_vol.label];
for v=1:length(leadfield_meeg_vol.leadfield)
    leadfield_meeg_vol.leadfield{v} = [leadfield_meg_vol.leadfield{v}; leadfield_eeg_vol.leadfield{v}];
end

%% %%%% SURFACE(cortex) Leadfields %%%%%%%%%%%%%%%%%%
%% MEG leadfield from MRI surface (cortex)
[headmodel_meg_cortex, leadfield_meg_cortex] = out_fieldtrip_headmodel(fullfile(datadir,'headmodel_surf_openmeeg_os_meg.mat'),chan_mat,meg_idx,include_megref);
%% EEG leadfield from MRI surface (cortex)
[headmodel_eeg_cortex, leadfield_eeg_cortex] = out_fieldtrip_headmodel(fullfile(datadir,'headmodel_surf_openmeeg_os_meg.mat'),chan_mat,eeg_idx,include_megref);


   
%% Converting Channels
[sens_eeg, sens_meg] = out_fieldtrip_channel(chan_mat, include_megref);
sens_eeg.type = 'eeg'; sens_meg.type = 'meg';



%% Loading MRI and meshes
mri = out_fieldtrip_mri(MriFile, 'anatomy');
mesh_volumes(1) = out_fieldtrip_tess(scalpFile);
mesh_volumes(2) = out_fieldtrip_tess(skullFile);
mesh_volumes(3) = out_fieldtrip_tess(brainFile);
mesh_volumes(4) = out_fieldtrip_tess(cortexFile);
mesh_volumes(5) = out_fieldtrip_tess(wmFile);
mesh_cortex = out_fieldtrip_tess(pialFile);


%% converting all to mni coordinate system needed for SimMEEG - shifting [X Y Z] to [Y X Z] 
% mri = ft_convert_coordsys(mri,'acpc');
sens_eeg.chanpos = sens_eeg.chanpos(:,[2 1 3]); sens_eeg.elecpos = sens_eeg.elecpos(:,[2 1 3]);
sens_meg.chanpos = sens_meg.chanpos(:,[2 1 3]); sens_meg.coilpos = sens_meg.coilpos(:,[2 1 3]); %sens_meg.coilori = sens_meg.coilori(:,[2 1 3]);
leadfield_meg_vol.pos = leadfield_meg_vol.pos(:,[2 1 3]);
leadfield_meg_cortex.pos = leadfield_meg_cortex.pos(:,[2 1 3]);
leadfield_eeg_vol.pos = leadfield_eeg_vol.pos(:,[2 1 3]);
leadfield_eeg_cortex.pos = leadfield_eeg_cortex.pos(:,[2 1 3]);
mesh=mesh_volumes; for a=1:5; mesh_volumes(a).pos=mesh(a).pos(:,[2 1 3]); end 

%% combing eeg and meg
sens_meeg = sens_meg;
sens_meeg.label = [sens_meg.label sens_eeg.label];
sens_meeg.chantype = [sens_meg.chantype repmat({'eeg'},1,length(sens_eeg.label))];
sens_meeg.chanpos = cat(1,sens_meg.chanpos,sens_eeg.chanpos);
sens_meeg.elecpos = sens_eeg.elecpos;
%% combine MEG+EEG leadfield from MRI surface (cortex)
headmodel_meeg_cortex = headmodel_eeg_cortex; 
leadfield_meeg_cortex = leadfield_meg_cortex; 
leadfield_meeg_cortex.label = [leadfield_meg_cortex.label leadfield_eeg_cortex.label];
for v=1:length(leadfield_meeg_cortex.leadfield)
    leadfield_meeg_cortex.leadfield{v} = [leadfield_meg_cortex.leadfield{v}; leadfield_eeg_cortex.leadfield{v}];
end



%% saving anatomy
sname = fullfile(datadir,'S001_SimMEEG_Anatomy_MEEG_v1.mat'); 
fprintf('Saving File: %s\n',sname); 
save(sname,'headmodel_*','leadfield_*','sens_*','mri','mesh_*'); 


