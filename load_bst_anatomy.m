function load_bst_anatomy(varargin)
global h
global bst

try
    bst_check_internet;
catch
    warndlg('Brainstorm must be open to load in Brainstorm Anatomy','Please open BrainStorm'); 
    return
end

bst = [];

%% Create popup figure
% delete(bst.popup_fig)
bst.popup_fig = figure('Position',[range(h.main_fig.Position([1 3]))/10 range(h.main_fig.Position([2 4]))/10 600 600],'Units','normalized','Name','Brainstorm Anatomy'); 
bst.popup_fig.MenuBar='none'; %bst.popup_fig.CloseRequestFcn = @close_bst_anatomy;

load_bst_anatomy_popup();

waitfor(bst.popup_fig);

h.waitfor_txt.String = sprintf('Converting Anatomy to SimMEEG format'); drawnow;
%% %%%%% CONVERT to SimMEEG format %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% convert all units to 'mm'
h.anatomy.mri = ft_convert_units(h.anatomy.mri,'mm');
for m=1:length(h.anatomy.mesh_volumes)
    h.anatomy.mesh_volumes(m) = ft_convert_units(h.anatomy.mesh_volumes(m),'mm');
end

h.anatomy.headmodel_eeg_cortex = ft_convert_units(h.anatomy.headmodel_eeg_cortex,'mm');
h.anatomy.headmodel_eeg_vol = ft_convert_units(h.anatomy.headmodel_eeg_vol,'mm');
h.anatomy.headmodel_meg_cortex = ft_convert_units(h.anatomy.headmodel_meg_cortex,'mm');
h.anatomy.headmodel_meg_vol = ft_convert_units(h.anatomy.headmodel_meg_vol,'mm');

h.anatomy.leadfield_eeg_cortex = ft_convert_units(h.anatomy.leadfield_eeg_cortex,'mm');
h.anatomy.leadfield_eeg_vol = ft_convert_units(h.anatomy.leadfield_eeg_vol,'mm');
h.anatomy.leadfield_meg_cortex = ft_convert_units(h.anatomy.leadfield_meg_cortex,'mm');
h.anatomy.leadfield_meg_vol = ft_convert_units(h.anatomy.leadfield_meg_vol,'mm');




%% make sens compatbile with SimMEEG
if ~isempty(h.anatomy.sens_meg)
    h.anatomy.sens_meg = ft_convert_units(h.anatomy.sens_meg,'mm');
    h.anatomy.sens_meg = sm_make_sens_compatible(h.anatomy.sens_meg,'meg');
    cfg.sens_type = 'meg'; cfg.sens_idx = 1:length(h.anatomy.sens_meg.label); h.anatomy.sens_meg = bs_select_meeg_sensors(cfg,h.anatomy.sens_meg);
end

if ~isempty(h.anatomy.sens_eeg)
    h.anatomy.sens_eeg = ft_convert_units(h.anatomy.sens_eeg,'mm');
    h.anatomy.sens_eeg = sm_make_sens_compatible(h.anatomy.sens_eeg,'eeg');
    cfg.sens_type = 'eeg'; cfg.sens_idx = 1:length(h.anatomy.sens_eeg.label); h.anatomy.sens_eeg = bs_select_meeg_sensors(cfg,h.anatomy.sens_eeg);
end

%% add BRANELab leadfield format
if ~isempty(h.anatomy.leadfield_eeg_cortex) 
x=cell2mat(h.anatomy.leadfield_eeg_cortex.leadfield(h.anatomy.leadfield_eeg_cortex.inside==1));
h.anatomy.leadfield_eeg_cortex.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
%leadfield.mesh_lf=reshape(x,[size(elec.chanpos,1) 3 size(x,2)/3]);
h.anatomy.leadfield_eeg_cortex.voxel_pos=h.anatomy.leadfield_eeg_cortex.pos(h.anatomy.leadfield_eeg_cortex.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
end

if ~isempty(h.anatomy.leadfield_meg_cortex)
x=cell2mat(h.anatomy.leadfield_meg_cortex.leadfield(h.anatomy.leadfield_meg_cortex.inside==1));
h.anatomy.leadfield_meg_cortex.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
%leadfield.mesh_lf=reshape(x,[size(elec.chanpos,1) 3 size(x,2)/3]);
h.anatomy.leadfield_meg_cortex.voxel_pos=h.anatomy.leadfield_meg_cortex.pos(h.anatomy.leadfield_meg_cortex.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
end

if ~isempty(h.anatomy.leadfield_eeg_vol)
x=cell2mat(h.anatomy.leadfield_eeg_vol.leadfield(h.anatomy.leadfield_eeg_vol.inside==1));
h.anatomy.leadfield_eeg_vol.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
%leadfield.mesh_lf=reshape(x,[size(elec.chanpos,1) 3 size(x,2)/3]);
h.anatomy.leadfield_eeg_vol.voxel_pos=h.anatomy.leadfield_eeg_vol.pos(h.anatomy.leadfield_eeg_vol.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
end

if ~isempty(h.anatomy.leadfield_meg_vol)
x=cell2mat(h.anatomy.leadfield_meg_vol.leadfield(h.anatomy.leadfield_meg_vol.inside==1));
h.anatomy.leadfield_meg_vol.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
%leadfield.mesh_lf=reshape(x,[size(elec.chanpos,1) 3 size(x,2)/3]);
h.anatomy.leadfield_meg_vol.voxel_pos=h.anatomy.leadfield_meg_vol.pos(h.anatomy.leadfield_meg_vol.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
end

%% set to existing channels (MEG default)
if ~isempty(h.anatomy.sens_meg)
    h.anatomy.sens = h.anatomy.sens_meg;
    h.anatomy.headmodel = h.anatomy.headmodel_meg_vol;
    h.anatomy.leadfield = h.anatomy.leadfield_meg_vol;
    h.anatomy.mesh_cortex = h.anatomy.mesh_volumes(5);
    h.menu_sens_type.Value = 1;
else
  h.menu_sens_type.Enable = 'inactive';
end

if ~isempty(h.anatomy.sens_eeg)
    h.anatomy.sens = h.anatomy.sens_eeg;
    h.anatomy.headmodel = h.anatomy.headmodel_eeg_vol;
    h.anatomy.leadfield = h.anatomy.leadfield_eeg_vol;
    h.anatomy.mesh_cortex = h.anatomy.mesh_volumes(5);
    h.menu_sens_type.Value = 2;
else
      h.menu_sens_type.Enable = 'inactive';
end
update_anatomy_fields;





