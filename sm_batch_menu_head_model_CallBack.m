function h = sm_batch_menu_head_model_CallBack(h)

h.listbox_chans.Value = 1; % resetting when switching between MEG and EEG sensors beccause of different # sensors

if h.menu_ori_normal.Value == 2 && h.menu_head_model.Value == 1
    h.menu_ori_normal.Value = 1; % Volume Headmodel can only have random orientations
    for v=1:3; h.edit_source_ori(v).ForegroundColor = h.src_clr(v,:); h.edit_source_ori(v).Enable = 'on'; end
end

h = sm_batch_menu_sens_type_CallBack(h);


if h.menu_head_model.Value == 1     % Volume
    if h.menu_sens_type.Value == 1 && ~isempty(h.anatomy.leadfield_meg_vol) && ~isempty(h.anatomy.headmodel_meg_vol)  % MEG
        h.anatomy.leadfield = h.anatomy.leadfield_meg_vol(h.menu_sens_montage.Value);
        h.anatomy.headmodel = h.anatomy.headmodel_meg_vol(h.menu_sens_montage.Value);
    elseif h.menu_sens_type.Value == 2 && ~isempty(h.anatomy.leadfield_eeg_vol) && ~isempty(h.anatomy.headmodel_eeg_vol) % EEG
        h.anatomy.leadfield = h.anatomy.leadfield_eeg_vol(h.menu_sens_montage.Value);
        h.anatomy.headmodel = h.anatomy.headmodel_eeg_vol(h.menu_sens_montage.Value);
    elseif h.menu_sens_type.Value == 3  % M/EEG
        h.anatomy.leadfield = [];  % will need to change this in the future for adding M/EEG
        h.anatomy.headmodel = [];
    end
    %     h.menu_ori_normal.Enable = 'inactive';  h.menu_ori_normal.ForegroundColor=[1 1 1]*0.5;   h.menu_ori_normal.Value = 1;    % Random orientations locked
elseif h.menu_head_model.Value == 2 && ...
        ( ~isempty(h.anatomy.leadfield_meg_cortex) && ~isempty(h.anatomy.headmodel_meg_cortex) || ... % Cortical Surface
        ~isempty(h.anatomy.leadfield_eeg_cortex) && ~isempty(h.anatomy.headmodel_eeg_cortex) )
    if h.menu_sens_type.Value == 1  % MEG
        h.anatomy.leadfield = h.anatomy.leadfield_meg_cortex(h.menu_sens_montage.Value);
        h.anatomy.headmodel = h.anatomy.headmodel_meg_cortex(h.menu_sens_montage.Value);
    elseif h.menu_sens_type.Value == 2  % EEG
        h.anatomy.leadfield = h.anatomy.leadfield_eeg_cortex(h.menu_sens_montage.Value);
        h.anatomy.headmodel = h.anatomy.headmodel_eeg_cortex(h.menu_sens_montage.Value);
        
    elseif h.menu_sens_type.Value == 3  % M/EEG
        h.anatomy.leadfield = [];  % will need to change this in the future for adding M/EEG
        h.anatomy.headmodel = [];
    end
    if ~isfield(h.anatomy.leadfield,'ori')
        h.anatomy.leadfield.ori = normals(h.anatomy.leadfield.voxel_pos, h.anatomy.mesh_volumes(4).tri);
    end
    %     h.menu_ori_normal.Enable = 'inactive';  h.menu_ori_normal.ForegroundColor=[1 1 1]*.5; h.menu_ori_normal.Value = 2;    % Random orientations locked
end
h = sm_batch_menu_sens_type_CallBack(h);
