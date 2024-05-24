function sm_ARM_set_source_ori(varargin)
global h

switch varargin{end}
    case 'Random'
        h.cfg.ARM_params.vx_ori =[];
        for v=1:length(h.cfg.ARM_params.vx_idx)
            az_el = deg2rad(randi([0 360],1,2)); [x,y,z] = sph2cart(az_el(1),az_el(2),1);
            h.cfg.ARM_params.vx_ori(v,:) = [x y z];  % source orientations (X, Y, Z)
        end
    case 'Cortical Surface'
        if ~isempty(h.anatomy.leadfield_eeg_cortex)
            vx_idx = find_nearest_voxel(h.anatomy.leadfield.voxel_pos(h.cfg.ARM_params.vx_idx,:),h.anatomy.leadfield_eeg_cortex(h.menu_sens_montage.Value).voxel_pos);
            norm_ori = normals(h.anatomy.leadfield_eeg_cortex(h.menu_sens_montage.Value).voxel_pos, h.anatomy.mesh_volumes(4).tri);
            h.cfg.ARM_params.vx_ori = norm_ori(vx_idx,:);
        elseif ~isempty(h.anatomy.leadfield_meg_cortex)
            vx_idx = find_nearest_voxel(h.anatomy.leadfield.voxel_pos(h.cfg.ARM_params.vx_idx,:),h.anatomy.leadfield_meg_cortex(h.menu_sens_montage.Value).voxel_pos);
            norm_ori = normals(h.anatomy.leadfield_meg_cortex(h.menu_sens_montage.Value).voxel_pos, h.anatomy.mesh_volumes(4).tri);
            h.cfg.ARM_params.vx_ori = norm_ori(vx_idx,:);
        else
            warndlg(sprintf('No Cortical Surfaces Exist in loaded EEG/MEG Headmodels.\n\nRandomizing orientations by default\n'));
            sm_ARM_set_source_ori('Random');
        end
        
end