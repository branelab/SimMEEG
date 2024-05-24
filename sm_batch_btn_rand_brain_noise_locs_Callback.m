function h = sm_batch_btn_rand_brain_noise_locs_Callback(h)

% Fully randomize brain noise locations, orientations, and amplitudes for each trial
num_noise = str2num(h.edit_num_noise_sources.String);
brn_vx = setdiff(1:size(h.anatomy.leadfield.voxel_pos,1),h.monte_params.cfg.source.vx_idx); % removing true source voxels from being included in possible brain noise voxels
% random noise amps
noise_amp = str2num(h.edit_brain_noise_amp.String);
h.monte_params.cfg.source.brain_noise_idx = zeros(num_noise,h.monte_params.cfg.study.num_trials);
h.monte_params.cfg.source.brain_noise_amp = h.monte_params.cfg.source.brain_noise_idx;
h.monte_params.cfg.source.brain_noise_ori=ones(num_noise,3,h.monte_params.cfg.study.num_trials);

% rn = randperm(length(brn_vx)); % same locs across trials
%     az = deg2rad(normrnd(0,180,num_noise,1));
%     el = deg2rad(normrnd(0,180,num_noise,1));
for t=1:h.monte_params.cfg.study.num_trials % new set of brain noise sources each trial
    % random locs across trials
    rn = randperm(length(brn_vx));
    h.monte_params.cfg.source.brain_noise_idx(:,t) = brn_vx(rn(1:num_noise)); % brain noise voxel locations
    h.monte_params.cfg.source.brain_noise_amp(:,t) = abs(normrnd(noise_amp(1),noise_amp(2),num_noise,1));  % randomly sampling ampltidues for noise sources from a mean = 20 nA and stdev = 10 <-- This si completely arbitirary right now.
    
    % random noise orientations
    if h.menu_ori_normal.Value==1  % Volume
        az = deg2rad(normrnd(0,180,num_noise,1));
        el = deg2rad(normrnd(0,180,num_noise,1));
        [x,y,z]=sph2cart(az,el,ones(size(az)));
        h.monte_params.cfg.source.brain_noise_ori(:,:,t) = [x,y,z];
    elseif h.menu_ori_normal.Value==2  % Cortically Constrained
        %         h.monte_params.cfg.source.brain_noise_ori(:,:,t) = h.anatomy.leadfield_eeg_cortex(h.menu_sens_montage.Value).ori(h.monte_params.cfg.source.brain_noise_idx(:,t),:);
        if ~isfield(h.anatomy.leadfield,'ori')    % cortical surface leadfields exist
            h.menu_head_model.Enable = 'on'; %h.menu_inv_headmodel.Enable = 'on';
            h.menu_head_model.BackgroundColor = [1 1 1];  h.menu_head_model.ForegroundColor = [1 1 1]*0;
            h.pwd = pwd;
            cd(fullfile(h.FieldTrip_dir,'private\'));
            %     for n=1:length(h.anatomy.leadfield_eeg_cortex); h.anatomy.leadfield_eeg_cortex(n).ori = normals(h.anatomy.mesh_cortex.pos,h.anatomy.mesh_cortex.tri,'vertex'); end% finding orientations normalized to cortical surface "cortically contrained"
            h.anatomy.leadfield.ori = normals(h.anatomy.mesh_cortex.pos,h.anatomy.mesh_cortex.tri,'vertex'); % finding orientations normalized to cortical surface "cortically contrained"
            cd(h.pwd);
            % else
            %     h.menu_head_model.Enable = 'inactive'; %h.menu_inv_headmodel.Enable = 'inactive';
            %     h.menu_head_model.Value = 1;
            %     h.menu_head_model.BackgroundColor = [1 1 1]*1;  h.menu_head_model.ForegroundColor = [1 1 1]*0;
        end
        
        try
            h.monte_params.cfg.source.brain_noise_ori(:,:,t) = h.anatomy.leadfield.ori(h.monte_params.cfg.source.brain_noise_idx(:,t),:);
        catch
            h.cfg.source.brain_noise_ori(:,:,t) = h.anatomy.leadfield.ori(h.cfg.source.brain_noise_idx(:,t),:);
        end
    end
end

