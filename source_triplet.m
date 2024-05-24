classdef source_triplet < handle
    % Class Dipole: create a dipole object for dipole fitting
    % USAGE:
    %   dip = dipole();
    %       or
    %   dip = dipole(num, mri_pos, mri_transform);
    %
    % Default Properties:

    properties  %% Accessible & Changeable
        vx_locs = nan(3,3); % positions X, Y, Z of slice positions of MRI slices of size(anatomy.mri.anatomy)
        vx_idx  = nan(1,3); % voxel indices in leadfield
        vx_ori  (3,3) double {mustBeInRange(vx_ori, -1, 1)} = [1 0 0; 0 1 0; 0 0 1]; % cartesian coordinates
        vx_ori_sph  (3,2) double {mustBeInRange(vx_ori_sph, -360, 360)} = [90 0; 0 90; 180 0]; % cartesian coordinates
        vx_amp  (1,3) double {mustBeInRange(vx_amp, 0, 1e10)} = [60 60 60]; % source amplitudes
        src_clr (3,3) double {mustBeInRange(src_clr, 0, 1)} = [0 .6 0; 0 0 1; 1 0 0]; % source colours
        sig_freqs = reshape(repmat([6 6],3, 1),[3 1 2]); % [source x tfr_roi x (start_freq end_freq)]; where tfr_roi increase with size
        sig_amp_perc = repmat(100,3,1)
        prepost_amp_perc = repmat(100,3,1)
        sig_amp_perc_std = zeros(3,1);
        prepost_amp_perc_std = zeros(3,1);
        sig_evoked_perc = repmat(50,3,1);
        prepost_evoked_perc = zeros(3,1);
        sig_durs  = zeros(3,1);
        sig_start = zeros(3,1);
        sig_win_type = zeros(3,1);
        sig_win_rise_time = zeros(3,1);
        sig_PLV_targets = zeros(3,1);
        prepost_PLV_targets = zeros(3,1);
        sig_PLI_targets = zeros(3,1);
        prepost_PLI_targets = zeros(3,1);
        sig_phase_lag = zeros(3,1);
        prepost_phase_lag = zeros(3,1);
        phase_amp_contrasts = [1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; % These are fixed --> sig_contrasts = cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
        sig_phase_amp_freq_idx = zeros(6, 1); %[0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; %(sig_contrasts x Nfreqs);  % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
        prepost_phase_amp_freq_idx = zeros(6, 1);   %(PLV_contrasts x Nfreqs);     % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
        sig_phase_amp_depth_perc = zeros(6, 1);
        prepost_phase_amp_depth_perc = zeros(6, 1);
        sig_phase_amp_depth_perc_range = zeros(6, 1);
        prepost_phase_amp_depth_perc_range = zeros(6, 1);
        sig_PLV_trials_est
        sig_PLV_evoked_est
        sig_PLI_trials_est
        sig_dPLI_trials_est
        prepost_PLV_trials_est
        prepost_PLV_evoked_est
        prepost_PLI_trials_est
        prepost_dPLI_trials_est
    end

    methods

        %% Add functions here
        %% Find Voxel indices in leadfield
        function obj = fcn_find_nearest_idx(obj, leadfield_pos)
            obj.vx_idx = find_nearest_voxel(obj.vx_locs, leadfield_pos);
        end
        %% convert ori to ori_sph
        function obj = fcn_ori2sph(obj)
            [az, el]= cart2sph(obj.vx_ori(:,1),obj.vx_ori(:,2),obj.vx_ori(:,3));
            obj.vx_ori_sph = rad2deg([az el]);
        end
        %% convert ori_sph to ori
        function obj = fcn_ori2cart(obj)
            sph = deg2rad(obj.vx_ori_sph);
            [obj.vx_ori(:,1),obj.vx_ori(:,2),obj.vx_ori(:,3)] = sph2cart( sph(:,1), sph(:,2), ones(size(obj.vx_ori_sph,1),1) );
        end
    end
end