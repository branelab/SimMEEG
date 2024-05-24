function run_source_modeling(varargin)
%% Run Source modeling
global h

h.run_inv_soln_flag = 1; % user pressed "run modeling" button
inv_soln = h.menu_inv_soln.String{h.menu_inv_soln.Value};


if strcmp(inv_soln,'Dipole')
    BL_DipoleModeling_SimMEEG(h); 
    return;
end


h.menu_invsoln_map_type.Value = 1; % need to set this to make sure peak found are based on standard P.img

% try
h.btn_3D_plot_peak_waves.Value=0;
h.next_inv_soln = length(h.inv_soln)+1;
h.current_inv_soln = h.next_inv_soln;

%% for Monte Carlo
if h.monte_carlo_flag == 1
    h.waitfor_txt.String = sprintf('Inverse Modeling using %s\n',inv_soln); drawnow;
else
    h.waitfor_panel.Visible='on';
    h.waitfor_txt.String = sprintf('Inverse Modeling using %s\n',inv_soln); drawnow;
end

%% Setting Leadfield and headmodels in h.inv_soln for future recall of saved files
if h.menu_inv_headmodel.Value==1    % Whole Brain
    if h.menu_sens_type.Value == 1 %MEG
        h.inv_soln(h.next_inv_soln).leadfield =  h.anatomy.leadfield_meg_vol(h.menu_sens_montage.Value);
        h.inv_soln(h.next_inv_soln).headmodel =  h.anatomy.headmodel_meg_vol(h.menu_sens_montage.Value);
        h.inv_soln(h.next_inv_soln).sens =  h.anatomy.sens_meg(h.menu_sens_montage.Value);
        h.inv_soln(h.current_inv_soln).headmodel_type = 'Whole Brain';
        %         h.anatomy.leadfield =  h.anatomy.leadfield_meg_vol;
        %         h.anatomy.headmodel =  h.anatomy.headmodel_meg_vol;
        %         h.anatomy.sens =  h.anatomy.sens_meg;
    elseif h.menu_sens_type.Value == 2 %EEG
        h.inv_soln(h.next_inv_soln).leadfield =  h.anatomy.leadfield_eeg_vol(h.menu_sens_montage.Value);
        h.inv_soln(h.next_inv_soln).headmodel =  h.anatomy.headmodel_eeg_vol(h.menu_sens_montage.Value);
        h.inv_soln(h.next_inv_soln).sens =  h.anatomy.sens_eeg(h.menu_sens_montage.Value);
        h.inv_soln(h.current_inv_soln).headmodel_type = 'Whole Brain';
        %         h.anatomy.leadfield =  h.anatomy.leadfield_eeg_vol;
        %         h.anatomy.headmodel =  h.anatomy.headmodel_eeg_vol;
        %         h.anatomy.sens =  h.anatomy.sens_eeg;
    end
    h.inv_soln(h.next_inv_soln).headmodel_mesh =  h.anatomy.mesh_volumes(3);

elseif h.menu_inv_headmodel.Value==2    % Cortical Surface
    if h.menu_sens_type.Value == 1 %MEG
        h.inv_soln(h.next_inv_soln).leadfield =  h.anatomy.leadfield_meg_cortex(h.menu_sens_montage.Value);
        h.inv_soln(h.next_inv_soln).headmodel =  h.anatomy.headmodel_meg_cortex(h.menu_sens_montage.Value);
        h.inv_soln(h.next_inv_soln).sens =  h.anatomy.sens_meg(h.menu_sens_montage.Value);
        h.inv_soln(h.current_inv_soln).headmodel_type = 'Cortical Surface';
        %         h.anatomy.leadfield =  h.anatomy.leadfield_meg_cortex;
        %         h.anatomy.headmodel =  h.anatomy.headmodel_meg_cortex;
        %         h.anatomy.sens =  h.anatomy.sens_meg;
    elseif h.menu_sens_type.Value == 2 %EEG

        h.inv_soln(h.next_inv_soln).leadfield =  h.anatomy.leadfield_eeg_cortex(h.menu_sens_montage.Value);
        h.inv_soln(h.next_inv_soln).headmodel =  h.anatomy.headmodel_eeg_cortex(h.menu_sens_montage.Value);
        h.inv_soln(h.next_inv_soln).sens =  h.anatomy.sens_eeg(h.menu_sens_montage.Value);
        h.inv_soln(h.current_inv_soln).headmodel_type = 'Cortical Surface';
        %         h.anatomy.leadfield =  h.anatomy.leadfield_eeg_cortex;
        %         h.anatomy.headmodel =  h.anatomy.headmodel_eeg_cortex;
        %         h.anatomy.sens =  h.anatomy.sens_eeg;
    end
    h.inv_soln(h.next_inv_soln).headmodel_mesh =  h.anatomy.mesh_volumes(4);

end
h.inv_soln(h.current_inv_soln).maxvectorori_Type =[];

h.current_3D_peak_idx =[]; h.current_3D_peak_voxels = []; h.current_norm_peak_swf =[]; h.current_peak_swf = [];
inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
h.cfg.study.bl_bmf.inside_idx = inside_idx;
mHref = []; mref = [];

%% checking variable exist
if ~isfield(h.sim_data,'sig_final'); h.sim_data.sig_final = nan(size(h.sim_data.sens_final,1), length(h.sim_data.cfg.source.vx_idx));  end

%% data
data = h.sim_data.sens_final; % [samples x channels x trials]
%% convert BrainSim data to Field Trip data format
ft_data = convert_bs2ft_data(data,h.anatomy,h.cfg);
h.ft_data = ft_data;

%% removing bad_chans from data and leadfields
if isfield(h.anatomy.sens,'good_sensors')
    if ~isempty(h.anatomy.sens.bad_sensors)
        data = data(:,h.anatomy.sens.good_sensors,:);
        h.inv_soln(h.next_inv_soln).leadfield.H=h.inv_soln(h.next_inv_soln).leadfield.H(h.anatomy.sens.good_sensors,:,:);
        fprintf('Rank = %.f\n',rank(squeeze(data(:,:,1))));
    end
end

%% getting Active & Control intervals
h.cfg.study.bl_bmf.act_int = str2num(h.edit_act_int.String); % active interval set by user
h.cfg.study.bl_bmf.ctrl_int = str2num(h.edit_ctrl_int.String);
act_s = round( (h.cfg.study.bl_bmf.act_int-h.cfg.study.lat_sim(1)) *h.cfg.study.srate);     act_s(act_s==0)=1;
ctrl_s = round( (h.cfg.study.bl_bmf.ctrl_int-h.cfg.study.lat_sim(1)) *h.cfg.study.srate);   ctrl_s(ctrl_s==0)=1;
h.cfg.study.bl_bmf.act_samps = act_s(1):act_s(end);
h.cfg.study.bl_bmf.ctrl_samps = ctrl_s(1):ctrl_s(end);
bs = round( (h.cfg.study.base_int-h.cfg.study.lat_sim(1))*h.cfg.study.srate);
h.cfg.study.h.cfg.study.base_samps = bs(1):bs(2); %           params.study.h.cfg.study.base_samps;     % -300 to 0 ms

h.cfg.study.bl_bmf.slice_orient = [1 1 1]; h.cfg.study.bl_bmf.sFaceAlpha = .25;
h.cfg.study.bl_bmf.vw_angle = [0 90];

h.inv_soln(h.current_inv_soln).params = h.cfg.study.bl_bmf;



%% Preprocesing Data for Field Trip
switch h.menu_inv_soln.String{h.menu_inv_soln.Value}
    case {'SIA' 'MIA' 'SPA'}
        h.cfg.study.bl_bmf.noise_alpha = str2num(h.edit_SPA_noise_alpha.String);
        %     h.cfg.study.bl_bmf.percCovNoise = str2num(h.edit_SPA_perc_whiten_data.String);

    case {'LCMV (FT)'  'MNE (FT)' 'sLORETA (FT)' 'eLORETA (FT)' 'dics (FT)' 'pcc (FT)' 'SAM (FT)'}   % for Field Trip Source modeling
        ft_cfg = [];
        ft_cfg.lambda         = str2num(h.edit_inv_lambda.String);
        ft_cfg.powmethod = h.menu_inv_powmethod.String{h.menu_inv_powmethod.Value}; %        = can be 'trace' or 'lambda1'
        ft_cfg.feedback = 'none'; %'         = give ft_progress indication, can be 'text', 'gui' or 'none' (default)
        ft_cfg.fixedori = 'yes'; %        = use fixed or free orientation,                   can be 'yes' or 'no'
        ft_cfg.projectnoise = 'yes'; %     = project noise estimate through filter,           can be 'yes' or 'no'
        ft_cfg.projectmom = 'yes'; %      = project the dipole moment timecourse on the direction of maximal power, can be 'yes' or 'no'
        ft_cfg.keepfilter = 'yes'; %      = remember the beamformer filter,                  can be 'yes' or 'no'
        ft_cfg.keepleadfield = 'yes'; %    = remember the forward computation,                can be 'yes' or 'no'
        ft_cfg.keepmom = 'no'; %          = remember the estimated dipole moment timeseries, can be 'yes' or 'no'
        ft_cfg.keepcov = 'no'; %          = remember the estimated dipole covariance,        can be 'yes' or 'no'
        ft_cfg.kurtosis = 'no'; %         = compute the kurtosis of the dipole timeseries,   can be 'yes' or 'no'
        % These options influence the forward computation of the leadfield
        ft_cfg.reducerank = h.menu_inv_reducerank.String{h.menu_inv_reducerank.Value}; %      = reduce the leadfield rank, can be 'no' or a number (e.g. 2)
        ft_cfg.normalize = h.menu_inv_normalize.String{h.menu_inv_normalize.Value}; %        = normalize the leadfield
        ft_cfg.normalizeparam = str2num(h.edit_inv_normalizeparam.String); %  = parameter for depth normalization (default = 0.5)
        ft_cfg.prewhiten = h.menu_inv_prewhiten.String{h.menu_inv_prewhiten.Value};  % = 'no' or 'yes', prewhiten the leadfield matrix with the noise covariance matrix C
        % ft_cfg.realfilter   = 'yes';
    case  {'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}  % for Brainstorm3 source modeling
        %% Brainstorm Parameters are gathered below when running inverse solution - see "process_inverse_2018.m" lines 541-573 and 632-649
        h.bst.gui_flag = 0; % export sens_final data to brainstorm to run Brainstorm's inverse solutions from Brainstorm

        [h] = sm_open_bst_study(h);  % opens study in brainstorm, if ProtocolName doesn't exist then it creates a new Protocol and Study.

        if isappdata(0, 'BrainstormRunning') && ~isempty(h.bst.Study)
            [h] = sm_bst_update_study(h);    % updating selected study in Brainstorm to run source modeling directly from brainstorm
            sm_bst_select_exported_data; % selecting previously exported SimMMEG data within Brainstorm if it exists

            SubjectNames = fileparts(h.bst.Study.BrainStormSubject);   % returns Subject name
            sFiles = [];
            sFiles = bst_process('CallProcess', 'process_select_files_data', sFiles, [], ...
                'subjectname',   SubjectNames, ...
                'condition',     'NewCondition', ...
                'tag',           'Avg: 01-sens_final', ...
                'includebad',    0, ...
                'includeintra',  0, ...
                'includecommon', 0);
        else
            sFiles = [];
        end

        if isempty(sFiles) || h.radio_overwrite_bst_data.Value==1  % already exported sens_final trials and created average
            [h] = sm_bst_update_study(h);    % updating selected study in Brainstorm to run source modeling directly from brainstorm
            sm_bst_select_exported_data; % selecting previously exported SimMMEG data within Brainstorm if it exists
            sm_export2brainstorm; % open Brainstorm and create database

        else
            fprintf('Using existing data in Brainstorm Subject: %s/NewCondition \n',SubjectNames)
        end

    case {'sMCMV' 'bRAPBeam' 'TrapMUSIC'}   % Alex Moiseev's sub-space multi-source scalar beamformer and TRAP-MUSIC
        if h.menu_inv_headmodel.Value==2   % Moiseev's MCMV & TRAP-MUSIC Inv Solutions
            msgbox('Only Volume Head Models are allowed for Moiseev Beamformers');
            return
        else
            % making transformation grid for Moiseev's programs
            if ~isfield(h.anatomy,'moiseev') % only Volume leadfields allowed
                vx_pos = h.inv_soln(h.next_inv_soln).leadfield.voxel_pos;

                %% creating grid "vx_grid" as per Alex's description "nVoxels = nX * nY * nZ, and the ordering is such that the last index (i.e. Z direction) changes most fast, and the first (i.e. X) - most slow"
                vx_res = diff(vx_pos(:,1)); vx_res = min(abs(vx_res(vx_res~=0)));   % voxel resolution
                % making 3D grid
                xg = min(vx_pos(:,1))-vx_res:vx_res:max(vx_pos(:,1))+vx_res;
                yg = min(vx_pos(:,2))-vx_res:vx_res:max(vx_pos(:,2))+vx_res;
                zg = min(vx_pos(:,3))-vx_res:vx_res:max(vx_pos(:,3))+vx_res;

                %% Moiseev's programs' (C++) box order
                v=0; h.anatomy.moiseev.vx_grid=[];
                for d1=1:length(xg)
                    for d2=1:length(yg)
                        for d3=1:length(zg)
                            v=v+1;
                            h.anatomy.moiseev.vx_grid(v,:) = [xg(d1) yg(d2) zg(d3)];
                        end
                    end
                end
                h.anatomy.moiseev.dims = [length(xg) length(yg) length(zg)];
                % iV = idx3Dto1D(vx_grid, dims, 1)
                h.anatomy.moiseev.v_idx = find_nearest_voxel(vx_pos,h.anatomy.moiseev.vx_grid);
                grid_idx = zeros(size(h.anatomy.moiseev.vx_grid,1),1);    % inside brain indices of leadfields
                grid_idx(h.anatomy.moiseev.v_idx)=1; h.anatomy.moiseev.lstFlag = grid_idx;
                h.anatomy.moiseev.inside_idx = find(h.anatomy.moiseev.lstFlag==1);
                %                 figure(998); clf; hold on; scatter3(h.anatomy.moiseev.vx_grid(grid_idx==0,1),h.anatomy.moiseev.vx_grid(grid_idx==0,2),h.anatomy.moiseev.vx_grid(grid_idx==0,3),'k.');
                %         scatter3(h.anatomy.moiseev.vx_grid(grid_idx==1,1),h.anatomy.moiseev.vx_grid(grid_idx==1,2),h.anatomy.moiseev.vx_grid(grid_idx==1,3),'ro');


                %% matlab box order seed in SimMEEG
                v=0; h.anatomy.moiseev.simmeeg_grid=[];
                for d1=1:length(zg)
                    for d2=1:length(yg)
                        for d3=1:length(xg)
                            v=v+1;
                            h.anatomy.moiseev.simmeeg_grid(v,:) = [xg(d3) yg(d2) zg(d1)];
                        end
                    end
                end
                h.anatomy.moiseev.simmeeg_idx = find_nearest_voxel(vx_pos,h.anatomy.moiseev.simmeeg_grid);
                grid_idx = zeros(size(h.anatomy.moiseev.simmeeg_grid,1),1);    % inside brain indices of leadfields
                grid_idx(h.anatomy.moiseev.simmeeg_idx)=1; h.anatomy.moiseev.simmeeg_inside = grid_idx;
                h.anatomy.moiseev.simmeeg_inside_idx = find(h.anatomy.moiseev.simmeeg_inside==1);

                %         figure(999); clf; hold on;
                %         scatter3(h.anatomy.moiseev.simmeeg_grid(h.anatomy.moiseev.simmeeg_inside==0,1),h.anatomy.moiseev.simmeeg_grid(h.anatomy.moiseev.simmeeg_inside==0,2),h.anatomy.moiseev.simmeeg_grid(h.anatomy.moiseev.simmeeg_inside==0,3),'k.');
                %         scatter3(h.anatomy.moiseev.simmeeg_grid(h.anatomy.moiseev.simmeeg_inside==1,1),h.anatomy.moiseev.simmeeg_grid(h.anatomy.moiseev.simmeeg_inside==1,2),h.anatomy.moiseev.simmeeg_grid(h.anatomy.moiseev.simmeeg_inside==1,3),'ro');

            else
            end
            if h.menu_inv_headmodel.Value==1
                %% Covariance matrix
                act_samps = h.cfg.study.bl_bmf.act_samps;
                ctrl_samps = h.cfg.study.bl_bmf.ctrl_samps;
                xs=round(size(ctrl_samps,2)/2);
                %% Covariance for Active vs. Control interval
                [R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(data,act_samps,ctrl_samps(xs+1:end));
                %% Covariance for Noise estimate from control samples
                % splitting ctrl interval in half to calculate null distribution of noise
                [nR,nN,nRbar,nNbar,nRinv,nNinv]=BRANELab_calc_cov(data,ctrl_samps(1:xs),ctrl_samps(xs+1:end));
                %     [nR,nN,nRbar,nNbar,nRinv,nNinv]=BRANELab_calc_cov(h.sim_data.sens_noise_scaled,ctrl_samps(1:xs),ctrl_samps(xs+1:end));
                %% reconfigure leadfields indices to match the moissev grid data
                %     H = permute(h.inv_soln(h.current_inv_soln).leadfield.H,[3 2 1]); arrH = nan(size(h.anatomy.moiseev.lstFlag,1),size(H,2),size(H,3));
                H = permute(h.inv_soln(h.next_inv_soln).leadfield.H,[3 2 1]); arrH = nan(size(h.anatomy.moiseev.lstFlag,1),size(H,2),size(H,3));
                arrH(h.anatomy.moiseev.v_idx,:,:) = H;
                % other parameters
                maxSrc = str2num(h.edit_SPA_max_sources.String);   % Maximum number of r sources to be found
                gap = 3;    % distance between voxels for finding peaks (should be eual to 15 mm) <-- imbedded function in Moiseev code
            else
                msgbox('Only Volume Head Models are allowed for Moiseev Beamformers');
                return
            end

        end
end

%% %%%%% Inverse Solutions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.tt0 = tic; % time to compute inv_soln
switch inv_soln
    case 'SPA'
        h.cfg.study.bl_bmf.loc_flag = h.menu_SPA_loc_flag.Value-1;
        h.cfg.study.bl_bmf.noise_alpha = str2num(h.edit_SPA_noise_alpha.String);
        [SPA]=BRANElab_LCMV_beamformer_SPA(h.inv_soln(h.next_inv_soln).leadfield.H,data,...
            h.cfg.study.bl_bmf.act_samps,h.cfg.study.bl_bmf.ctrl_samps,h.cfg.study.bl_bmf.loc_flag,h.cfg.study.bl_bmf.noise_alpha);
        h.inv_soln(h.next_inv_soln).Type = inv_soln;
        h.inv_soln(h.next_inv_soln).soln = SPA;
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.menu_SPA_loc_flag.String{h.menu_SPA_loc_flag.Value} ,h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case 'SIA'
        h.cfg.study.bl_bmf.loc_flag = h.menu_SPA_loc_flag.Value-1;
        h.cfg.study.bl_bmf.plot_flag=0; h.cfg.study.bl_bmf.text_flag=1;
        h.cfg.study.bl_bmf.perc_crit = 2;
        max_sources = str2num(h.edit_SPA_max_sources.String); %floor(size(h.inv_soln(h.next_inv_soln).leadfield.H,1)/5);   % 1 dipole source has 5 degrees of freedom (3 spatial and 2 orientations), thus dividing number of channels/sensors by 5
        [SIA]=BRANElab_MCMV_beamformer_SIA(h.inv_soln(h.next_inv_soln).leadfield.H,[],[],data,...
            h.cfg.study.bl_bmf.act_samps,h.cfg.study.bl_bmf.ctrl_samps,...
            h.inv_soln(h.next_inv_soln).leadfield.voxel_pos,h.cfg.study.bl_bmf.loc_flag,h.cfg.study.bl_bmf.plot_flag,h.cfg.study.bl_bmf.perc_crit,h.cfg.study.bl_bmf.noise_alpha,max_sources);
        h.inv_soln(h.next_inv_soln).Type = inv_soln;
        h.inv_soln(h.next_inv_soln).soln = SIA;
        h.inv_soln(h.next_inv_soln).soln.P.img_org = h.inv_soln(h.next_inv_soln).soln.P.img;
        h.inv_soln(h.next_inv_soln).soln.wts_org = h.inv_soln(h.next_inv_soln).soln.wts;
        h.inv_soln(h.next_inv_soln).soln.ori_org = h.inv_soln(h.next_inv_soln).soln.ori;
        %         h.inv_soln(h.next_inv_soln).soln.P.img = h.inv_soln(h.next_inv_soln).soln.P.nulled_img; % shifting to nulled img
        h.inv_soln(h.next_inv_soln).soln.wts = h.inv_soln(h.next_inv_soln).soln.nulled_wts;
        h.inv_soln(h.next_inv_soln).soln.ori = h.inv_soln(h.next_inv_soln).soln.nulled_ori;
        h.inv_soln(h.next_inv_soln).soln.null_thresh(end+1) = 0;    % making this 0 so that all peaks are found during plotting.

        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.menu_SPA_loc_flag.String{h.menu_SPA_loc_flag.Value} ,h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case 'MIA'
        h.cfg.study.bl_bmf.loc_flag = h.menu_SPA_loc_flag.Value-1;
        h.cfg.study.bl_bmf.plot_flag=0; h.cfg.study.bl_bmf.text_flag=1;
        h.cfg.study.bl_bmf.perc_crit = 2;
        max_sources = str2num(h.edit_SPA_max_sources.String); %floor(size(h.inv_soln(h.next_inv_soln).leadfield.H,1)/5);   % 1 dipole source has 5 degrees of freedom (3 spatial and 2 orientations), thus dividing number of channels/sensors by 5
        [MIA]=BRANElab_MCMV_beamformer_MIA(h.inv_soln(h.next_inv_soln).leadfield.H,[],[],data,...
            h.cfg.study.bl_bmf.act_samps,h.cfg.study.bl_bmf.ctrl_samps,...
            h.inv_soln(h.next_inv_soln).leadfield.voxel_pos,h.cfg.study.bl_bmf.loc_flag,h.cfg.study.bl_bmf.plot_flag,h.cfg.study.bl_bmf.perc_crit,h.cfg.study.bl_bmf.noise_alpha,h.cfg.study.bl_bmf.text_flag,h.anatomy,max_sources);
        h.inv_soln(h.next_inv_soln).Type = inv_soln;
        h.inv_soln(h.next_inv_soln).soln = MIA;
        h.inv_soln(h.next_inv_soln).soln.P.img_org = h.inv_soln(h.next_inv_soln).soln.P.img;
        h.inv_soln(h.next_inv_soln).soln.wts_org = h.inv_soln(h.next_inv_soln).soln.wts;
        h.inv_soln(h.next_inv_soln).soln.ori_org = h.inv_soln(h.next_inv_soln).soln.ori;
        h.inv_soln(h.next_inv_soln).soln.wts = h.inv_soln(h.next_inv_soln).soln.nulled_wts;
        %         h.inv_soln(h.next_inv_soln).soln.P.img = h.inv_soln(h.next_inv_soln).soln.P.nulled_img; % shifting to nulled img
        h.inv_soln(h.next_inv_soln).soln.ori = h.inv_soln(h.next_inv_soln).soln.nulled_ori;
        h.inv_soln(h.next_inv_soln).soln.null_thresh(end+1) = 0;    % making this 0 so that all peaks are found during plotting.

        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.menu_SPA_loc_flag.String{h.menu_SPA_loc_flag.Value} ,h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case {'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}
        [h] = sm_bst_update_study(h);
        % --- Run within SimMEEG --- BUT doesn't return same MNE/sLORETA solutions as if running within Brainstorm - Not sure why but something to do with the Covariance matrix calculations.
        %% Calculate Noise and Data Covariances --> BST uses averaged (evoked Data) to calculate covariances
        %       %% testing to see if we can get same as brainstorm
        %         x = load('C:\BRANELab\matlab_progs\general_progs\brainstorm3\inv_testing_bst_data.mat');
        %         data = permute(repmat(x.avg_data.F,[1 1 90]),[2 1 3]);
        %
        sm_bst_noisecov(data);  % create Noise and Data covariances from Brainstorm function (changed by ATH only for entering data - calculations are those from Brainstorm)
        % replacing with BRANE Lab's coavraince to see what is happening with Brainstorm's LCMV that is not returning m,uch that is meaningful.
        %         [R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(data,h.cfg.study.bl_bmf.act_samps,h.cfg.study.bl_bmf.ctrl_samps);
        %         h.inv_soln(h.current_inv_soln).soln.bst_cov.DataCovMat.NoiseCov
        %
        bst_cov = h.inv_soln(h.current_inv_soln).soln.bst_cov;
        %         bst_cov.NoiseCovMat = noise_cov;
        %         bst_cov.DataCovMat = data_cov;
        %% ChannelMat - creating BST ChannelMat format
        ChannelMat.Channel = struct('Name','','Comment',[],'Type','','Group',[],'Loc',[],'Orient',[],'Weight',[]);
        ChannelMat.Projector = [];

        for v=1:length(h.anatomy.sens.label)
            ChannelMat.Channel(v).Name = h.anatomy.sens.label{v};
            if any(strcmpi(h.anatomy.sens.type, {'EEG'}))
                ChannelMat.Channel(v).Type = 'EEG';
            elseif any(strcmpi(h.anatomy.sens.type, {'MEG'}))
                ChannelMat.Channel(v).Type = 'MEG';
            else
                ChannelMat.Channel(v).Type = h.anatomy.sens.type;
            end
            ChannelMat.Channel(v).Loc = h.anatomy.sens.chanpos(v,:)*1e-3;  % convert back to meters
            ChannelMat.Channel(v).Loc = h.anatomy.sens.chanpos(v,:)*1e-3;  % convert back to meters
            ChannelMat.Channel(v).Weight = 1;
        end

        GoodChannel = h.anatomy.sens.good_sensors;
        ChannelFlag = zeros(length(h.anatomy.sens.label),1);
        ChannelFlag(GoodChannel) = 1;
        %% convert leadfields into BST's "HeadModel"
        HeadModel = struct('MEGMethod','','EEGMethod','','ECOGMethod','','SEEGMethod','','Gain','','Comment','',...
            'HeadModelType','','GridLoc',[],'GridOrient','','GridAtlas','','GridOptions','','SurfaceFile','','Param','');
        if h.menu_sens_type.Value ==1 % MEG
            HeadModel.MEGMethod = h.inv_soln(h.current_inv_soln).headmodel.type;
        elseif h.menu_sens_type.Value ==2 % EEG
            HeadModel.EEGMethod = h.inv_soln(h.current_inv_soln).headmodel.type;
        else
            HeadModel.MEGMethod = h.inv_soln(h.current_inv_soln).headmodel.type;
            HeadModel.EEGMethod = h.inv_soln(h.current_inv_soln).headmodel.type;
        end
        %% converting leadfield to BST's vector leadfield "Gain" matrix
        HeadModel.Gain = cell2mat(h.inv_soln(h.current_inv_soln).leadfield.leadfield(h.inv_soln(h.current_inv_soln).leadfield.inside==1));
        HeadModel.Comment = h.menu_head_model.String{h.menu_head_model.Value};
        HeadModel.Type = h.menu_head_model.String{h.menu_head_model.Value};
        HeadModel.GridLoc = h.inv_soln(h.current_inv_soln).leadfield.pos/1000; % converting grid locs from 'mm' to 'm'
        %         HeadModel.GridOrient = h.inv_soln(h.current_inv_soln).leadfield./1000; % converting grid locs from 'mm' to 'm'
        HeadModel.History(1) = {datetime};
        HeadModel.History(2) = {'compute'};
        HeadModel.History(3) = {'Converted from Field Trip format to Brainstorm Format within SimMEEG'};
        %% HeadModel - Need to adjust headmodel as what is done in BST's "process_inverse_2018.m" for average ref
        HeadModel.Gain = HeadModel.Gain(GoodChannel, :);
        % Apply average reference: separately SEEG, ECOG, EEG
        if any(ismember(unique({ChannelMat.Channel.Type}), {'EEG','ECOG','SEEG'})) && isappdata(0, 'BrainstormRunning')
            % Create average reference montage
            sMontage = panel_montage('GetMontageAvgRef', [], ChannelMat.Channel(GoodChannel), ChannelFlag(GoodChannel), 0);
            HeadModel.Gain = sMontage.Matrix * HeadModel.Gain;
            % Apply average reference operator on both sides of the noise covariance matrix
            bst_cov.NoiseCovMat.NoiseCov(GoodChannel, GoodChannel) = sMontage.Matrix * bst_cov.NoiseCovMat.NoiseCov(GoodChannel, GoodChannel) * sMontage.Matrix';
            if ~isempty(bst_cov.DataCovMat)
                bst_cov.DataCovMat.NoiseCov(GoodChannel, GoodChannel) = sMontage.Matrix * bst_cov.DataCovMat.NoiseCov(GoodChannel, GoodChannel) * sMontage.Matrix';
            end
        end
        %% create OPTIONS --> BST style
        %% Creating "OPTIONS" for running Brainstorm's "bst_process_inverse_2018.m"
        opt = struct('InverseMethod','','InverseMeasure','','SourceOrient','','DataTypes','','Comment','','DisplayMessages',1,...
            'ComputeKernel',1,'Loose',0.2,'UseDepth',[],'WeightExp',[],'WeightLimit',[],'NoiseMethod','','NoiseReg','','SnrMethod','',...
            'SnrRms','','SnrFixed','','NoiseCovMat','','DataCovMat','','ChannelTypes','','FunctionName','');
        switch inv_soln
            case 'LCMV (BST)'
                opt.InverseMethod = 'lcmv';
                opt.InverseMeasure = 'nai'; % as of Oct-2022, Brainstorm only has the NAI option enabled
                opt.DataTypes = h.menu_sens_type.String(h.menu_sens_type.Value);
                opt.Comment = sprintf('PNAI: %s(Unconstr)',opt.DataTypes{:});
                opt.UseDepth = 0;
                opt.WeightExp = 0;
                opt.WeightLimit = 0;
                opt.FunctionName = 'lcmvp';
                opt.SnrMethod = 'rms';
                opt.SnrRms = 1.0000e-06; % not sure where this comes from

            case {'MNE (BST)' 'sLORETA (BST)'}
                opt.InverseMethod = 'minnorm';
                opt.DataTypes = h.menu_sens_type.String(h.menu_sens_type.Value);
                switch inv_soln
                    case 'MNE (BST)'
                        opt.InverseMeasure = 'amplitude';
                        opt.Comment = sprintf('MN: %s(Unconstr)',opt.DataTypes{:});
                        opt.UseDepth = h.radio_inv_bst_depth_weight.Value;
                        opt.FunctionName = 'mn';
                    case 'sLORETA (BST)'
                        opt.InverseMeasure = 'sloreta';
                        opt.Comment = sprintf('sLORETA: %s(Unconstr)',opt.DataTypes{:});
                        opt.UseDepth = 0;
                        opt.FunctionName = 'sloreta';
                end
                opt.WeightExp = str2num(h.edit_inv_bst_depth_order.String);
                opt.WeightLimit = str2num(h.edit_inv_bst_depth_max.String);
                opt.SnrMethod = 'fixed';
                opt.SnrRms = 1.0000e-06; % not sure where this comes from
        end
        opt.SourceOrient = {'free'}; % SimMEEG currently only allows for this option
        opt.NoiseMethod = h.cfg.study.bst_params.inv_NoiseMethod{h.menu_inv_bst_reg_noise_type.Value};
        if h.menu_inv_bst_reg_noise_type.Value==1; opt.NoiseReg = str2double(h.edit_inv_bst_reg_noise_factor.String); else opt.NoiseReg = 0; end
        opt.SnrFixed = 3;  % fixed to 3 dipole orientations
        %         opt.NoiseCovMat = h.inv_soln(h.current_inv_soln).soln.bst_cov.NoiseCovMat;
        %         opt.DataCovMat = h.inv_soln(h.current_inv_soln).soln.bst_cov.DataCovMat;
        opt.NoiseCovMat = bst_cov.NoiseCovMat;
        opt.DataCovMat = bst_cov.DataCovMat;
        opt.ChannelTypes = repmat(h.menu_sens_type.String(h.menu_sens_type.Value),1,size(opt.NoiseCovMat.NoiseCov,1));
        %% adjust covariance based on BST functions
        % NoiseCov: keep only the good channels
        opt.NoiseCovMat = bst_cov.NoiseCovMat;
        opt.NoiseCovMat.NoiseCov = opt.NoiseCovMat.NoiseCov(GoodChannel, GoodChannel);
        if isfield(opt.NoiseCovMat, 'FourthMoment') && ~isempty(opt.NoiseCovMat.FourthMoment)
            opt.NoiseCovMat.FourthMoment = opt.NoiseCovMat.FourthMoment(GoodChannel, GoodChannel);
        end
        if isfield(opt.NoiseCovMat, 'nSamples') && ~isempty(opt.NoiseCovMat.nSamples)
            opt.NoiseCovMat.nSamples = opt.NoiseCovMat.nSamples(GoodChannel, GoodChannel);
        end
        % DataCov: keep only the good channels
        if ~isempty(bst_cov.DataCovMat)
            opt.DataCovMat = bst_cov.DataCovMat;
            opt.DataCovMat.NoiseCov     = opt.DataCovMat.NoiseCov(GoodChannel, GoodChannel);
            opt.DataCovMat.FourthMoment = opt.DataCovMat.FourthMoment(GoodChannel, GoodChannel);
            opt.DataCovMat.nSamples     = opt.DataCovMat.nSamples(GoodChannel, GoodChannel);
        end
        % Get channels types
        opt.ChannelTypes = {ChannelMat.Channel(GoodChannel).Type};
        %% testing to see if we can get same as brainstorm
        %         x = load('C:\BRANELab\matlab_progs\general_progs\brainstorm3\inv_testing_bst_data.mat');
        %         opt.NoiseCovMat = x.lcmv.Options.NoiseCovMat;
        %         opt.DataCovMat = x.lcmv.Options.DataCovMat;
        %         opt=x.lcmv.Options;

        %% Testing with inputs directly from Brainstorm. Results from "bst_inverse_linear_2018.m" are close but not exact even when using all same - not sure why?
        %%      might have to do with something Brainstorm does with calling in params from BST GUI while running "bst_inverse_linear_2018.m" but A.Herdman cannot figure out where or how.
        %%      Caution! using current setting the way they are so results are likely to vary slightly between those calulated in Brainstorm vs. SimMEEG and
        %         load C:\BRANELab\Software\SimMEEG\nCov
        %         opt.NoiseCovMat = nCov;
        %         opt.DataCovMat = dCov;
        %         load C:\BRANELab\Software\SimMEEG\hdm3;
        %         load C:\BRANELab\Software\SimMEEG\lcmv3;
        %         opt = lcmv3.Options;
        %         HeadModel = hdm3;

        %% Running inverse solution
        %         load mne;
        %         opt = OPTIONS

        %% Tony's code re-written from Brainstorm <-- Not returning same ImagingKernel due to covariance matrix calculations above.
        % [BST_inv, opt] = bst_inverse_linear_2018(HeadModel, opt);

        %% Running source modeling directly from brainstorm - because Tony wasn't able to figure out how covariance calculations were being differently calculated and modified within brainstorm
        [BST_inv, opt] = sm_bst_source_modeling;


        %         BST_inv = x.lcmv;

        BST_inv.ImageGridAmp = BST_inv.ImagingKernel * squeeze(nanmean(data,3))';

        %         load C:\BRANELab\Software\SimMEEG\lcmv4;
        %         BST_inv.ImageGridAmp = lcmv4.ImageGridAmp;

        %         dims = size(BST_inv.ImageGridAmp);
        %         BST_inv.swf = permute(reshape(BST_inv.ImageGridAmp,[3 dims(1)/3 dims(2)]),[3 2 1]); % [sens x samps x 3dipole_vecotrs]
        dims = size(BST_inv.ImagingKernel);
        BST_inv.wts = permute(reshape(BST_inv.ImagingKernel,[3 dims(1)/3 dims(2)]),[3 2 1]); % [sens x samps x 3dipole_vecotrs]
        %
        %         %% calculating image power from filter weights;
        %         swf = []; swf_rms=[]; swf_pwr=[];BST_inv.noise_pwr=[];
        %         for ox=1:3
        %             swf(:,ox,:)=squeeze(nanmean(data,3))*squeeze(BST_inv.wts(:,:,ox));
        %             %    swf_snr(:,ox,:)=20*log10(bsxfun(@rdivide,swf(:,ox,:),nanmean(swf(h.cfg.study.base_samps,ox,:),1)));
        %             %   swf_snr(:,ox,:)=(bsxfun(@rdivide,swf(:,ox,:),nanmean(swf(h.cfg.study.base_samps,ox,:),1))); % SNR in percent from baseline
        %             %            swf_snr(:,ox,:)=(bsxfun(@rdivide,swf(:,ox,:),nanstd(swf(h.cfg.study.base_samps,ox,:),1))); % normalized to baseline
        %             swf_pwr(ox,:)=rms(swf(h.cfg.study.bl_bmf.act_samps,ox,:),1)./rms(swf(h.cfg.study.bl_bmf.ctrl_samps,ox,:),1);
        %             swf_rms(ox,:)=rms(swf(h.cfg.study.bl_bmf.act_samps,ox,:),1);
        %             %            BST_inv.noise_pwr(ox,:)=rms(swf(h.cfg.study.bl_bmf.ctrl_samps,ox,:),1)./rms(swf(h.cfg.study.bl_bmf.ctrl_samps,ox,:),1);
        %         end


        %% Vector Summation
        iVertSource = 1:size(BST_inv.ImageGridAmp,1);
        vec_act = BST_inv.ImageGridAmp(:,h.cfg.study.bl_bmf.act_samps);
        vec_ctrl = BST_inv.ImageGridAmp(:,h.cfg.study.bl_bmf.ctrl_samps);
        idx1 = 1:3:length(iVertSource); idx2 = 2:3:length(iVertSource); idx3 = 3:3:length(iVertSource);

        switch h.menu_inv_bst_maxvectorori.String{h.menu_inv_bst_maxvectorori.Value}
            case 'RMS'
                vec_img  = rms(vec_act,2) ./ rms(vec_ctrl,2);   % RMS and simple division by ctrl interval values to yield signal-to-noise image
                %                 vec_img  = rms(vec_act,2);   % raw
                img = squeeze(sqrt(vec_img(idx1,:,:,:).^2 + vec_img(idx2,:,:,:).^2 + vec_img(idx3,:,:,:).^2)); % RMS

                [BST_inv.P.img, max_idx] = max(abs(img),[],2); % maximum across time
                % Approximating orientations using RMS img
                x = nan(length(idx1),1)'; y=x; z=x;
                for v=1:length(idx1)
                    x(v)=vec_img(idx1(v),max_idx(v))./BST_inv.P.img(v);
                    y(v)=vec_img(idx2(v),max_idx(v))./BST_inv.P.img(v);
                    z(v)=vec_img(idx3(v),max_idx(v))./BST_inv.P.img(v);
                end

                BST_inv.ori=[x;y;z]';
                BST_inv.vectori_type = 'rms';

            case 'Max'    % vector dipole orientation with maximum img power between active and control interval
                %                 img = squeeze(max(abs(BST_inv.ImageGridAmp(:,h.cfg.study.bl_bmf.act_samps)),[],2));
                %                 vec_img = reshape(img,[3 size(img,1)/3]); % [3dipole_vecotrs x grid_locs]
                %% took same functions as in Brainstorm "bst_source_orient.m" to get image values
                %                 iVertSource = 1:size(img,1);
                %                 idx1 = 1:3:length(iVertSource); idx2 = 2:3:length(iVertSource); idx3 = 3:3:length(iVertSource);
                %                 [BST_inv.P.img, BST_inv.max_ori] = max(abs(squeeze(cat(2, img(idx1), img(idx2), img(idx3)))), [], 2);
                %                 vec_img  = max(abs(vec_act),[],2) ./ max(abs(vec_ctrl),[],2);   % RMS and simple division by ctrl interval values to yield signal-to-noise image
                vec_img  = max(abs(vec_act),[],2);   % RMS and simple division by ctrl interval values to yield signal-to-noise image
                [BST_inv.P.img, BST_inv.max_ori] = max(abs(squeeze(cat(2, vec_img(idx1), vec_img(idx2), vec_img(idx3)))), [], 2);
                ori = zeros(size(BST_inv.max_ori,1),3); for v=1:size(ori,1); ori(v,BST_inv.max_ori(v))=1; end
                BST_inv.ori = ori;
                BST_inv.vectori_type = 'maxabs';
        end

        nd=sort(BST_inv.P.img); BST_inv.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;
        % sotring variables in SimMEEG format
        h.inv_soln(h.next_inv_soln).Type = inv_soln;
        h.inv_soln(h.next_inv_soln).soln = BST_inv;
        % adding covariance in BRANELab format
        h.inv_soln(h.next_inv_soln).soln.RNcov.R = bst_cov.DataCovMat.NoiseCov;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Rbar = bst_cov.DataCovMat.FourthMoment;
        h.inv_soln(h.next_inv_soln).soln.RNcov.N = bst_cov.NoiseCovMat.NoiseCov;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Nbar = bst_cov.NoiseCovMat.FourthMoment;

        % adding name
        h.inv_soln(h.next_inv_soln).maxvectorori_Type = h.menu_inv_bst_maxvectorori.String{h.menu_inv_bst_maxvectorori.Value};
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.menu_inv_bst_maxvectorori.String{h.menu_inv_bst_maxvectorori.Value}, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case 'LCMV (FT)'
        act_int = str2num(h.edit_act_int.String);
        ctrl_int = str2num(h.edit_ctrl_int.String);
        inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);

        %% averaging data and getting covariance matrix
        cfg = [];
        cfg.channel = h.anatomy.sens.label(h.anatomy.sens.good_sensors);
        cfg.baseline = h.sim_data.cfg.study.base_int;
        ft_data = ft_timelockbaseline(cfg, ft_data);
        cfg = [];
        cfg.channel = h.anatomy.sens.label(h.anatomy.sens.good_sensors);
        cfg.covariance='yes';
        cfg.covariancewindow = [ctrl_int(1) act_int(2)]; % calculating covariance across entire ctrl and act interval just as stated in Field Trip's Tutorial
        avg_data = ft_timelockanalysis(cfg,ft_data);

        %% Calculating LCMV
        cfg = [];
        cfg.lcmv = ft_cfg;
        cfg.keepleadfield = cfg.lcmv.keepleadfield;
        cfg.reducerank = cfg.lcmv.reducerank;
        cfg.normalize = cfg.lcmv.normalize;
        cfg.normalizeparam = cfg.lcmv.normalizeparam;
        cfg.lcmv=rmfield(cfg.lcmv,'keepleadfield');  cfg.lcmv=rmfield(cfg.lcmv,'reducerank');
        cfg.lcmv=rmfield(cfg.lcmv,'normalize');  cfg.lcmv=rmfield(cfg.lcmv,'normalizeparam');

        cfg.weightnorm          = h.menu_inv_weightnorm.String{h.menu_inv_weightnorm.Value}; % ''; %'nai'; %'unitnoisegain'; %'unitgain'; %'arraygain'; % 'nai';  % based on equation 4.47 from Sekihara & Nagarajan (2008)
        cfg.lcmv.fixedori       = 'yes'; %  'yes' = get single orientation using svd
        cfg.method              = 'lcmv';
        cfg.grad                = h.inv_soln(h.next_inv_soln).sens;
        cfg.sourcemodel         = h.inv_soln(h.next_inv_soln).leadfield;
        cfg.headmodel           = h.inv_soln(h.next_inv_soln).headmodel;

        LCMV = ft_sourceanalysis(cfg, avg_data);
        LCMV = ft_sourcedescriptives([], LCMV); % to get the neural-activity-index

        h.LCMV = LCMV;

        %% adding covariance in BRANELab format
        h.inv_soln(h.next_inv_soln).soln.RNcov.R = avg_data.cov;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Rbar = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.N = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.Nbar = [];

        %% creating LCMV power image and estimating orientations
        LCMV.P.img = LCMV.avg.nai; %LCMV.avg.pow(inside_idx)./LCMV.avg.noise(inside_idx);   % simple ratio

        nd=sort(LCMV.P.img); LCMV.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;
        ori=cell2mat(LCMV.avg.ori(inside_idx)); LCMV.ori=ori';
        wts = cell2mat(LCMV.avg.filter(inside_idx)); LCMV.wts = wts';

        %                 img=LCMV.P.img; ori=LCMV.ori; null_thresh=max(img)*.55;
        %                 min_max=[min(img) max(img)*.95];  vol_types=1;
        %                 seed_idx = 1:3;  ln_wdth = 1; ln_wdth2 = 2;
        %                 figure(1004); clf; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT_new(img,ori,null_thresh,8,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,jet(255),min_max,h.anatomy.vol.bnd(3),[],...
        %                     h.cfg.study.bl_bmf.vw_angle,1,vol_types,h.inv_soln(h.current_inv_soln).leadfield.pos,inside_idx,h.cfg.study.bl_bmf.slice_orient,h.cfg.study.bl_bmf.sFaceAlpha);
        %                 hold on; mrk_size=150; s1=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'k+','sizedata',mrk_size,'linewidth',ln_wdth); s2=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'ks','sizedata',mrk_size,'linewidth',ln_wdth2);
        %                 title('FT LCMV'); colorbar; caxis(min_max);

        % removing 'avg' to reduce storage space
        LCMV = rmfield(LCMV,{'avg'}); %MNE.avg = rmfield(MNE.avg,{'mom','filter','noisecov'});

        h.inv_soln(h.next_inv_soln).Type = inv_soln;
        h.inv_soln(h.next_inv_soln).soln = LCMV;
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, cfg.weightnorm ,h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case 'eLORETA (FT)'
        act_int = str2num(h.edit_act_int.String);
        ctrl_int = str2num(h.edit_ctrl_int.String);
        inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);

        %% averaging data and getting covariance matrix
        cfg = []; cfg.baseline = h.sim_data.cfg.study.base_int;
        cfg.channel = h.anatomy.sens.label(h.anatomy.sens.good_sensors);
        ft_data = ft_timelockbaseline(cfg, ft_data);
        cfg = [];
        cfg.channel = h.anatomy.sens.label(h.anatomy.sens.good_sensors);
        cfg.covariance='yes';
        cfg.covariancewindow = [ctrl_int(1) act_int(2)]; % calculating covariance across entire ctrl and act interval just as stated in Field Trip's Tutorial
        avg_data = ft_timelockanalysis(cfg,ft_data);

        %% FT eLORETA
        cfg=[];
        cfg.method = 'eloreta';
        cfg.eloreta = ft_cfg; % using same additional options for all inverse solns
        cfg.keepleadfield = cfg.eloreta.keepleadfield;
        cfg.reducerank = cfg.eloreta.reducerank;
        cfg.normalize = cfg.eloreta.normalize;
        cfg.normalizeparam = cfg.eloreta.normalizeparam;
        cfg.eloreta=rmfield(cfg.eloreta,'keepleadfield');  cfg.eloreta=rmfield(cfg.eloreta,'reducerank');
        cfg.eloreta=rmfield(cfg.eloreta,'normalize');  cfg.eloreta=rmfield(cfg.eloreta,'normalizeparam');

        cfg.grad                = h.inv_soln(h.next_inv_soln).sens;
        cfg.sourcemodel         = h.inv_soln(h.next_inv_soln).leadfield;
        cfg.headmodel           = h.inv_soln(h.next_inv_soln).headmodel;
        eloreta  = ft_sourceanalysis(cfg,avg_data);

        %% calculating image power from filter weights;
        base_samps=h.sim_data.cfg.study.base_samps;     % -300 to 0 ms
        act_samps=h.cfg.study.bl_bmf.act_samps;
        ctrl_samps=h.cfg.study.bl_bmf.ctrl_samps;
        ctrl_samps1 = ctrl_samps(1:round(length(ctrl_samps)/2));
        ctrl_samps2 = ctrl_samps(round(length(ctrl_samps)/2)+1:end);

        y=cell2mat(eloreta.avg.filter); eloreta.wts=permute(reshape(y,[size(y,1) size(y,2)/length(inside_idx) length(inside_idx)]),[3 1 2]);
        for ox = 1:size(eloreta.wts,2)
            swf(:,ox,:)=avg_data.avg'*squeeze(eloreta.wts(:,ox,:))';
            swf_pwr(ox,:)=rms(swf(h.cfg.study.bl_bmf.act_samps,ox,:),1)./rms(swf(h.cfg.study.bl_bmf.ctrl_samps,ox,:),1);
        end
        %         eloreta.wts = y';
        eloreta.wts = permute(eloreta.wts,[3 1 2]);
        switch h.menu_inv_maxvectorori.String{h.menu_inv_maxvectorori.Value}
            case 'RMS'
                %                 eloreta.P.img = abs(squeeze(nanmean(rms(swf(h.cfg.study.bl_bmf.act_samps,:,:),2))));
                eloreta.P.img = abs(squeeze(nanmean(swf_pwr)))';
                % Approximating orientations using RMS swf
                x=swf_pwr(1,:)./max(swf_pwr,[],1);
                y=swf_pwr(2,:)./max(swf_pwr,[],1);
                z=swf_pwr(3,:)./max(swf_pwr,[],1);
                eloreta.ori=[x;y;z]';
            case 'Max'    % vector dipole orientation with maximum swf power between active and control interval
                [img,max_ori]=max(swf_pwr); eloreta.P.img = img';
                ori = zeros(size(max_ori,2),3); for v=1:size(ori,1); ori(v,max_ori(v))=1; end
                eloreta.ori = ori;
            case 'avg.pow'
                eloreta.P.img = squeeze(eloreta.avg.pow(inside_idx));  % source image based on recalcuated swf data
                ori = cell2mat(eloreta.avg.ori(inside_idx));
                eloreta.ori = ori';
        end

        nd=sort(eloreta.P.img); eloreta.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;

        pow = eloreta.avg.pow;
        eloreta = rmfield(eloreta,{'avg'}); %MNE.avg = rmfield(MNE.avg,{'mom','filter','noisecov'});
        eloreta.avg.pow = single(pow);

        h.inv_soln(h.next_inv_soln).Type = inv_soln;
        h.inv_soln(h.next_inv_soln).soln = eloreta;
        h.inv_soln(h.next_inv_soln).maxvectorori_Type = h.menu_inv_maxvectorori.String{h.menu_inv_maxvectorori.Value};
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type,h.menu_inv_maxvectorori.String{h.menu_inv_maxvectorori.Value}, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);

        %% adding covariance in BRANELab format
        h.inv_soln(h.next_inv_soln).soln.RNcov.R = avg_data.cov;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Rbar = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.N = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.Nbar = [];

        %                 %% Plotting eLORETA on new figure
        %                     img=eloreta.P.img; ori=eloreta.ori; null_thresh=max(img)*.15;
        %                 min_max=[min(img) max(img)*.85]; %min_max = [0 5];
        %                 sFaceAlpha=0; vol_types=1;
        %                 figure(1006); clf; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT_new(img,ori,null_thresh,15,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,jet(255),min_max,h.anatomy.vol.bnd(3),[],...
        %                     h.cfg.study.bl_bmf.vw_angle,1,vol_types,h.inv_soln(h.current_inv_soln).leadfield.pos,h.cfg.study.bl_bmf.inside_idx,h.cfg.study.bl_bmf.slice_orient,h.cfg.study.bl_bmf.sFaceAlpha);
        %                 seed_idx = 1:3;  ln_wdth = 1; ln_wdth2 = 2;
        %                 hold on; mrk_size=150; s1=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'k+','sizedata',mrk_size,'linewidth',ln_wdth); s2=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'ks','sizedata',mrk_size,'linewidth',ln_wdth2);
        %
        %                 title('FT eLORETA'); colorbar; caxis(map_scale);
    case 'sLORETA (FT)'
        act_int = str2num(h.edit_act_int.String);
        ctrl_int = str2num(h.edit_ctrl_int.String);
        inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);

        %% averaging data and getting covariance matrix
        cfg = []; cfg.baseline = h.sim_data.cfg.study.base_int;
        cfg.channel = h.anatomy.sens.label(h.anatomy.sens.good_sensors);
        ft_data = ft_timelockbaseline(cfg, ft_data);
        cfg = [];
        cfg.channel = h.anatomy.sens.label(h.anatomy.sens.good_sensors);
        cfg.covariance='yes';
        cfg.covariancewindow = [ctrl_int(1) act_int(2)]; % calculating covariance across entire ctrl and act interval just as stated in Field Trip's Tutorial
        avg_data = ft_timelockanalysis(cfg,ft_data);

        %% FT sLORETA
        cfg=[];
        cfg.method = 'sloreta';
        cfg.sloreta = ft_cfg; % using same additional options for all inverse solns
        cfg.keepleadfield = cfg.sloreta.keepleadfield;
        cfg.reducerank = cfg.sloreta.reducerank;
        cfg.normalize = cfg.sloreta.normalize;
        cfg.normalizeparam = cfg.sloreta.normalizeparam;
        cfg.sloreta=rmfield(cfg.sloreta,'keepleadfield');  cfg.sloreta=rmfield(cfg.sloreta,'reducerank');
        cfg.sloreta=rmfield(cfg.sloreta,'normalize');  cfg.sloreta=rmfield(cfg.sloreta,'normalizeparam');

        try
            sourcemodel = h.inv_soln(h.next_inv_soln).leadfield;
            sourcemodel.leadfield = sourcemodel.leadfield';  % transpose for sLORETA - bug in Field Trips program. sLORETA is the only invSoln that needs to do this
            cfg.grad                = h.inv_soln(h.next_inv_soln).sens;
            cfg.sourcemodel         = sourcemodel;
            cfg.headmodel           = h.inv_soln(h.next_inv_soln).headmodel;
            sloreta  = ft_sourceanalysis(cfg,avg_data);
            sloreta = ft_sourcedescriptives([], sloreta); % This is gets the nai to get sloreta.P.img
        catch
            sourcemodel = h.inv_soln(h.next_inv_soln).leadfield;
            sourcemodel.leadfield = sourcemodel.leadfield;  % transpose for sLORETA - bug in Field Trips program. sLORETA is the only invSoln that needs to do this
            cfg.grad                = h.inv_soln(h.next_inv_soln).sens;
            cfg.sourcemodel         = sourcemodel;
            cfg.headmodel           = h.inv_soln(h.next_inv_soln).headmodel;
            sloreta  = ft_sourceanalysis(cfg,avg_data);
            sloreta = ft_sourcedescriptives([], sloreta); % This is gets the nai to get sloreta.P.img
        end

        %% calculating image power from filter weights;
        base_samps=h.sim_data.cfg.study.base_samps;     % -300 to 0 ms
        act_samps=h.cfg.study.bl_bmf.act_samps;
        ctrl_samps=h.cfg.study.bl_bmf.ctrl_samps;
        ctrl_samps1 = ctrl_samps(1:round(length(ctrl_samps)/2));
        ctrl_samps2 = ctrl_samps(round(length(ctrl_samps)/2)+1:end);

        y=cell2mat(sloreta.avg.filter); %sloreta.wts=permute(reshape(y,[3 size(y,2)/length(inside_idx) length(inside_idx)]),[3 1 2]);
        sloreta.wts = y';
        %         sloreta.wts = permute(sloreta.wts,[3 1 2]);

        %                 sloreta.P.img=squeeze(sloreta.avg.pow(inside_idx));  % source image based on recalcuated swf data
        sloreta.P.img=sloreta.avg.nai; % nai exactly same as -->  squeeze(sloreta.avg.pow(inside_idx)./sloreta.avg.noise(inside_idx));  % source image based on recalculated swf data

        nd=sort(sloreta.P.img); sloreta.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;

        ori = cell2mat(sloreta.avg.ori(inside_idx));
        sloreta.ori = ori';

        sloreta = rmfield(sloreta,{'avg'}); %MNE.avg = rmfield(MNE.avg,{'mom','filter','noisecov'});

        h.inv_soln(h.next_inv_soln).Type = 'sLORETA (FT)';
        h.inv_soln(h.next_inv_soln).soln = sloreta;
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
        %% adding covariance in BRANELab format
        h.inv_soln(h.next_inv_soln).soln.RNcov.R = avg_data.cov;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Rbar = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.N = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.Nbar = [];

        %         %% PLotting sloreta on new figure
        %             img=sloreta.P.img; ori=sloreta.ori; null_thresh=sloreta.null_thresh*2;
        %         min_max=[min(img) max(img)*.15]; %min_max = [0 5];
        %         sFaceAlpha=0; vol_types=1;
        %         figure(1006); clf; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT_new(img,ori,null_thresh,15,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,jet(255),min_max,h.anatomy.vol.bnd(3),[],...
        %             h.cfg.study.bl_bmf.vw_angle,1,vol_types,h.inv_soln(h.current_inv_soln).leadfield.pos,h.cfg.study.bl_bmf.inside_idx,h.cfg.study.bl_bmf.slice_orient,h.cfg.study.bl_bmf.sFaceAlpha);
        %         seed_idx = 1:3;  ln_wdth = 1; ln_wdth2 = 2;
        %         hold on; mrk_size=150; s1=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'k+','sizedata',mrk_size,'linewidth',ln_wdth); s2=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'ks','sizedata',mrk_size,'linewidth',ln_wdth2);
        %
        %         title('FT sLORETA'); colorbar; caxis(map_scale);
    case 'MNE (FT)'
        %         act_int = str2num(h.edit_act_int.String);
        %         ctrl_int = str2num(h.edit_ctrl_int.String);
        %         inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);

        %% averaging data and getting covariance matrix
        cfg = []; cfg.baseline = h.sim_data.cfg.study.base_int;
        cfg.channel = h.anatomy.sens.label(h.anatomy.sens.good_sensors);
        ft_data = ft_timelockbaseline(cfg, ft_data);
        cfg = [];
        cfg.channel = h.anatomy.sens.label(h.anatomy.sens.good_sensors);
        cfg.covariance='yes';
        %         cfg.covariancewindow = [ctrl_int(1) act_int(2)]; % calculating covariance across entire ctrl and act interval just as stated in Field Trip's Tutorial
        cfg.covariancewindow = [-inf 0];
        avg_data = ft_timelockanalysis(cfg,ft_data);

        %% FT MNE
        cfg = [];
        cfg.mne = ft_cfg;
        cfg.keepleadfield = cfg.mne.keepleadfield;
        cfg.reducerank = cfg.mne.reducerank;
        cfg.normalize = cfg.mne.normalize;
        cfg.normalizeparam = cfg.mne.normalizeparam;
        cfg.mne=rmfield(cfg.mne,'keepleadfield');  cfg.mne=rmfield(cfg.mne,'reducerank');
        cfg.mne=rmfield(cfg.mne,'normalize');  cfg.mne=rmfield(cfg.mne,'normalizeparam');

        cfg.method              = 'mne';
        cfg.mne.fixedori        = 'yes'; %  'yes' = get single orientation using svd
        cfg.grad                = h.inv_soln(h.next_inv_soln).sens;
        cfg.sourcemodel         = h.inv_soln(h.next_inv_soln).leadfield;
        cfg.headmodel           = h.inv_soln(h.next_inv_soln).headmodel;
        cfg.mne.scalesourcecov = h.menu_inv_scalesourcecov.String{h.menu_inv_scalesourcecov.Value};

        MNE = ft_sourceanalysis(cfg,avg_data);
        MNE = ft_sourcedescriptives([], MNE); % This is required to baseline MNE data to get it to fit better for deeper sources

        % calculating image power from filter weights;
        y=cell2mat(permute(MNE.avg.filter,[2 1])); MNE.wts=permute(reshape(y,[3 size(y,1)/3 size(y,2)]),[2 1 3]);
        swf=[];swf_snr=[];swf_pwr=[];MNE.noise_pwr=[];
        for ox=1:3
            swf(:,ox,:)=avg_data.avg'*squeeze(MNE.wts(:,ox,:))';
            %    swf_snr(:,ox,:)=20*log10(bsxfun(@rdivide,swf(:,ox,:),nanmean(swf(h.cfg.study.base_samps,ox,:),1)));
            %   swf_snr(:,ox,:)=(bsxfun(@rdivide,swf(:,ox,:),nanmean(swf(h.cfg.study.base_samps,ox,:),1))); % SNR in percent from baseline
            %            swf_snr(:,ox,:)=(bsxfun(@rdivide,swf(:,ox,:),nanstd(swf(h.cfg.study.base_samps,ox,:),1))); % normalized to baseline
            swf_pwr(ox,:)=rms(swf(h.cfg.study.bl_bmf.act_samps,ox,:),1)./rms(swf(h.cfg.study.bl_bmf.ctrl_samps,ox,:),1);
            %            MNE.noise_pwr(ox,:)=rms(swf(h.cfg.study.bl_bmf.ctrl_samps,ox,:),1)./rms(swf(h.cfg.study.bl_bmf.ctrl_samps,ox,:),1);
        end
        MNE.wts = permute(MNE.wts,[3 1 2]);
        switch h.menu_inv_maxvectorori.String{h.menu_inv_maxvectorori.Value}
            case 'RMS'
                %                 MNE.P.img = abs(squeeze(nanmean(rms(swf(h.cfg.study.bl_bmf.act_samps,:,:),2))));
                MNE.P.img = abs(squeeze(nanmean(swf_pwr)))';
                % Approximating orientations using RMS swf
                x=swf_pwr(1,:)./max(swf_pwr,[],1);
                y=swf_pwr(2,:)./max(swf_pwr,[],1);
                z=swf_pwr(3,:)./max(swf_pwr,[],1);
                MNE.ori=[x;y;z]';
            case 'Max'    % vector dipole orientation with maximum swf power between active and control interval
                [img,max_ori]=max(swf_pwr); MNE.P.img = img';
                ori = zeros(size(max_ori,2),3); for v=1:size(ori,1); ori(v,max_ori(v))=1; end
                MNE.ori = ori;
                %        MNE.P.img=squeeze(max(abs(MNE.avg.pow(h.cfg.study.bl_bmf.inside_idx,h.cfg.study.bl_bmf.act_samps))'))';
                %    MNE.P.img=squeeze(max(max(swf_snr(h.cfg.study.bl_bmf.act_samps,:,:))));
                %    MNE.P.img=squeeze(rms(swf_pwr,1))';   % averaging across orientations
                %    MNE.P.img=squeeze(rms(swf_pwr./MNE.noise_pwr,1))';   % averaging across orientations

            case 'avg.pow'
                MNE.P.img=squeeze(max(abs(MNE.avg.pow(h.cfg.study.bl_bmf.inside_idx,h.cfg.study.bl_bmf.act_samps))'))';
                [~,max_ori]=max(swf_pwr);
                ori = zeros(size(max_ori,2),3); for v=1:size(ori,1); ori(v,max_ori(v))=1; end
                MNE.ori = ori;
        end
        nd=sort(MNE.P.img); MNE.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;
        %        MNE.P.img = MNE.avg.pow(inside_idx,s);


        %% removing swf and swf_pwr to save memory storage
        %         MNE = rmfield(MNE,{'swf','swf_pwr'});
        MNE.avg = rmfield(MNE.avg,{'mom','filter','noisecov'});
        MNE.avg.pow = single(MNE.avg.pow);
        h.inv_soln(h.next_inv_soln).Type = inv_soln;
        h.inv_soln(h.next_inv_soln).soln = MNE;
        %% adding covariance in BRANELab format
        h.inv_soln(h.next_inv_soln).soln.RNcov.R = avg_data.cov;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Rbar = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.N = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.Nbar = [];

        h.inv_soln(h.next_inv_soln).maxvectorori_Type = h.menu_inv_maxvectorori.String{h.menu_inv_maxvectorori.Value};
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.menu_inv_maxvectorori.String{h.menu_inv_maxvectorori.Value}, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
        %% Plotting MNE
        %         img=MNE.P.img; ori=MNE.ori; null_thresh=max(img)*.25; %MNE.null_thresh*2;
        %
        %         min_max=[min(img) max(img)*.15]; %min_max = [0 5];
        %         sFaceAlpha=0; vol_types=1;
        %         figure(1006); clf; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT_new(img,ori,null_thresh,15,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,jet(255),min_max,h.anatomy.vol.bnd(3),[],...
        %             h.cfg.study.bl_bmf.vw_angle,1,vol_types,h.inv_soln(h.current_inv_soln).leadfield.pos,h.cfg.study.bl_bmf.inside_idx,h.cfg.study.bl_bmf.slice_orient,h.cfg.study.bl_bmf.sFaceAlpha);
        %         seed_idx = 1:3;  ln_wdth = 1; ln_wdth2 = 2;
        %         hold on; mrk_size=150; s1=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'k+','sizedata',mrk_size,'linewidth',ln_wdth); s2=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'ks','sizedata',mrk_size,'linewidth',ln_wdth2);
        %        title('FT MNE'); colorbar; caxis(min_max);
        %
    case 'SAM (FT)'
        act_int = str2num(h.edit_act_int.String);
        ctrl_int = str2num(h.edit_ctrl_int.String);
        inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);

        cfg.method              = 'sam';
        cfg.grad                = h.inv_soln(h.next_inv_soln).sens;
        cfg.sourcemodel         = h.inv_soln(h.next_inv_soln).leadfield;
        cfg.headmodel           = h.inv_soln(h.next_inv_soln).headmodel;

        %% calculating directly - using BRANELab's covariance
        xs=round(size(h.cfg.study.bl_bmf.ctrl_samps,2)/2);
        [R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(data, h.cfg.study.bl_bmf.act_samps, h.cfg.study.bl_bmf.ctrl_samps(xs+1:end));
        [avg] = ft_inverse_sam(cfg.sourcemodel, cfg.grad, cfg.headmodel, nanmean(permute(data,[2 1 3]),3), R);
        SAM.avg = avg;
        h.SAM = SAM;

        %% adding covariance in BRANELab format
        h.inv_soln(h.next_inv_soln).soln.RNcov.R = R;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Rbar = Rbar;
        h.inv_soln(h.next_inv_soln).soln.RNcov.N = N;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Nbar = Nbar;

        %% creating SAM power image and estimating orientations
        SAM.P.img = SAM.avg.pseudoZ'; %
        SAM.P.img = SAM.avg.pow'; %

        nd=sort(SAM.P.img); SAM.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;
        ori=cell2mat(SAM.avg.ori(inside_idx)); SAM.ori=ori';
        wts = cell2mat(SAM.avg.filter(inside_idx));
        wts = reshape(wts, size(SAM.avg.filter{1},2), size(SAM.avg.filter,2));
        SAM.wts = wts;

        %                 img=SAM.P.img; ori=SAM.ori; null_thresh=max(img)*.55;
        %                 min_max=[min(img) max(img)*.95];  vol_types=1;
        %                 seed_idx = 1:3;  ln_wdth = 1; ln_wdth2 = 2;
        %                 figure(1004); clf; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT_new(img,ori,null_thresh,8,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,jet(255),min_max,h.anatomy.vol.bnd(3),[],...
        %                     h.cfg.study.bl_bmf.vw_angle,1,vol_types,h.inv_soln(h.current_inv_soln).leadfield.pos,inside_idx,h.cfg.study.bl_bmf.slice_orient,h.cfg.study.bl_bmf.sFaceAlpha);
        %                 hold on; mrk_size=150; s1=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'k+','sizedata',mrk_size,'linewidth',ln_wdth); s2=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'ks','sizedata',mrk_size,'linewidth',ln_wdth2);
        %                 title('FT SAM'); colorbar; caxis(min_max);

        % removing 'avg' to reduce storage space
        SAM = rmfield(SAM,{'avg'}); %MNE.avg = rmfield(MNE.avg,{'mom','filter','noisecov'});

        h.inv_soln(h.next_inv_soln).Type = inv_soln;
        h.inv_soln(h.next_inv_soln).soln = SAM;
        %         h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, cfg.weightnorm ,h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s pseudoZ (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case {'dics (FT)'}   % requires Time-Freq do have ben computed
        act_int = str2num(h.edit_act_int.String);
        ctrl_int = str2num(h.edit_ctrl_int.String);
        inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
        if isempty(h.fieldtrip_tfr_data)
            h.inv_soln(h.next_inv_soln).Type = 'dics';
            h.inv_soln(h.next_inv_soln).soln = [];
            h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s Not Computed',h.inv_soln(h.next_inv_soln).Type);
        else
            %% FT Source modeling configuration parameters
            cfg=[];
            cfg.method = inv_soln;
            cfg.dics = ft_cfg; % using same additional options for all inverse solns
            cfg.keepleadfield = cfg.dics.keepleadfield;
            cfg.reducerank = cfg.dics.reducerank;
            cfg.normalize = cfg.dics.normalize;
            cfg.normalizeparam = cfg.dics.normalizeparam;
            cfg.dics=rmfield(cfg.dics,'keepleadfield');  cfg.dics=rmfield(cfg.dics,'reducerank');
            cfg.dics=rmfield(cfg.dics,'normalize');  cfg.dics=rmfield(cfg.dics,'normalizeparam');

            cfg.grad                = h.inv_soln(h.next_inv_soln).sens;
            cfg.sourcemodel         = h.inv_soln(h.next_inv_soln).leadfield;
            cfg.headmodel           = h.inv_soln(h.next_inv_soln).headmodel;
            %% running source modeling
            dics = ft_sourceanalysis(cfg, h.fieldtrip_tfr_data);
            dics = ft_sourcedescriptives([], dics); % to get the neural-activity-index
            dics.P.img = dics.avg.nai;

            %% dics power image and orientations

            nd=sort(dics.P.img); dics.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;
            ori=cell2mat(dics.avg.ori(inside_idx)); dics.ori=ori';
            wts = cell2mat(dics.avg.filter(inside_idx)); dics.wts = wts';

            %         img=dics.P.img; ori=dics.ori; null_thresh=max(img)*.05;
            %         min_max=[min(img) max(img)*.95];  vol_types=1;
            %         seed_idx = 1:3;  ln_wdth = 1; ln_wdth2 = 2;
            %         figure(1004); clf; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT_new(img,ori,null_thresh,8,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,jet(255),min_max,h.anatomy.mesh_volumes(3),[],...
            %             h.cfg.study.bl_bmf.vw_angle,1,vol_types,h.inv_soln(h.current_inv_soln).leadfield.pos,inside_idx,h.cfg.study.bl_bmf.slice_orient,h.cfg.study.bl_bmf.sFaceAlpha);
            %         hold on; mrk_size=150; s1=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'k+','sizedata',mrk_size,'linewidth',ln_wdth); s2=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'ks','sizedata',mrk_size,'linewidth',ln_wdth2);
            %         title('FT dics'); colorbar; caxis(min_max);

            % removing 'avg' to reduce storage space
            dics = rmfield(dics,{'avg'}); %MNE.avg = rmfield(MNE.avg,{'mom','filter','noisecov'});

            h.inv_soln(h.next_inv_soln).Type = 'dics';
            h.inv_soln(h.next_inv_soln).soln = dics;
            h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
        end
        %% adding covariance in BRANELab format
        h.inv_soln(h.next_inv_soln).soln.RNcov.R = avg_data.cov;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Rbar = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.N = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.Nbar = [];
        %                 %% Plotting dics on new figure
        %                     img=dics.P.img; ori=dics.ori; null_thresh=max(img)*.15;
        %                 min_max=[min(img) max(img)*.85]; %min_max = [0 5];
        %                 sFaceAlpha=0; vol_types=1;
        %                 figure(1006); clf; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT_new(img,ori,null_thresh,15,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,jet(255),min_max,h.anatomy.vol.bnd(3),[],...
        %                     h.cfg.study.bl_bmf.vw_angle,1,vol_types,h.inv_soln(h.current_inv_soln).leadfield.pos,h.cfg.study.bl_bmf.inside_idx,h.cfg.study.bl_bmf.slice_orient,h.cfg.study.bl_bmf.sFaceAlpha);
        %                 seed_idx = 1:3;  ln_wdth = 1; ln_wdth2 = 2;
        %                 hold on; mrk_size=150; s1=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'k+','sizedata',mrk_size,'linewidth',ln_wdth); s2=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'ks','sizedata',mrk_size,'linewidth',ln_wdth2);
        %
        %                 title('FT dics'); colorbar; caxis(map_scale);
    case {'pcc (FT)'}
        act_int = str2num(h.edit_act_int.String);
        ctrl_int = str2num(h.edit_ctrl_int.String);
        inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);

        %% averaging data and getting covariance matrix
        if h.menu_inv_datatype.Value == 1     % use averaged data
            cfg = []; cfg.baseline = h.sim_data.cfg.study.base_int;
            cfg.channel = h.anatomy.sens.label(h.anatomy.sens.good_sensors);
            ft_data = ft_timelockbaseline(cfg, ft_data);
            cfg = [];
            cfg.channel = h.anatomy.sens.label(h.anatomy.sens.good_sensors);
            cfg.covariance='yes';
            cfg.covariancewindow = [ctrl_int(1) act_int(2)]; % calculating covariance across entire ctrl and act interval just as stated in Field Trip's Tutorial
            avg_data = ft_timelockanalysis(cfg,ft_data);
        end

        %% FT Source modeling configuration parameters
        cfg=[];
        cfg.method = inv_soln;
        cfg.pcc = ft_cfg; % using same additional options for all inverse solns
        cfg.keepleadfield = cfg.pcc.keepleadfield;
        cfg.reducerank = cfg.pcc.reducerank;
        cfg.normalize = cfg.pcc.normalize;
        cfg.normalizeparam = cfg.pcc.normalizeparam;
        cfg.pcc=rmfield(cfg.pcc,'keepleadfield');  cfg.pcc=rmfield(cfg.pcc,'reducerank');
        cfg.pcc=rmfield(cfg.pcc,'normalize');  cfg.pcc=rmfield(cfg.pcc,'normalizeparam');

        cfg.grad                = h.inv_soln(h.next_inv_soln).sens;
        cfg.sourcemodel         = h.inv_soln(h.next_inv_soln).leadfield;
        cfg.headmodel           = h.inv_soln(h.next_inv_soln).headmodel;

        %% running source modeling
        if h.menu_inv_datatype.Value == 1     % use averaged data
            pcc  = ft_sourceanalysis(cfg,avg_data);
            pcc = ft_sourcedescriptives([], pcc); % to get the neural-activity-index
            pcc.P.img = pcc.avg.nai; %pcc.avg.pow(inside_idx)./pcc.avg.noise(inside_idx);   % simple ratio
        elseif h.menu_inv_datatype.Value == 2  % conduct source modeling on time-Freq data
            pcc = ft_sourceanalysis(cfg, h.fieldtrip_tfr_data);
            pcc = ft_sourcedescriptives([], pcc); % to get the neural-activity-index
            pcc.P.img = pcc.avg.nai;
        end

        %% pcc power image and orientations

        nd=sort(pcc.P.img); pcc.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;
        ori=cell2mat(pcc.avg.ori(inside_idx)); pcc.ori=ori';
        wts = cell2mat(pcc.avg.filter(inside_idx)); pcc.wts = wts';

        %         img=pcc.P.img; ori=pcc.ori; null_thresh=max(img)*.05;
        %         min_max=[min(img) max(img)*.95];  vol_types=1;
        %         seed_idx = 1:3;  ln_wdth = 1; ln_wdth2 = 2;
        %         figure(1004); clf; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT_new(img,ori,null_thresh,8,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,jet(255),min_max,h.anatomy.mesh_volumes(3),[],...
        %             h.cfg.study.bl_bmf.vw_angle,1,vol_types,h.inv_soln(h.current_inv_soln).leadfield.pos,inside_idx,h.cfg.study.bl_bmf.slice_orient,h.cfg.study.bl_bmf.sFaceAlpha);
        %         hold on; mrk_size=150; s1=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'k+','sizedata',mrk_size,'linewidth',ln_wdth); s2=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'ks','sizedata',mrk_size,'linewidth',ln_wdth2);
        %         title('FT pcc'); colorbar; caxis(min_max);

        % removing 'avg' to reduce storage space
        pcc = rmfield(pcc,{'avg'}); %MNE.avg = rmfield(MNE.avg,{'mom','filter','noisecov'});

        h.inv_soln(h.next_inv_soln).Type = 'pcc';
        h.inv_soln(h.next_inv_soln).soln = pcc;
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);

        %% adding covariance in BRANELab format
        h.inv_soln(h.next_inv_soln).soln.RNcov.R = avg_data.cov;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Rbar = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.N = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.Nbar = [];
        %                 %% Plotting pcc on new figure
        %                     img=pcc.P.img; ori=pcc.ori; null_thresh=max(img)*.15;
        %                 min_max=[min(img) max(img)*.85]; %min_max = [0 5];
        %                 sFaceAlpha=0; vol_types=1;
        %                 figure(1006); clf; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT_new(img,ori,null_thresh,15,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,jet(255),min_max,h.anatomy.vol.bnd(3),[],...
        %                     h.cfg.study.bl_bmf.vw_angle,1,vol_types,h.inv_soln(h.current_inv_soln).leadfield.pos,h.cfg.study.bl_bmf.inside_idx,h.cfg.study.bl_bmf.slice_orient,h.cfg.study.bl_bmf.sFaceAlpha);
        %                 seed_idx = 1:3;  ln_wdth = 1; ln_wdth2 = 2;
        %                 hold on; mrk_size=150; s1=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'k+','sizedata',mrk_size,'linewidth',ln_wdth); s2=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'ks','sizedata',mrk_size,'linewidth',ln_wdth2);
        %
        %                 title('FT pcc'); colorbar; caxis(map_scale);
    case 'sMCMV'        % Alex Moiseev's sub-space MCMV

        % subspace MCMV
        sInSMCMV.beamType = h.menu_SPA_loc_flag.String{h.menu_SPA_loc_flag.Value};
        sInSMCMV.R = R;
        sInSMCMV.arrN = nN; %N;
        sInSMCMV.arrH = arrH;
        sInSMCMV.lstFlag = h.anatomy.moiseev.lstFlag;
        sInSMCMV.dims = h.anatomy.moiseev.dims;
        sInSMCMV.nSrc = maxSrc;
        sInSMCMV.Cavg = Rbar;        % This is C-bar - only needed for MER
        sInSMCMV.bMCMVS = true;    % Should be ALWAYS set to true
        sInSMCMV.pVal = 1;        % p-value (kind of) for the peak to be considered real. Set to 1 to get all peaks
        sInSMCMV.gap = gap;
        sInSMCMV.bVerbose = true;
        sInSMCMV.bPlotLambda = false;
        sInSMCMV.bRAPBeam = false;   % Choose RAP (true) or SMCMV (false). You can try both
        sInSMCMV.bDoTRAP = false;             % Don't ask - just set to false :)
        sOutSMCMV = doSMCMV(sInSMCMV);

        %% just finding max peaks in each images run -- assuming that each image in fImg corresponds to a single peak
        mcmv_voxel=[]; img_max=[]; pidx=[];
        for v=1:maxSrc
            [img_max(v),pidx(v)] = max(sOutSMCMV.fImg(v,h.anatomy.moiseev.inside_idx)');
            mcmv_voxel(v,:) = h.anatomy.moiseev.vx_grid(h.anatomy.moiseev.inside_idx(pidx(v)),:);
        end

        nd=sort(reshape(sOutSMCMV.fImg,numel(sOutSMCMV.fImg),1)); nd = nd(~isnan(nd));
        h.inv_soln(h.next_inv_soln).soln.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;

        % finding voxel locations within original anatomical leadfield grid
        m_idx = find_nearest_voxel(mcmv_voxel,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.inv_soln(h.current_inv_soln).leadfield.inside==1,:));
        %         true_idx = h.cfg.source.vx_idx;

        %         figure(998); clf; set(gcf,'color','w'); hold on;
        %         opt.vol_nums=1; vol = h.anatomy.mesh_volumes(3);
        %         [p1]=bl_plot_mesh(vol,opt); view(-90,90);
        %         scatter3(vx_pos(true_idx,1),vx_pos(true_idx,2),vx_pos(true_idx,3),'ks','SizeData',220,'linewidth',2);
        %         scatter3(vx_pos(m_idx,1),vx_pos(m_idx,2),vx_pos(m_idx,3),'r.','SizeData',220,'linewidth',2);

        % weights wts=signal    wts_noise = noise
        h.inv_soln(h.next_inv_soln).soln.wts = zeros(size(H,3),size(H,1));
        h.inv_soln(h.next_inv_soln).soln.wts(:,m_idx) = sOutSMCMV.Wr;
        %         h.inv_soln(h.next_inv_soln).soln.wts_noise = zeros(size(H,3),size(H,1));
        h.inv_soln(h.next_inv_soln).soln.wts_noise = sOutSMCMV.Wr;

        % P.img
        h.inv_soln(h.next_inv_soln).soln.P.img = zeros(size(H,1),1);
        h.inv_soln(h.next_inv_soln).soln.P.img(m_idx) = real(img_max);
        h.inv_soln(h.next_inv_soln).soln.type = sInSMCMV.beamType;
        % MCMV indices & Orientations
        h.inv_soln(h.next_inv_soln).soln.MCMV_idx = m_idx;
        h.inv_soln(h.next_inv_soln).soln.ori = zeros(size(H,1),3);
        h.inv_soln(h.next_inv_soln).soln.ori(m_idx,:) = sOutSMCMV.U3D;
        % add in Covariance data
        h.inv_soln(h.next_inv_soln).soln.RNcov.R = R;
        h.inv_soln(h.next_inv_soln).soln.RNcov.N = nN;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Rbar = Rbar;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Nbar = [];

        sOutSMCMV = rmfield(sOutSMCMV,'arrH');  % to reduce storage
        h.inv_soln(h.next_inv_soln).Type = inv_soln;
        h.inv_soln(h.next_inv_soln).soln.sOutSMCMV = sOutSMCMV;
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case 'bRAPBeam' % Alex Moiseev's bRAP MCMV beamformer

        % subspace MCMV
        sInSMCMV.beamType = h.menu_SPA_loc_flag.String{h.menu_SPA_loc_flag.Value};
        sInSMCMV.R = R;
        sInSMCMV.arrN = N;
        sInSMCMV.arrH = arrH;
        sInSMCMV.lstFlag = h.anatomy.moiseev.lstFlag;
        sInSMCMV.dims = h.anatomy.moiseev.dims;
        sInSMCMV.nSrc = maxSrc;
        sInSMCMV.Cavg = Rbar;        % This is C-bar - only needed for MER
        sInSMCMV.bMCMVS = true;    % Should be ALWAYS set to true
        sInSMCMV.pVal = 1;        % p-value (kind of) for the peak to be considered real. Set to 1 to get all peaks
        sInSMCMV.gap = gap;
        sInSMCMV.bVerbose = true;
        sInSMCMV.bPlotLambda = false;
        sInSMCMV.bRAPBeam = true;   % Choose RAP (true) or SMCMV (false). You can try both
        sInSMCMV.bDoTRAP = false;             % Don't ask - just set to false :)
        sOutSMCMV = doSMCMV(sInSMCMV);

        %% just finding max peaks in each images run -- assuming that each image in fImg corresponds to a single peak
        mcmv_voxel=[]; img_max=[]; pidx=[];
        for v=1:maxSrc
            [img_max(v),pidx(v)] = max(sOutSMCMV.fImg(v,h.anatomy.moiseev.inside_idx)');
            mcmv_voxel(v,:) = h.anatomy.moiseev.vx_grid(h.anatomy.moiseev.inside_idx(pidx(v)),:);
        end

        nd=sort(reshape(sOutSMCMV.fImg,numel(sOutSMCMV.fImg),1)); nd = nd(~isnan(nd));
        h.inv_soln(h.next_inv_soln).soln.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;


        % finding voxel locations within original anatomical leadfield grid
        m_idx = find_nearest_voxel(mcmv_voxel,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.inv_soln(h.current_inv_soln).leadfield.inside==1,:));

        % weights wts=signal    wts_noise = noise
        h.inv_soln(h.next_inv_soln).soln.wts = zeros(size(H,3),size(H,1));
        h.inv_soln(h.next_inv_soln).soln.wts(:,m_idx) = sOutSMCMV.Wr;
        %         h.inv_soln(h.next_inv_soln).soln.wts_noise = zeros(size(H,3),size(H,1));
        h.inv_soln(h.next_inv_soln).soln.wts_noise = sOutSMCMV.Wr;

        % P.img
        h.inv_soln(h.next_inv_soln).soln.P.img = zeros(size(H,1),1);
        h.inv_soln(h.next_inv_soln).soln.P.img(m_idx) = real(img_max);
        h.inv_soln(h.next_inv_soln).soln.type = sInSMCMV.beamType;
        % MCMV indices & Orientations
        h.inv_soln(h.next_inv_soln).soln.MCMV_idx = m_idx;
        h.inv_soln(h.next_inv_soln).soln.ori = zeros(size(H,1),3);
        h.inv_soln(h.next_inv_soln).soln.ori(m_idx,:) = sOutSMCMV.U3D;
        % add in Covariance data
        h.inv_soln(h.next_inv_soln).soln.RNcov.R = R;
        h.inv_soln(h.next_inv_soln).soln.RNcov.N = N;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Rbar = Rbar;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Nbar = [];

        sOutSMCMV = rmfield(sOutSMCMV,'arrH');  % to reduce storage
        h.inv_soln(h.next_inv_soln).Type = inv_soln;
        h.inv_soln(h.next_inv_soln).soln.sOutSMCMV = sOutSMCMV;
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case 'TrapMUSIC'    % Alex Moiseev's code for TRAP MUSIC

        % TRAP MUSIC
        sInTRAP.R = R;            % Full covariance
        sInTRAP.arrN = N;    % Noise covariance
        sInTRAP.arrH = arrH;    % Lead fields
        sInTRAP.lstFlag = h.anatomy.moiseev.lstFlag; % Flags
        sInTRAP.dims = h.anatomy.moiseev.dims;    % ROI dimensions in voxels
        sInTRAP.nSrc = maxSrc;    % Max number of sources to extract
        sInTRAP.gap = gap;        % Min allowed number nodes between the peaks (default 2)
        sInTRAP.bPlotLambda = false;    % Set to false
        sInTRAP.bVerbose = true;    % Enables more printouts

        sOutTRAP = trapMUSIC(sInTRAP);


        %% just finding max peaks in each images run -- assuming that each image in fImg corresponds to a single peak
        mcmv_voxel=[]; img_max=[]; pidx=[];
        for v=1:maxSrc
            [img_max(v),pidx(v)] = max(sOutTRAP.fImg(v,h.anatomy.moiseev.inside_idx)');
            mcmv_voxel(v,:) = h.anatomy.moiseev.vx_grid(h.anatomy.moiseev.inside_idx(pidx(v)),:);
        end

        nd=sort(reshape(sOutTRAP.fImg,numel(sOutTRAP.fImg),1)); nd = nd(~isnan(nd));
        h.inv_soln(h.next_inv_soln).soln.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;

        % finding voxel locations within original anatomical leadfield grid
        m_idx = find_nearest_voxel(mcmv_voxel,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.inv_soln(h.current_inv_soln).leadfield.inside==1,:));

        % weights wts=signal    wts_noise = noise
        h.inv_soln(h.next_inv_soln).soln.wts = zeros(size(H,3),size(H,1));
        h.inv_soln(h.next_inv_soln).soln.wts(:,m_idx) = sOutTRAP.Wr;
        %         h.inv_soln(h.next_inv_soln).soln.wts_noise = zeros(size(H,3),size(H,1));
        h.inv_soln(h.next_inv_soln).soln.wts_noise = sOutTRAP.Wr;

        % P.img
        h.inv_soln(h.next_inv_soln).soln.P.img = zeros(size(H,1),1);
        h.inv_soln(h.next_inv_soln).soln.P.img(m_idx) = real(img_max);
        h.inv_soln(h.next_inv_soln).soln.type = inv_soln;
        % MCMV indices & Orientations
        h.inv_soln(h.next_inv_soln).soln.MCMV_idx = m_idx;
        h.inv_soln(h.next_inv_soln).soln.ori = zeros(size(H,1),3);
        h.inv_soln(h.next_inv_soln).soln.ori(m_idx,:) = sOutTRAP.U3D;
        % add in Covariance data
        h.inv_soln(h.next_inv_soln).soln.RNcov.R = R;
        h.inv_soln(h.next_inv_soln).soln.RNcov.N = N;
        h.inv_soln(h.next_inv_soln).soln.RNcov.Rbar = [];
        h.inv_soln(h.next_inv_soln).soln.RNcov.Nbar = [];

        sOutTRAP = rmfield(sOutTRAP,'arrH');  % to reduce storage
        h.inv_soln(h.next_inv_soln).Type = inv_soln;
        h.inv_soln(h.next_inv_soln).soln.sOutTRAP = sOutTRAP;

        %                 figure(998); clf; set(gcf,'color','w'); hold on;
        %         opt.vol_nums=1; vol = h.anatomy.mesh_volumes(3);
        %         true_idx = h.cfg.source.vx_idx
        %         [p1]=bl_plot_mesh(vol,opt); view(-90,90);
        %         scatter3(vx_pos(true_idx,1),vx_pos(true_idx,2),vx_pos(true_idx,3),'ks','SizeData',220,'linewidth',2);
        %         scatter3(vx_pos(m_idx,1),vx_pos(m_idx,2),vx_pos(m_idx,3),'r.','SizeData',220,'linewidth',2);
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
end

sm_calc_post_source_modeling();


