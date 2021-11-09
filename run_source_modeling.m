function run_source_modeling(varargin)
%% Run Source modeling
global h

% try
h.btn_3D_plot_peak_waves.Value=0;
h.next_inv_soln = length(h.inv_soln)+1;
h.current_inv_soln = h.next_inv_soln;
inv_soln = h.menu_inv_soln.String{h.menu_inv_soln.Value};
% fprintf('Inverse Modeling using %s\n',inv_soln);
% hm = msgbox(sprintf('Running Inverse Source Modeling\n\n          %s',inv_soln));
if h.monte_carlo_flag == 1
    h.waitfor_txt.String = sprintf('Inverse Modeling using %s\n',inv_soln); drawnow;
else
    h.waitfor_panel.Visible='on';
    h.waitfor_txt.String = sprintf('Inverse Modeling using %s\n',inv_soln); drawnow;
end

%% make-shift paramaters right now to test simulations
if h.menu_inv_headmodel.Value==1    % Whole Brain
    if h.menu_sens_type.Value == 1 %MEG
        h.inv_soln(h.next_inv_soln).leadfield =  h.anatomy.leadfield_meg_vol;
        h.inv_soln(h.next_inv_soln).headmodel =  h.anatomy.headmodel_meg_vol;
        h.inv_soln(h.next_inv_soln).sens =  h.anatomy.sens_meg;
        h.inv_soln(h.current_inv_soln).headmodel_type = 'Whole Brain';
%         h.anatomy.leadfield =  h.anatomy.leadfield_meg_vol;
%         h.anatomy.headmodel =  h.anatomy.headmodel_meg_vol;
%         h.anatomy.sens =  h.anatomy.sens_meg;
    elseif h.menu_sens_type.Value == 2 %EEG
        h.inv_soln(h.next_inv_soln).leadfield =  h.anatomy.leadfield_eeg_vol;
        h.inv_soln(h.next_inv_soln).headmodel =  h.anatomy.headmodel_eeg_vol;
        h.inv_soln(h.next_inv_soln).sens =  h.anatomy.sens_eeg;
        h.inv_soln(h.current_inv_soln).headmodel_type = 'Whole Brain';
%         h.anatomy.leadfield =  h.anatomy.leadfield_eeg_vol;
%         h.anatomy.headmodel =  h.anatomy.headmodel_eeg_vol;
%         h.anatomy.sens =  h.anatomy.sens_eeg;
    end
    h.inv_soln(h.next_inv_soln).headmodel_mesh =  h.anatomy.mesh_volumes(3);
    
elseif h.menu_inv_headmodel.Value==2    % Cortical Surface
    if h.menu_sens_type.Value == 1 %MEG
        h.inv_soln(h.next_inv_soln).leadfield =  h.anatomy.leadfield_meg_cortex;
        h.inv_soln(h.next_inv_soln).headmodel =  h.anatomy.headmodel_meg_cortex;
        h.inv_soln(h.next_inv_soln).sens =  h.anatomy.sens_meg;
        h.inv_soln(h.current_inv_soln).headmodel_type = 'Cortical Surface';
%         h.anatomy.leadfield =  h.anatomy.leadfield_meg_cortex;
%         h.anatomy.headmodel =  h.anatomy.headmodel_meg_cortex;
%         h.anatomy.sens =  h.anatomy.sens_meg;
    elseif h.menu_sens_type.Value == 2 %EEG
        
        h.inv_soln(h.next_inv_soln).leadfield =  h.anatomy.leadfield_eeg_cortex;
        h.inv_soln(h.next_inv_soln).headmodel =  h.anatomy.headmodel_eeg_cortex;
        h.inv_soln(h.next_inv_soln).sens =  h.anatomy.sens_eeg;
        h.inv_soln(h.current_inv_soln).headmodel_type = 'Cortical Surface';
%         h.anatomy.leadfield =  h.anatomy.leadfield_eeg_cortex;
%         h.anatomy.headmodel =  h.anatomy.headmodel_eeg_cortex;
%         h.anatomy.sens =  h.anatomy.sens_eeg;
    end
    h.inv_soln(h.next_inv_soln).headmodel_mesh =  h.anatomy.mesh_volumes(4);
    
end
h.current_3D_peak_idx =[]; h.current_3D_peak_voxels = []; h.current_norm_peak_swf =[]; h.current_peak_swf = [];
inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
h.cfg.study.bl_bmf.inside_idx = inside_idx;
mHref = []; mref = [];

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


%% checking rank of data. If it doesn't match the num_chans then white noise will be added until rank sufficient.
xs=round(size(h.cfg.study.bl_bmf.ctrl_samps,2)/2);
R=[]; t=0;
[R,N,Rbar,~,~,~]=BRANELab_calc_cov(data,h.cfg.study.bl_bmf.act_samps,h.cfg.study.bl_bmf.ctrl_samps(xs+1:end));
r2=rank(R);
xstd=min(min(std(data)))/size(data,2);
while r2<size(h.inv_soln(h.next_inv_soln).leadfield.H,1)
    t=t+1;
    data=data+(2*xstd*randn(size(data)));
    [R,N,Rbar,~,~,~]=BRANELab_calc_cov(data,h.cfg.study.bl_bmf.act_samps,h.cfg.study.bl_bmf.ctrl_samps(xs+1:end));
    r2=rank(R);
    fprintf('Rank = %.f\tChannel Count=%.f\n',r2,size(data,2));
    if 100*(2*t*xstd)/max(max(std(data)))>10; fprintf('WARNING! Insufficient Rank in Data or Data is empty\n'); end % stopping infinite loop
    % fprintf('Rank = %.f\n',r2);
    fprintf('Percent white noise added for sufficient rank = %.3f %%\n\n',100*(2*t*xstd)/max(max(std(data))));
end
h.cfg.study.bl_bmf.noise_alpha = str2num(h.edit_SPA_noise_alpha.String);


%% Preprocesing Data for Field Trip
switch h.menu_inv_soln.String{h.menu_inv_soln.Value}
    case {'SIA' 'MIA' 'SPA'}
        
        %     if h.menu_inv_soln.Value>3 && h.menu_inv_soln.Value<=7  % field Trip's Inv Solutions
    case {'LCMV' 'eLORETA' 'sLORETA' 'MNE' 'dics' 'pcc'}
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
tt0 = tic; % time to compute inv_soln
switch inv_soln
    case 'SPA'
        %% SPA
        
        h.cfg.study.bl_bmf.loc_flag = h.menu_SPA_loc_flag.Value-1;
        h.cfg.study.bl_bmf.noise_alpha = str2num(h.edit_SPA_noise_alpha.String);
        [SPA]=BRANElab_LCMV_beamformer_SPA(h.inv_soln(h.next_inv_soln).leadfield.H,data,...
            h.cfg.study.bl_bmf.act_samps,h.cfg.study.bl_bmf.ctrl_samps,h.cfg.study.bl_bmf.loc_flag,h.cfg.study.bl_bmf.noise_alpha);
        h.inv_soln(h.next_inv_soln).Type = 'SPA';
        h.inv_soln(h.next_inv_soln).soln = SPA;
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.menu_SPA_loc_flag.String{h.menu_SPA_loc_flag.Value} ,h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case 'SIA'
        h.cfg.study.bl_bmf.loc_flag = h.menu_SPA_loc_flag.Value-1;
        h.cfg.study.bl_bmf.plot_flag=0; h.cfg.study.bl_bmf.text_flag=1;
        h.cfg.study.bl_bmf.perc_crit = 2;
        
        
        %% SIA
        max_sources = str2num(h.edit_SPA_max_sources.String); %floor(size(h.inv_soln(h.next_inv_soln).leadfield.H,1)/5);   % 1 dipole source has 5 degrees of freedom (3 spatial and 2 orientations), thus dividing number of channels/sensors by 5
        [SIA]=BRANElab_MCMV_beamformer_SIA(h.inv_soln(h.next_inv_soln).leadfield.H,[],[],data,...
            h.cfg.study.bl_bmf.act_samps,h.cfg.study.bl_bmf.ctrl_samps,...
            h.inv_soln(h.next_inv_soln).leadfield.voxel_pos,h.cfg.study.bl_bmf.loc_flag,h.cfg.study.bl_bmf.plot_flag,h.cfg.study.bl_bmf.perc_crit,h.cfg.study.bl_bmf.noise_alpha,max_sources);
        h.inv_soln(h.next_inv_soln).Type = 'SIA';
        h.inv_soln(h.next_inv_soln).soln = SIA;
        h.inv_soln(h.next_inv_soln).soln.P.img_org = h.inv_soln(h.next_inv_soln).soln.P.img;
        h.inv_soln(h.next_inv_soln).soln.wts_org = h.inv_soln(h.next_inv_soln).soln.wts;
        h.inv_soln(h.next_inv_soln).soln.ori_org = h.inv_soln(h.next_inv_soln).soln.ori;
        h.inv_soln(h.next_inv_soln).soln.P.img = h.inv_soln(h.next_inv_soln).soln.P.nulled_img; % shifting to nulled img
        h.inv_soln(h.next_inv_soln).soln.wts = h.inv_soln(h.next_inv_soln).soln.nulled_wts;
        h.inv_soln(h.next_inv_soln).soln.ori = h.inv_soln(h.next_inv_soln).soln.nulled_ori;
        
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.menu_SPA_loc_flag.String{h.menu_SPA_loc_flag.Value} ,h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case 'MIA'
        h.cfg.study.bl_bmf.loc_flag = h.menu_SPA_loc_flag.Value-1;
        h.cfg.study.bl_bmf.plot_flag=0; h.cfg.study.bl_bmf.text_flag=1;
        h.cfg.study.bl_bmf.perc_crit = 2;
        
        %% MIA
        max_sources = str2num(h.edit_SPA_max_sources.String); %floor(size(h.inv_soln(h.next_inv_soln).leadfield.H,1)/5);   % 1 dipole source has 5 degrees of freedom (3 spatial and 2 orientations), thus dividing number of channels/sensors by 5
        [MIA]=BRANElab_MCMV_beamformer_MIA(h.inv_soln(h.next_inv_soln).leadfield.H,[],[],data,...
            h.cfg.study.bl_bmf.act_samps,h.cfg.study.bl_bmf.ctrl_samps,...
            h.inv_soln(h.next_inv_soln).leadfield.voxel_pos,h.cfg.study.bl_bmf.loc_flag,h.cfg.study.bl_bmf.plot_flag,h.cfg.study.bl_bmf.perc_crit,h.cfg.study.bl_bmf.noise_alpha,h.cfg.study.bl_bmf.text_flag,h.anatomy,max_sources);
        h.inv_soln(h.next_inv_soln).Type = 'MIA';
        h.inv_soln(h.next_inv_soln).soln = MIA;
        h.inv_soln(h.next_inv_soln).soln.P.img_org = h.inv_soln(h.next_inv_soln).soln.P.img;
        h.inv_soln(h.next_inv_soln).soln.wts_org = h.inv_soln(h.next_inv_soln).soln.wts;
        h.inv_soln(h.next_inv_soln).soln.ori_org = h.inv_soln(h.next_inv_soln).soln.ori;
        h.inv_soln(h.next_inv_soln).soln.wts = h.inv_soln(h.next_inv_soln).soln.nulled_wts;
        h.inv_soln(h.next_inv_soln).soln.P.img = h.inv_soln(h.next_inv_soln).soln.P.nulled_img; % shifting to nulled img
        h.inv_soln(h.next_inv_soln).soln.ori = h.inv_soln(h.next_inv_soln).soln.nulled_ori;

        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.menu_SPA_loc_flag.String{h.menu_SPA_loc_flag.Value} ,h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case 'LCMV'
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
        
        h.inv_soln(h.next_inv_soln).Type = 'LCMV';
        h.inv_soln(h.next_inv_soln).soln = LCMV;
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, cfg.weightnorm ,h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
    case 'eLORETA'
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
                eloreta.P.img = abs(squeeze(nanmean(rms(swf(h.cfg.study.bl_bmf.act_samps,:,:),2))));
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
        
        pow = eloreta.avg.pow;
        eloreta = rmfield(eloreta,{'avg'}); %MNE.avg = rmfield(MNE.avg,{'mom','filter','noisecov'});
        eloreta.avg.pow = single(pow);
        
        h.inv_soln(h.next_inv_soln).Type = 'eLORETA';
        h.inv_soln(h.next_inv_soln).soln = eloreta;
        h.inv_soln(h.next_inv_soln).maxvectorori_Type = h.menu_inv_maxvectorori.String{h.menu_inv_maxvectorori.Value};
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s %s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type,h.menu_inv_maxvectorori.String{h.menu_inv_maxvectorori.Value}, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
        
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
    case 'sLORETA'
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
        
        cfg.grad                = h.inv_soln(h.next_inv_soln).sens;
        cfg.sourcemodel         = h.inv_soln(h.next_inv_soln).leadfield;
        cfg.headmodel           = h.inv_soln(h.next_inv_soln).headmodel;
        sloreta  = ft_sourceanalysis(cfg,avg_data);
        sloreta = ft_sourcedescriptives([], sloreta); % This is gets the nai to get sloreta.P.img

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
        
        ori = cell2mat(sloreta.avg.ori(inside_idx));
        sloreta.ori = ori';
        
        sloreta = rmfield(sloreta,{'avg'}); %MNE.avg = rmfield(MNE.avg,{'mom','filter','noisecov'});
        
        h.inv_soln(h.next_inv_soln).Type = 'sLORETA';
        h.inv_soln(h.next_inv_soln).soln = sloreta;
        h.inv_soln(h.next_inv_soln).ListBox_name = sprintf('%s (%.f-%.f ms) %s',h.inv_soln(h.next_inv_soln).Type, h.inv_soln(h.current_inv_soln).params.act_int*1000,h.inv_soln(h.current_inv_soln).headmodel_type);
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
    case 'MNE'
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
                MNE.P.img = abs(squeeze(nanmean(rms(swf(h.cfg.study.bl_bmf.act_samps,:,:),2))));
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
        h.inv_soln(h.next_inv_soln).Type = 'MNE';
        h.inv_soln(h.next_inv_soln).soln = MNE;
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
    case {'dics'}   % requires Time-Freq do have been computed
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
    case {'pcc'}
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
 end

h.inv_soln(h.next_inv_soln).soln.compute_time = toc(tt0); % time to compute inv_soln

%% Listbox Name

%% original map that gets changed with spatiotemporal mapping
h.inv_soln(h.current_inv_soln).org_img = h.inv_soln(h.current_inv_soln).soln.P.img;

min_max = round([min(h.inv_soln(h.current_inv_soln).soln.P.img) max(h.inv_soln(h.current_inv_soln).soln.P.img)],3,'significant');
if min_max(1)==min_max(2); min_max(2)=min_max(1)+1;
elseif min_max(1)>min_max(2); min_max(2)=min_max(1)+1;
end

h.inv_soln(h.current_inv_soln).soln.plot_min_max = real(min_max);
h.inv_soln(h.current_inv_soln).soln.plot_thresh = round(h.inv_soln(h.current_inv_soln).soln.plot_min_max(2)*.5,3,'significant'); %h.inv_soln(h.current_inv_soln).soln.null_thresh;
h.edit_3D_min_max.String = sprintf('%.3f %.3f',h.inv_soln(h.current_inv_soln).soln.plot_min_max);

% updating threshold slider
h.slider_3D_image_thresh.Max = h.inv_soln(h.current_inv_soln).soln.plot_min_max(2);
h.slider_3D_image_thresh.Min = h.inv_soln(h.current_inv_soln).soln.plot_min_max(1);
h.slider_3D_image_thresh_text_max.String = num2str(h.slider_3D_image_thresh.Max);
if h.slider_3D_image_thresh.Value>h.slider_3D_image_thresh.Max
    h.slider_3D_image_thresh.Value = h.slider_3D_image_thresh.Max;
elseif h.slider_3D_image_thresh.Value < h.slider_3D_image_thresh.Min
    h.slider_3D_image_thresh.Value = h.slider_3D_image_thresh.Min;
end

bs_update_3D_listbox();
bs_plot_inv_soln;
fprintf('Inverse Modeling Completed for %s\n',inv_soln);

%% converting to single precision for reduce memory storage
h.inv_soln(h.next_inv_soln).soln.wts = single(h.inv_soln(h.next_inv_soln).soln.wts);
h.inv_soln(h.next_inv_soln).soln.ori = single(h.inv_soln(h.next_inv_soln).soln.ori);
h.inv_soln(h.next_inv_soln).soln.ori = single(h.inv_soln(h.next_inv_soln).soln.ori);
h.inv_soln(h.next_inv_soln).soln.P.img = single(h.inv_soln(h.next_inv_soln).soln.P.img);
% reducing redundancy to save storage space
if isfield(h.inv_soln(h.next_inv_soln).leadfield,'cfg')
    h.inv_soln(h.next_inv_soln).leadfield = rmfield(h.inv_soln(h.next_inv_soln).leadfield,{'cfg'});
end
% if isfield(h.inv_soln(h.next_inv_soln).leadfield,'leadfield')
%     h.inv_soln(h.next_inv_soln).leadfield = rmfield(h.inv_soln(h.next_inv_soln).leadfield,{'leadfield'});
% end
% if isfield(h.inv_soln(h.next_inv_soln).leadfield,'H')
%     h.inv_soln(h.next_inv_soln).leadfield = rmfield(h.inv_soln(h.next_inv_soln).leadfield,{'H'});
% end

h.inv_soln(h.next_inv_soln).leadfield = double2single(h.inv_soln(h.next_inv_soln).leadfield);
h.inv_soln(h.current_inv_soln).TFR_results =[];

% catch me
%     display(me)
%     errordlg('ERROR! Inverse Source Modeling Failed');
% end
if h.monte_carlo_flag ~= 1
    h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
end
% close(hm);


