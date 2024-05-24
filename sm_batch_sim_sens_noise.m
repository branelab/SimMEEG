function h = sm_batch_sim_sens_noise(h)

h.sim_data.sens_noise_type = h.menu_monte_synthetic_real_data.String{h.menu_monte_synthetic_real_data.Value};


if h.menu_noise_projection.Value == 1   % Simulate Sensor noise for all sensors
    %% NOTE: Simulating sensor noise is quite complicated because of spatial and temporal covariance needs to be modeled --> see reading in GUI directory: DeMunckEtAl_StatDip_KronPrud_IEEE_SP_2002_spatiotemporal_noise_covariance_simulation.pdf
    % Although not that accurate at representing "real" spatiotemporal noise - the following code attempts to do so.
    % this gives at least a better approximation that uncorrelated gaussian noise that is often used in simulations.
    if h.menu_synthetic_noise_cov_type.Value==1 % no shaping by spatial or temporal covariance
        h.sim_data.sens_noise = sim_noise_wav(h.cfg,size(h.anatomy.leadfield.H,1));
        h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
        %         for t=randperm(h.monte_params.cfg.study.num_trials); figure(3); clf; surf(cov(squeeze(h.sim_data.sens_noise_scaled(:,:,t(1))))); view(0,90); shading interp; axis tight; caxis([-.045 .045]); pause; end
    elseif h.menu_synthetic_noise_cov_type.Value==2     % Spatial Covariance shaping only
        
        noiseVec2 = sim_noise_wav(h.cfg,size(h.anatomy.leadfield.H,1));
        h.monte_params.cfg.study.sensor_noise_cov_exp = str2num(h.edit_synthetic_noise_cov_exp.String);  %  see citation below
        %         noiseVec2 = noiseVec2/max(max(max(abs(noiseVec2))));
        noiseVec = permute(noiseVec2,[2 1 3]);
        for t=1:size(noiseVec,3)     % randomizing covariance structure across trials
            c_cov = random('exp',h.monte_params.cfg.study.sensor_noise_cov_exp,size(noiseVec,1),size(noiseVec,1)); % spatial covariances 'exp' (1/3.57^1.02) not accounting for distance between channels so loosely based on Huizenga et al., 2002 IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. 49, NO. 6, pp. 533:539; see Huizenga_2002_noise_cov_randomization_exp_values.pdf
            noiseVec(:,:,t) = c_cov*squeeze(noiseVec(:,:,t));
            r_cov(:,:,t) = c_cov;
        end
        h.sim_data.sens_noise = permute(noiseVec,[2 1 3]);
        h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
    elseif h.menu_synthetic_noise_cov_type.Value==3  || h.menu_synthetic_noise_cov_type.Value==4   % Temporal or Spatiotemporal shaping using ARM
        %         answ = questdlg(sprintf('Temporal & Spatiotemporal shaping of Noise using ARM can take several minutes depending on ARM order and # interactions?\n\nYou can reduce ARM order which speeds up computations but at a cost of modeling.\n\nDo you want to Continue?\n'),'Simulate Spatiotemporal Noise?','Yes','No','No');
        answ = 'Yes';
        switch answ
            case 'No'
            case 'Yes'
                
                %% Implemented using ARM simulations but see --> "spatio-temporal noise covariance is a Kronecker product of a spatial and a temporal covariance matrix"
                num_samps = h.monte_params.cfg.study.num_samps;
                num_chans = size(h.anatomy.leadfield.H,1);
                num_trials = h.monte_params.cfg.study.num_trials;
                
                if h.menu_synthetic_noise_cov_type.Value==3 % Only Temporally shaped
                    num_interactions = 0;  % only one interaction among sensors thus no spatial correlations by design. There still could be some but these should be random.
                    blk_trials = num_trials; % running blocks of trials so that ARM interactions aren't always among same electrodes for every trial.
                    shape_type = 'Temporally';
                elseif h.menu_synthetic_noise_cov_type.Value==4 % Spatiotemporally shaped
                    shape_type = 'SpatioTemporally';
                    num_interactions = str2num(h.edit_synthetic_noise_ARM_interaction.String);
                    %                     num_blks = ceil(num_trials/num_interactions);
                    num_blks = ceil(num_trials/10);
                    blk_trials = ceil(num_trials/num_blks); % running blocks of trials so that ARM interactions aren't always among same electrodes for every trial.
                    %                     blk_trials = num_trials; % running blocks of trials so that ARM interactions aren't always among same electrodes for every trial.
                end
                
                ARM_order = str2num(h.edit_synthetic_noise_ARM_order.String);
                noise_wav = [];
                %                 hm = msgbox(sprintf('Using ARM order = %.f to reduce modeling time for %.f sensors',ARM_order,num_chans),'');
                hw = waitbar(0,sprintf('Creating %s-Shaped Noise',shape_type)); hw.Units = 'normalized'; hw.Position(2) = hw.Position(2)+.2;
                %% adding noise
                %                 for blk=1:num_blks
                k=0;
                while k==0  % running until neough trials within noise_wav after removing trials that may have ARM artefacts
                    cfg = h.cfg;
                    pad_samps = num_samps*10; % initial padded samples needed for ARM calculations - added as additional trials
                    cfg.study.num_trials = blk_trials+10; % padding trials for passing samples because matrix is reshaped num_samps*num_trials
                    % changing cfg based on noise params from Panel "Study Parameters" not from Panel "Generate Sensor Noise".
                    freqs = str2num(h.edit_synthetic_noise_freqs.String);
                    sig_freqs = repmat(freqs,[num_chans 1 1]); cfg.source.sig_freqs = permute(sig_freqs,[1 3 2]);
                    cfg.study.synthetic_noise_flag = h.menu_synthetic_noise_flag.Value;
                    cfg.study.synthetic_noise_freqs = str2num(h.edit_synthetic_noise_freqs.String);
                    cfg.study.synthetic_pink_noise_slope = h.monte_params.cfg.study.synthetic_pink_noise_slope;
                    
                    noiseVec2 = sim_noise_wav(cfg, num_chans);
                    
                    dims = size(noiseVec2);
                    noiseVec2 = reshape(permute(noiseVec2,[2 1 3]),[dims(2) dims(1)*dims(3)]);
                    [noiseVec2] = arm_generate_signal(num_chans, num_samps, ARM_order, num_interactions, pad_samps, noiseVec2);
                    % removing trials that have > overall variance than the rest
                    % will be removed;
                    noiseVec2 = permute(reshape(noiseVec2,[num_chans, num_samps, blk_trials]),[2 1 3]);
                    xvar = squeeze(var(var(noiseVec2))); good_trials = ~isoutlier(xvar);
                    fprintf('%.f of %.f good trials\n',sum(good_trials),length(good_trials));
                    noiseVec2 = noiseVec2(:,:,good_trials);
                    
                    noise_wav = cat(3,noise_wav,noiseVec2);
                    waitbar(size(noise_wav,3)/num_trials,hw, sprintf('Creating %s-Shaped Noise (%.f of %.f) Trials',shape_type,size(noise_wav,3),num_trials));
                    if size(noise_wav,3)>num_trials; k=1; break; end
                    
                end
                dims = size(noise_wav);
                noise_wav = noise_wav(1:num_samps,1:num_chans,1:num_trials);
                close (hw);
                h.sim_data.sens_noise = noise_wav;
                h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
                %                 for t=1:size(h.sim_data.sens_noise_scaled,3); figure(3); clf;
                %                     subplot(1,3,1); plot(cfg.study.lat_sim, squeeze(h.sim_data.sens_noise_scaled(:,:,t(1))),'k'); title(sprintf('noise waves for trial %.f'),t)
                %                     subplot(1,3,2); surf(cov(squeeze(h.sim_data.sens_noise_scaled(:,:,t)))); view(0,90); shading interp; axis tight; caxis([-.045 .045]); title('Spatial (Sensor) Covariance');
                %                     subplot(1,3,3); surf(cov(squeeze(h.sim_data.sens_noise_scaled(:,:,t)'))); view(0,90); shading interp; axis tight; caxis([-.15 .15]); title('Temporal Covariance');
                %                     pause;
                %                 end
        end
        
    elseif h.menu_synthetic_noise_cov_type.Value==5     % Real Spatial Covariance shaping
        
        fprintf('Real Spatial Covariance shaping not implemented in batch script yet. PLease choose another shaping type\n');
        %% Real Spatial Covariance shaping
        % [cov_file, fpath] = uigetfile('*.mat','Load Real Data for Covariance shaping');
        %         load(fullfile(fpath,cov_file),'data'); % loading data from pre-created file with data struct = [chan x samps x trial];
        %
        %         noiseVec2 = sim_noise_wav(h.cfg,size(h.anatomy.leadfield.H,1));
        %         if size(data,1)<size(noiseVec2,2) || size(data,2)<size(noiseVec2,1) || size(data,3)<size(noiseVec2,3)
        %             dims = size(data);
        %             dims2 = size(noiseVec2);
        %             txt = sprintf('Size of Loaded Data has insufficient data to simulate noise\nDimensions:\n   Data = [%.f chans x %.f samps x %.f trials]\n   Noise = [%.f chans x %.f samps x %.f trials]',dims,dims2([2 1 3]));
        %             warndlg(txt);
        %             h.sim_data.sens_noise = [];
        %             h.sim_data.sens_noise_scaled = [];
        %         else
        %             %         noiseVec2 = noiseVec2/max(max(max(abs(noiseVec2))));
        %             noiseVec = permute(noiseVec2,[2 1 3]);
        %             n_rand = randperm(size(data,3)); % randomly selecting data trials
        %             fprintf('Randomly selected %.f trials from data with %.f trials\n',size(noiseVec2,3),size(data,3))
        %             data = data(:,:,n_rand(1:size(noiseVec2,3)));
        %             for t=1:size(noiseVec,3)     % randomizing covariance structure across trials
        %                 %             c_cov = random('exp',h.monte_params.cfg.study.sensor_noise_cov_exp,size(noiseVec,1),size(noiseVec,1)); % spatial covariances 'exp' (1/3.57^1.02) not accounting for distance between channels so loosely based on Huizenga et al., 2002 IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. 49, NO. 6, pp. 533:539; see Huizenga_2002_noise_cov_randomization_exp_values.pdf
        %                 c_cov = cov(squeeze(data(1:length(h.anatomy.sens.good_sensors),h.monte_params.cfg.study.base_samps,t))');
        %                 noiseVec(:,:,t) = c_cov*squeeze(noiseVec(:,:,t));
        %                 r_cov(:,:,t) = c_cov;
        %             end
        %             h.sim_data.sens_noise = permute(noiseVec,[2 1 3]);
        %             h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
        %             %         for t=randperm(h.monte_params.cfg.study.num_trials); figure(3); clf; surf(cov(squeeze(h.sim_data.sens_noise_scaled(:,:,t(1))))); view(0,90); shading interp; axis tight; caxis([-.00045 .00045]); pause; end
        %         end
    end
    
    
elseif h.menu_noise_projection.Value == 2   % Simulate Brain noise
    if ~isfield(h.monte_params.cfg.source,'brain_noise_idx') % in case "Random locs" button has yet to be pressed
        btn_rand_brain_noise_locs_Callback([],[]);
    end
    %% code below replaces "sim_brain_noise_wav.";
    num_chans = str2num(h.edit_num_noise_sources.String);
    % broadband noise
    h.cfg.study.noise_flag = h.menu_synthetic_noise_flag.Value;    % 1=Broadband, 2=NarrowBand, 3=Notched Band, 4=Pink
    h.sim_data.brain_noise = sim_noise_wav(h.monte_params.cfg,num_chans);
    h.sim_data.brain_noise = h.sim_data.brain_noise/max(max(max(abs(h.sim_data.brain_noise))));

    
    
    %% projecting brain source locations across trials
    h.sim_data.sens_noise=[];
    for t=1:h.monte_params.cfg.study.num_trials
        h.sim_data.sens_noise(:,:,t) = sm_batch_project_SimSignals(h.sim_data.brain_noise(:,:,t),h.anatomy.leadfield,h.monte_params.cfg.source.brain_noise_idx(:,t),h.monte_params.cfg.source.brain_noise_amp(:,t),h.monte_params.cfg.source.brain_noise_ori(:,:,t),h);
    end
    h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
    % validating code
    %     figure(1); clf; hold on; for t=1:h.monte_params.cfg.study.num_trials; plot(h.monte_params.cfg.study.lat_sim,squeeze(h.sim_data.sens_noise_scaled(:,:,t)),'k'); end
    %     figure(2); clf; calc_fft(squeeze(h.sim_data.sens_noise(:,48,:)),256,1,'k');
    %         for t=randperm(h.monte_params.cfg.study.num_trials); figure(3); clf; surf(cov(squeeze(h.sim_data.sens_noise_scaled(:,:,t(1))))); view(0,90); shading interp; axis tight; caxis([-.025 .025]); pause; end
    %     scatter3(h.axes_anatomy,h.anatomy.leadfield.voxel_pos(h.monte_params.cfg.source.brain_noise_idx,1),h.anatomy.leadfield.voxel_pos(h.monte_params.cfg.source.brain_noise_idx,2),h.anatomy.leadfield.voxel_pos(h.monte_params.cfg.source.brain_noise_idx,3),'k.');
elseif h.menu_noise_projection.Value == 3       % Generative Adversarial Networks (GAN) generated noise <-- currently must have srate=256 and duration -1 to 1 sec with chans = 1:66
    implemented_flag = 0;    % Not implmented yet
    if implemented_flag ==1
        rest_dir = 'C:\BRANELab\matlab_progs\general_progs\EEG_sim\SimSignals_GUI\GANsimEEG\Trained_GANs';
        sname = fullfile(rest_dir,sprintf('Tujillo_Resting_State_%s_GAN_trained.mat',blk_name{h.menu_resting_state.Value}));
        fprintf('GAN generated noise using:  %s\n',sname);
        load(sname,'dlnetGenerator');
        
        %% Generate EEG sensor noise using pre-trained GAN
        
        executionEnvironment = "auto";
        ZZNew = randn(1,1,1000,h.monte_params.cfg.study.num_trials,'single');
        dlZZNew = dlarray(ZZNew,'SSCB');
        if (executionEnvironment == "auto" && canUseGPU) || executionEnvironment == "gpu"
            dlZZNew = gpuArray(dlZZNew);
        end
        
        dlXGeneratedNew = predict(dlnetGenerator,dlZZNew);
        y = gather(extractdata(dlXGeneratedNew));
        ydata = permute(squeeze(y(:,:,:,:)),[2 1 3]);
        ydata = ydata(1:size(h.anatomy.sens.label),:,:);
        ydata = bl_reref(ydata,[1:size(ydata,1)],[],1);    % average referencing
        h.sim_data.sens_noise = permute(ydata,[2 1 3]);
        h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
    end
elseif h.menu_noise_projection.Value == 4        % Real Sensor Noise
    
end
fprintf('Finished Simulating "%s" Noise\n',h.menu_noise_projection.String{h.menu_noise_projection.Value});
%% simulate Brain Noise
