% function sm_batch_load_gavg_results(study_name,data_dir,inv_soln_fieldnames)

sim_data_fieldnames = {'cfg.source.vx_idx' 'cfg.source.vx_locs' ...
    };

inv_soln_fieldnames = {'ori_error_mse' 'wave_error_mse_norm_act' 'wave_error_mse_norm_ctrl' 'half_width' ...
    'soln.P.img' ...
    'classifier_metrics.true_idx' 'classifier_metrics.loc_error' 'classifier_metrics.Hits' 'classifier_metrics.FA' 'classifier_metrics.Miss' 'classifier_metrics.num_CR'...
     'classifier_performance.MCC' 'classifier_performance.TPR' 'classifier_performance.FPR' ...
    'classifier_results.stepwise_thresholded_MCC' 'classifier_results.stepwise_thresholded_TPR' 'classifier_results.stepwise_thresholded_FPR' ...
    'classifier_results.stepwise_thresholded_num_hits' 'classifier_results.stepwise_thresholded_num_miss' 'classifier_results.stepwise_thresholded_num_fa' 'classifier_results.stepwise_thresholded_num_CR' ...
    'classifier_results.stepwise_thresholded_loc_error' 'classifier_results.stepwise_thresholded_ori_error_mse' 'classifier_results.stepwise_thresholded_wave_error_mse_norm_act' 'classifier_results.stepwise_thresholded_wave_error_mse_norm_ctrl' ...
    'classifier_results.stepwise_thresh_vals' };


data_dir = 'C:\Data\SimMEEG_studies\Exp3_Synth_VEP_brain_noise\';
study_name = 'Exp3_Synth_VEP_brain_noise';
anat_dir = 'C:\Data\SimMEEG_studies\';
anat_name = 'SimMEEG_studies_default_Anatomy.mat';
anatomy_file = fullfile(anat_dir,anat_name);

study_file = [study_name '_study_parameters.mat'];

%% load parameter file
load(fullfile(data_dir, study_file),'monte_params');
%% load sens_meg sens_eeg from anatomy file
load(anatomy_file,'sens_meg','sens_eeg'); 


%% load datafiles
% datafiles = dir(fullfile(data_dir,'*Run*.mat'));
num_runs = monte_params.num_sims;
num_trials = monte_params.num_trials;
num_reruns = monte_params.num_reruns;
start_run = str2num(monte_params.h.edit_monte_run_num.String);

for nr = 1:num_reruns
    for r = start_run:start_run+num_runs     % Looping through number of runs
        
        for nt = 1:length(num_trials) % numbr of trials
            for p = 1:length(monte_params.plv_range)  % this level needs to run "Sim Source Data" btn and simulate new source waveforms
                %% LOCATION Change from initial source locations based on random norm distribution for parameters set in Monte Carlo Tab by "Location Range (min:step:max)" and "Std Dev"
                for lc = 1:size(monte_params.source_loc_range_X,2) % looping through "Location Range (min:step:max)"
                    %% ORIENTATION Change from initial source locations based on random norm distribution for parameters set in Monte Carlo Tab by "Orientation Range (min:step:max)" and "Std Dev"
                    for m = 1:size(monte_params.source_ori_range_Az,2) % looping through "Location Range (min:step:max)"
                        %% AMPLITUDE Values based on random norm distribution for parameters set in Monte Carlo Tab by "Amplitude Range (min:step:max)" and "Std Dev"
                        for a=1:size(monte_params.source_amp_range,2) % looping through "Amplitude Range (min:step:max)"
                            %% SNR based on random norm distribution for parameters set in Monte Carlo Tab by "Location Range (min:step:max)" and "Std Dev"
                            for s = 1:length(monte_params.SNR_range) % looping through "SNR Range (min:step:max)"
                                
                                %% SENSOR Type
                                for st = 1:length(monte_params.sens_type) % SENSOR Type
                                     if monte_params.sens_type(st)==1 % MEG sensors selected
                                         mont_str = cellstr(monte_params.h.listbox_monte_MEG_sens_montage.String);
                                        sens_montage = mont_str(monte_params.h.listbox_monte_MEG_sens_montage.Value');
                                        
                                    elseif monte_params.sens_type(st)==2 % EEG sensors selected
                                        mont_str = cellstr(monte_params.h.listbox_monte_EEG_sens_montage.String);
                                        sens_montage = mont_str(monte_params.h.listbox_monte_EEG_sens_montage.Value');
                                     end
                                    
                                     for sm = 1:length(sens_montage) % number of sensors
                                         num_sens = strtrim(char(sens_montage{sm}));                                         
                                         amp_mu    = monte_params.source_amp_range(:,a);   % "Amplitude Range (min:step:max)"
                                         locX_mu    = monte_params.source_loc_range_X(:,lc);   % "Location Range (min:step:max)"
                                         locY_mu    = monte_params.source_loc_range_Y(:,lc);   % "Location Range (min:step:max)"
                                         locZ_mu    = monte_params.source_loc_range_Z(:,lc);   % "Location Range (min:step:max)"
                                         oriAz_mu    = monte_params.source_ori_range_Az(:,m);   % "Orientation Range (min:step:max)"
                                         oriEl_mu    = monte_params.source_ori_range_El(:,m);   % "Orientation (min:step:max)"
                                         snr_mu    = monte_params.SNR_range(s);   % "SNR Range (min:step:max)"
                                         
                                        noise_types = {'SynthNoise' 'BrainNoise' 'GANNoise' 'RealSensorNoise' 'RealRestingNoise'};
                                        if monte_params.h.menu_monte_synthetic_real_data.Value <5  % 'Synthetic Data' - simulate data from user-defined singals+prepost
                                            xsname = sprintf('%s_%s%s_%s_SourceAmps_%.f_%.f_%.f_LocX_%.f_%.f_%.f_LocY_%.f_%.f_%.f_LocZ_%.f_%.f_%.f_Az_%.f_%.f_%.f_El_%.f_%.f_%.f_SNR_%.3f_PLV_%.3f_Trials_%s_Run_%.f.mat',...
                                                monte_params.h.edit_study_name.String,monte_params.h.menu_sens_type.String{monte_params.sens_type(st)},num_sens,noise_types{monte_params.h.menu_monte_synthetic_real_data.Value},...
                                                amp_mu,locX_mu,locY_mu,locZ_mu,oriAz_mu,oriEl_mu,snr_mu,monte_params.plv_range(p),monte_params.h.edit_num_trials.String,r);
                                            sname = fullfile(data_dir,xsname);
%                                         elseif h.menu_monte_synthetic_real_data.Value == 5  % 'Real Sensors or Sources' - load real sensors data and randomize phase + load in source waveform data
%                                             [~,real_fname,~]= fileparts(h.monte_real_sensor_files{r});
%                                             xsname = sprintf('%s_%s%s_SourceAmps_%.f_%.f_%.f_LocX_%.f_%.f_%.f_LocY_%.f_%.f_%.f_LocZ_%.f_%.f_%.f_Az_%.f_%.f_%.f_El_%.f_%.f_%.f_SNR_%.3f_PLV_%.3f_Trials_%s_Run_%.f.mat',...
%                                                 real_fname,h.menu_sens_type.String{h.monte_params.sens_type(st)},num_sens,noise_types{h.menu_monte_synthetic_real_data.Value},amp_mu,locX_mu,locY_mu,locZ_mu,oriAz_mu,oriEl_mu,snr_mu,h.monte_params.plv_range(p),h.edit_num_trials.String,h.monte_params.sim_run_num);
%                                             sname = fullfile(h.data_dir,xsname);
                                        else
                                        end
                                        
                                        %% Loading dataset
                                        if exist(sname,'file')
                                            fprintf('Loading:    %s\n',xsname);
                                            load(sname);
                                            
                                            % only loading data based on fieldnames
                                            gavg_results(nr, r, nt, p, lc, m, a, s, st, sm).sim_data = sm_batch_load_results(sim_data,sim_data_fieldnames);
                                            gavg_results(nr, r, nt, p, lc, m, a, s, st, sm).inv_soln = sm_batch_load_results(inv_soln,inv_soln_fieldnames);
                                          
                                            
                                        else
                                            fprintf('File does NOT exist:    %s\n',xsname);
                                        end
                                        
                                     end
                                end
                            end
                        end
                    end
                end
            end
        end     % num_trials
    end
end



%% 

