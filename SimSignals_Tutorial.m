% SimSignals_Tutorial_v2 -includes amp-phase interactions

%% %%% Example #1 Simulation Parameters for two frequency (10 & 18 Hz) signals, no Phase-Ampltidue Couplings %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Study Parameters
clear cfg;
cfg.study.srate=256;
cfg.study.dur=[-2 2]; % (sec) start and end times for whole trial
cfg.study.lat_sim=[cfg.study.dur(1):1/cfg.study.srate:cfg.study.dur(2)-(1/cfg.study.srate)]; % latency of each trial
cfg.study.num_samps=length(cfg.study.lat_sim);
cfg.study.num_trials=90;
cfg.study.max_perm_plv=9; % integer of cfg.study.num_trials/10 = maximum number of permutations to search for PLV_trials (must be <10 or memory will fail on most computers) .
cfg.study.noise_flag=3;  % Whitening noise to be added to each sources
cfg.study.noise_amp_perc=100; % percent of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
cfg.study.noise_freqs=[1 100]; % [1 100] Frequency (Hz) of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
cfg.study.plv_thresh=.05;   % stoppping criterion when search for best PLV/PLI matched to sig_PLV_targets, sig_PLI_targets, etc. (e.g., 0.05).
cfg.study.plot_sim_flag=0;  % plots evoked, trial, and PLV results. Note: This can take time because of filtering.
cfg.study.plot_time_int=[-0.4 1];   % time interval to plot
cfg.study.plot_freq_int=[5 30]; % frequencies for calculating and plotting PLV/PLI
cfg.study.base_int=[-0.4 0];    % base line interval for plotting
%% Signal Parameters
cfg.source.sig_freqs(1,:,:) = [10 10; 18 18];  % dipole#1 (Nfreqs x [minfreq maxfreq]); sinusoids will be randomly assigned between [minfreq maxfreq] across num_trials;
cfg.source.sig_freqs(2,:,:) = [10 10; 18 18];  % dipole#2 (Nfreqs x [minfreq maxfreq]); sinusoids will be randomly assigned between [minfreq maxfreq] across num_trials;
cfg.source.sig_freqs(3,:,:) = [10 10; 18 18];  % dipole#3  (Nfreqs x [minfreq maxfreq]); sinusoids will be randomly assigned between [minfreq maxfreq] across num_trials;
cfg.source.sig_amp_perc = [70 30; 70 30; 70 30];   % 3sigs x Nfreqs % defines amount of amplitude within the signal duration (sig_dur) for the peak of the Hanning window across sig_dur.
cfg.source.prepost_amp_perc = [30 70; 30 70; 30 70];   % 3sigs x Nfreqs % defines amount of amplitude within the signal duration (sig_dur) for the peak of the Hanning window across sig_dur.
cfg.source.sig_amp_perc_std = [0 0; 0 0; 0 0];   % 3sigs x Nfreqs % standard deviation (as a percent value) of the signal amplitude (sig_amp_perc) across trials
cfg.source.prepost_amp_perc_std = [0 0; 0 0; 0 0];   % 3sigs x Nfreqs % standard deviation (as a percent value) of the prepost amplitude (sig_amp_perc) across trials

cfg.source.sig_evoked_perc = [50 50; 50 50; 50 50]; %(3signals x Nfreqs); % based on signal phases so that the evoked (averaged amplitude) signal
cfg.source.prepost_evoked_perc = [0 0; 0 0; 0 0];  %(3signals x Nfreqs); % based on signal phases so that the evoked (averaged amplitude) signal
cfg.source.sig_durs = [.6 .6; .6 .6; .6 .6];  %(3signals x Nfreqs);    % duration of Hanning windowing of signals amplitude defined by sig_perc .
cfg.source.sig_start = [.2 .2; .2 .2; .2 .2];  %(3signals x Nfreqs);   % start (onset) time of the signals relative to 0=sample(1).
cfg.source.sig_win_type = [3 3; 3 3; 3 3];  %(3signals x Nfreqs) % type of windowing function (1)='Hann' (2)='Gauss' (3)='Triang'  (4)='Blackman' ;
cfg.source.sig_win_rise_time=[.1 .1; .1 .1; .1 .1]; %(3signals x Nfreqs) % rise_time of windowing function. If less than 1/2 sig_duration then window function will have a plateau of the difference in duration between the rise/fall time and the signal duration;

cfg.source.sig_PLV_targets = [.4 0; .4 0; .4 0];  %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
cfg.source.prepost_PLV_targets = [0 0.4; 0 0.4; 0 0.4];  %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
cfg.source.sig_PLI_targets = [.4 0; .4 0; -.4 0];   %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
cfg.source.prepost_PLI_targets = [0 0.4; 0 0.4; 0 -.4]; %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
cfg.source.sig_phase_lag = ([0 0; 0 0; 0 0]/360)*2*pi;  % phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval relative to first sample --> (phase_lag/360)*2*pi);    cos(phase_lag) = correlation of signal relative to zero-phase onset
cfg.source.prepost_phase_lag =([0 0; 0 0; 0 0]/360)*2*pi;  % phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval relative to first sample --> (phase_lag/360)*2*pi);    cos(phase_lag) = correlation of signal relative to zero-phase onset

cfg.source.phase_amp_contrasts = [1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; % sig_contrasts = cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
cfg.source.sig_phase_amp_freq_idx = [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; %(sig_contrasts x Nfreqs);  % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
cfg.source.prepost_phase_amp_freq_idx = [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ];   %(PLV_contrasts x Nfreqs);     % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
cfg.source.sig_phase_amp_depth_perc = [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; % amplitude-modulation depth as a percentage of signal's amplitude (sig_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
cfg.source.prepost_phase_amp_depth_perc = [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; % depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
cfg.source.sig_phase_amp_depth_perc_range = [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; % +/- range of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: sig_phase_amp_depth_perc +/- sig_phase_amp_depth_perc_range must be within [0 100]
cfg.source.prepost_phase_amp_depth_perc_range = [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; % +/- range devitaion of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: prepost_phase_amp_depth_perc +/- prepost_phase_amp_depth_perc_range must be within [0 100]
%% Run simulation 
run_multiple=1; % (0) run simulation only once and replicate simulation num_iter (1) run simulation multiple times = num_iter
num_iter=10;

% replicate same data to get many trials so that noise phase can be added in other frequency bands and cancel out across more trials for TFR analyses.
[sig_final,sig_wav,prepost_wav,noise_wav,cfg,prepost_win,sig_win] = SimSignals(cfg);
sfg=repmat(sig_final,[1 1 1 num_iter]);
sig_phase=repmat(cfg.sig_phase,[1 1 1 num_iter]);
prepost_phase=repmat(cfg.prepost_phase,[1 1 1 num_iter]);
sig_final=reshape(sfg,[size(sfg,1) size(sfg,2) size(sfg,3)*size(sfg,4)]); sig_final_org=sig_final; % reshape back to single iteration run

% OR
if run_multiple==1
    % Run multiple iterations to yield more trials and so that noise phase can be added in other frequency bands and cancel out across more trials for TFR analyses
    for tt=1:num_iter
        t1=tic;
        fprintf('Iteration #%.f\n',tt);
        [sig_final,sig_wav,prepost_wav,noise_wav,cfg,prepost_win,sig_win] = SimSignals(cfg);
        p_time(tt)=toc(t1);
        sfg(:,:,:,tt)=sig_final;
        cfg_tt(tt)=cfg;
        sig_phase(:,:,:,tt)=cfg_tt(tt).sig_phase; prepost_phase(:,:,:,tt)=cfg_tt(tt).prepost_phase;
    end
    dims=size(sig_phase); sg=permute(reshape(permute(sig_phase,[1 3 2 4]),[dims(1) dims(3) dims(2)*dims(4)]),[1 3 2]);
    dims=size(prepost_phase); pg=permute(reshape(permute(prepost_phase,[1 3 2 4]),[dims(1) dims(3) dims(2)*dims(4)]),[1 3 2]);
end


%% %%% Example #2 Simulation Parameters for three frequency (6, 22, & 18 Hz) signals, with Phase-Ampltidue Couplings %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Study Parameters
clear cfg;
cfg.study.srate=256;
cfg.study.dur=[-2 2]; % (sec) start and end times for whole trial
cfg.study.lat_sim=[cfg.study.dur(1):1/cfg.study.srate:cfg.study.dur(2)-(1/cfg.study.srate)]; % latency of each trial
cfg.study.num_samps=length(cfg.study.lat_sim);
cfg.study.num_trials=90;
cfg.study.max_perm_plv=9; % integer of cfg.study.num_trials/10 = maximum number of permutations to search for PLV_trials (must be <10 or memory will fail on most computers) .
cfg.study.noise_flag=0;  % Whitening noise to be added to each sources
cfg.study.noise_amp_perc=[]; % percent of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
cfg.study.noise_freqs=[]; % [1 100] (or [1:100] for brain noise) Frequency (Hz) of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
cfg.study.plv_thresh=.05;   % stoppping criterion when search for best PLV/PLI matched to sig_PLV_targets, sig_PLI_targets, etc. (e.g., 0.05).
cfg.study.plot_sim_flag=0;  % plots evoked, trial, and PLV results. Note: This can take time because of filtering.
cfg.study.plot_time_int=[-0.4 1];   % time interval to plot
cfg.study.plot_freq_int=[3 56]; % frequencies for calculating and plotting PLV/PLI
cfg.study.base_int=[-0.4 0];    % base line interval for plotting
%% Signal Parameters
cfg.source.sig_freqs(1,:,:) = [6 6; 22 22; 40 40]; % dipole#1 (Nfreqs x [minfreq maxfreq]); sinusoids will be randomly assigned between [minfreq maxfreq] across num_trials;
cfg.source.sig_freqs(2,:,:) = [6 6; 22 22; 40 40]; % dipole#2 (Nfreqs x [minfreq maxfreq]); sinusoids will be randomly assigned between [minfreq maxfreq] across num_trials;
cfg.source.sig_freqs(3,:,:) = [6 6; 22 22; 40 40]; % dipole#3  (Nfreqs x [minfreq maxfreq]); sinusoids will be randomly assigned between [minfreq maxfreq] across num_trials;
cfg.source.sig_amp_perc = [70 20 70; 70 20 70; 70 20 70];   % 3sigs x Nfreqs % defines amount of amplitude within the signal duration (sig_dur) for the peak of the Hanning window across sig_dur.
cfg.source.prepost_amp_perc = [30 60 30;30 60 30; 30 60 30];   % 3sigs x Nfreqs % defines amount of amplitude within the signal duration (sig_dur) for the peak of the Hanning window across sig_dur.
cfg.source.sig_amp_perc_std = [5 5 5; 5 5 5; 5 5 5];   % 3sigs x Nfreqs % standard deviation (as a percent value) of the signal amplitude (sig_amp_perc) across trials
cfg.source.prepost_amp_perc_std = [5 5 5; 5 5 5; 5 5 5];   % 3sigs x Nfreqs % standard deviation (as a percent value) of the prepost amplitude (sig_amp_perc) across trials

cfg.source.sig_evoked_perc=[0 0 0; 0 0 0; 0 0 0]; %(3signals x Nfreqs); % based on signal phases so that the evoked (averaged amplitude) signal
cfg.source.prepost_evoked_perc=[0 0 0; 0 0 0; 0 0 0]; %(3signals x Nfreqs); % based on signal phases so that the evoked (averaged amplitude) signal
cfg.source.sig_durs=[.6 .6 .6; .6 .6 .6; .6 .6 .6]; %(3signals x Nfreqs);    % duration of Hanning windowing of signals amplitude defined by sig_perc .
cfg.source.sig_start= [.2 .2 .2; .2 .2 .2; .2 .2 .2]; %(3signals x Nfreqs);   % start (onset) time of the signals relative to 0=sample(1).
cfg.source.sig_win_type =[3 3 3; 3 3 3; 3 3 3];  %(3signals x Nfreqs) % type of windowing function (1)='Hann' (2)='Gauss' (3)='Triang'  (4)='Blackman' ;
cfg.source.sig_win_rise_time=[.1 .1 .1; .1 .1 .1; .1 .1 .1]; %(3signals x Nfreqs) % rise_time of windowing function. If less than 1/2 sig_duration then window function will have a plateau of the difference in duration between the rise/fall time and the signal duration;

cfg.source.sig_PLV_targets=[.4 .4 .4; .4 .4 .4; .4 .4 .4]; %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
cfg.source.prepost_PLV_targets=[0 0 0; 0 0 0; 0 0 0]; %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
cfg.source.sig_PLI_targets=[.4 .4 .4; .4 .4 .4; .4 .4 .4]; %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
cfg.source.prepost_PLI_targets=[0 0 0; 0 0 0; 0 0 0]; %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
cfg.source.sig_phase_lag=([0 0 0; 0 0 0; 0 0 0]/360)*2*pi;  % phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval relative to first sample --> (phase_lag/360)*2*pi);    cos(phase_lag) = correlation of signal relative to zero-phase onset
cfg.source.prepost_phase_lag=([0 0 0; 0 0 0; 0 0 0]/360)*2*pi;  % phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval relative to first sample --> (phase_lag/360)*2*pi);    cos(phase_lag) = correlation of signal relative to zero-phase onset

cfg.source.phase_amp_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; % sig_contrasts = cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
cfg.source.sig_phase_amp_freq_idx=[0 0 1; 0 0 0; 0 0 0; 0 1 0; 0 0 0; 0 0 0]; %(sig_contrasts x Nfreqs);  % index in the matrix = the carrier frequency and the number at that index = the modulattion frequency --> this example will create a 40-Hz carrier amplitude modulated (index 3) by a 6-Hz modulator phase (number 1 at index).
cfg.source.prepost_phase_amp_freq_idx=[0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0];  %(PLV_contrasts x Nfreqs);     % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
cfg.source.sig_phase_amp_depth_perc=[0 0 75; 0 0 0; 0 0 0; 0 75 0; 0 0 0; 0 0 0]; % amplitude-modulation depth as a percentage of signal's amplitude (sig_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
cfg.source.prepost_phase_amp_depth_perc=[0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]; % depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
cfg.source.sig_phase_amp_depth_perc_range=[0 0 25; 0 0 0; 0 0 0; 0 25 0; 0 0 0; 0 0 0]; % +/- range of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: sig_phase_amp_depth_perc +/- sig_phase_amp_depth_perc_range must be within [0 100]
cfg.source.prepost_phase_amp_depth_perc_range=[0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]; % +/- range devitaion of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: prepost_phase_amp_depth_perc +/- prepost_phase_amp_depth_perc_range must be within [0 100]
%% Run simulation 
run_multiple=0; % (0) run simulation only once and replicate simulation num_iter (1) run simulation multiple times = num_iter
num_iter=1;

% replicate same data to get many trials so that noise phase can be added in other frequency bands and cancel out across more trials for TFR analyses.
[sig_final,sig_wav,prepost_wav,noise_wav,cfg,prepost_win,sig_win] = SimSignals(cfg);
sfg=repmat(sig_final,[1 1 1 num_iter]);
sig_phase=repmat(cfg.sig_phase,[1 1 1 num_iter]);
prepost_phase=repmat(cfg.prepost_phase,[1 1 1 num_iter]);
sig_final=reshape(sfg,[size(sfg,1) size(sfg,2) size(sfg,3)*size(sfg,4)]); sig_final_org=sig_final; % reshape back to single iteration run

% OR
if run_multiple==1
    % Run multiple iterations to yield more trials and so that noise phase can be added in other frequency bands and cancel out across more trials for TFR analyses
    for tt=1:num_iter
        t1=tic;
        fprintf('Iteration #%.f\n',tt);
        [sig_final,sig_wav,prepost_wav,noise_wav,cfg,prepost_win,sig_win] = SimSignals(cfg);
        p_time(tt)=toc(t1);
        sfg(:,:,:,tt)=sig_final;
        cfg_tt(tt)=cfg;
        sig_phase(:,:,:,tt)=cfg_tt(tt).sig_phase; prepost_phase(:,:,:,tt)=cfg_tt(tt).prepost_phase;
    end
    dims=size(sig_phase); sg=permute(reshape(permute(sig_phase,[1 3 2 4]),[dims(1) dims(3) dims(2)*dims(4)]),[1 3 2]);
    dims=size(prepost_phase); pg=permute(reshape(permute(prepost_phase,[1 3 2 4]),[dims(1) dims(3) dims(2)*dims(4)]),[1 3 2]);
end


%% %%% Example #3 Simulation Parameters for simulated real human EEG data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Real ERP data - Loading data from LetterTask LT014
xdata=load('M:\Data\LetterTask_mat\LT014_LetterDetect_MCMV_data.mat','erp_reref','srate','params', 't1_idx2', 't2_idx2');
erp_wvt_chans=[22 97 19 11 40 43]; % 'POOz' 'FCC1h' 'Pz' 'PO10h' 'PO9h' 'P8'
vv=[11 40 19]; % 'PO10h' 'PO9h' 'Pz' 
erp_final=double(squeeze(xdata.erp_reref(:,vv,xdata.t1_idx2(1:90)))); % only 90 trials for letters --> reduce memory overhead
erp_chan_names={'PO10h' 'PO9h' 'Pz' };
% lat=xdata.params.study.lat;
erp_final_org=erp_final(1:size(sig_final,1),:,:);
max_erp=max(max(max(abs(erp_final_org(base_samps,:,:)))));
erp_final=100*(erp_final_org/max_erp);    % normalizing to max for TFR 
%% Study Parameters
clear cfg;
cfg.study.srate=512;
cfg.study.dur=[-2.5 2.5]; % (sec) start and end times for whole trial
cfg.study.lat_sim=[cfg.study.dur(1):1/cfg.study.srate:cfg.study.dur(2)-(1/cfg.study.srate)]; % latency of each trial
cfg.study.num_samps=length(cfg.study.lat_sim);
cfg.study.num_trials=90;
cfg.study.max_perm_plv=9; % integer of cfg.study.num_trials/10 = maximum number of permutations to search for PLV_trials (must be <10 or memory will fail on most computers) .
cfg.study.noise_flag=0;  % Whitening noise to be added to each sources
cfg.study.noise_amp_perc=100; % percent of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
cfg.study.noise_freqs=[.5 100]; % [1 100] (or [1:100] for brain noise) Frequency (Hz) of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
cfg.study.plv_thresh=.05;   % stoppping criterion when search for best PLV/PLI matched to sig_PLV_targets, sig_PLI_targets, etc. (e.g., 0.05).
cfg.study.plot_sim_flag=0;  % plots evoked, trial, and PLV results. Note: This can take time because of filtering.
cfg.study.plot_time_int=[-0.2 .5];   % time interval to plot
cfg.study.plot_freq_int=[3 36]; % frequencies for calculating and plotting PLV/PLI
cfg.study.base_int=[-0.4 0];    % base line interval for plotting
%% Signal Parameters 
cfg.source.sig_freqs(1,:,:) =   [3.75 3.75; 4.5 4.5; 5.5 5.5; 7.25 7.25; 8 8; 8.5 8.5; 9.25 9.25; 10.25 10.25; 11 11; 11.75 11.75; 12.75 12.75; 13.25 13.25; 14.5 14.5; 15.75 15.75; ... % ERS
                                 10 10; 10.5 10.5; 10.75 10.75; 13 13; 19.5 19.5; 21.25 21.25; 22.75 22.75; 15.75 15.75]; % ERD
cfg.source.sig_freqs(2,:,:) =   [3.75 3.75; 4.5 4.5; 5.5 5.5; 7.25 7.25; 8 8; 8.5 8.5; 9.25 9.25; 10.25 10.25; 11 11; 11.75 11.75; 12.75 12.75; 13.25 13.25; 14.5 14.5; 15.75 15.75; ... % ERS
                                 10 10; 10.5 10.5; 10.75 10.75; 13 13; 19.5 19.5; 21.25 21.25; 22.75 22.75; 15.75 15.75]; % ERD
cfg.source.sig_freqs(3,:,:) =   [3.75 3.75; 4.5 4.5; 5.5 5.5; 7.25 7.25; 8 8; 8.5 8.5; 9.25 9.25; 10.25 10.25; 11 11; 11.75 11.75; 12.75 12.75; 13.25 13.25; 14.5 14.5; 15.75 15.75; ... % ERS
                                 10 10; 10.5 10.5; 10.75 10.75; 13 13; 19.5 19.5; 21.25 21.25; 22.75 22.75; 15.75 15.75]; % ERD

cfg.source.sig_amp_perc =       [   .30 .24 .32 .23 .27 .35 .33 .32 .24 .15 .15 .15 .32 .31 ...   
                                    .06 .06 .26 .26 .03 .05 .02 .02; ...
                                    .30 .24 .32 .23 .27 .25 .23 .32 .24 .15 .15 .15 .32 .31 ...
                                    .08 .13 .33 .32 .03 .05 .02 .02; ...
                                    .11 .11 .11 .11 .11 .11 .11 .11 .11 .11 .11 .11 .11 .11 ...
                                    .08 .12 .32 .12 .10 .10 .04 .06; ...
                                    ]*100; 
cfg.source.sig_amp_perc(1,1:4)=cfg.source.sig_amp_perc(1,1:4)-.15;
cfg.source.sig_amp_perc(2,1:4)=cfg.source.sig_amp_perc(2,1:4)-.15;
 
cfg.source.prepost_amp_perc =   [   .08 .11 .07 .13 .15 .07 .04 .07 .04 .06 .02 .04 .09 .09 ...
                                    .18 .18 .38 .54 .08 .17 .06 .17; ...
                                    .08 .11 .07 .03 .05 .07 .04 .07 .04 .06 .02 .04 .09 .08 ...
                                    .33 .39 .65 .68 .08 .17 .06 .17; ...
                                    .12 .12 .12 .12 .12 .12 .12 .12 .12 .12 .12 .12 .12 .12 ...
                                    .10 .13 .32 .13 .16 .18 .06 .03; ...
                                    ]*100; 

cfg.source.sig_amp_perc_std = repmat(5,3,size(cfg.source.sig_freqs,2));   % 3sigs x Nfreqs % standard deviation (as a percent value) of the signal amplitude (sig_amp_perc) across trials
cfg.source.prepost_amp_perc_std = repmat(5,3,size(cfg.source.sig_freqs,2));   % 3sigs x Nfreqs % standard deviation (as a percent value) of the prepost amplitude (sig_amp_perc) across trials
                                
% reducing overall ERS/ERD
cfg.source.sig_amp_perc(:,1:14)=cfg.source.sig_amp_perc(:,1:14)*.6; 
cfg.source.prepost_amp_perc(:,1:14)=cfg.source.prepost_amp_perc(:,1:14)*1; 
cfg.source.sig_amp_perc(:,15:end)=cfg.source.sig_amp_perc(:,15:end)*1; 
cfg.source.prepost_amp_perc(:,15:end)=cfg.source.prepost_amp_perc(:,15:end)*.8; 
                                
                                
cfg.source.sig_evoked_perc =        [repmat(.4,3,14) repmat(.1,3,8)]*100;
cfg.source.prepost_evoked_perc =    [repmat(.1,3,14) repmat(.1,3,8)]*100;

cfg.source.sig_start=           [repmat(.075,3,14)  repmat(.1,3,8)]; 
cfg.source.sig_durs=            [repmat(.4,3,1) repmat(.15,3,13)   repmat(.5,3,8)];
cfg.source.sig_win_rise_time=   [repmat(.015,3,4) repmat(.025,3,10) repmat(.1,3,8)];
cfg.source.sig_win_type =       repmat(3,3,22);

cfg.source.sig_PLV_targets=     [repmat(0.5,3,14) repmat(.05,3,8)];
cfg.source.prepost_PLV_targets= [repmat(.5,3,14) repmat(.05,3,8)];
cfg.source.sig_PLI_targets=     [repmat(.5,3,14) repmat(.05,3,8)];
cfg.source.prepost_PLI_targets= [repmat(.5,3,14) repmat(.05,3,8)];

cfg.source.sig_phase_lag =   (  repmat(0,3,22)/360)*2*pi;  % phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval relative to first sample --> (phase_lag/360)*2*pi);    cos(phase_lag) = correlation of signal relative to zero-phase onset
cfg.source.prepost_phase_lag=(  repmat(0,3,22)/360)*2*pi;  % phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval relative to first sample --> (phase_lag/360)*2*pi);    cos(phase_lag) = correlation of signal relative to zero-phase onset

cfg.source.sig_phase_lag(1:2,1:22)=-pi;
cfg.source.sig_start(3,19:22)=.3; % adjusting source 3 beta activity
cfg.source.sig_durs(3,19:22)=.4; % adjusting source 3 beta activity
cfg.source.sig_win_rise_time(3,19:22)=.1; % adjusting source 3 beta activity
cfg.source.sig_durs(3,1:18)=.9; % adjusting source 3 theta & alpha activity
cfg.source.sig_evoked_perc(3,1:22)=0; 
cfg.source.prepost_evoked_perc(3,1:22)=0; 

% matching PLVs top real PLVs
% 1-2 early broad ERS
xf=[1:4 7 8 13 14];
cfg.source.sig_PLV_targets(1,xf)=.6; cfg.source.prepost_PLV_targets(1,xf)=.2;
cfg.source.sig_PLI_targets(1,xf)=.6; cfg.source.prepost_PLI_targets(1,xf)=.2;
cfg.source.sig_PLV_targets(2:3,xf)=.4; cfg.source.prepost_PLV_targets(2:3,xf)=.35;
cfg.source.sig_PLI_targets(2:3,xf)=.4; cfg.source.prepost_PLI_targets(2:3,xf)=.35;

% 10-Hz ERS for 1-2, 1-3, & 2-3 
cfg.source.sig_start(1,17:18)=.15;
xf=15:18;
cfg.source.sig_start(1,15:16)=.15; cfg.source.sig_start(2,15:16)=.1; cfg.source.sig_start(3,15:16)=.15; 
cfg.source.sig_durs(1:3,xf)=.5; 
cfg.source.sig_win_rise_time(1:3,xf)=.05; 
cfg.source.sig_PLV_targets(1,xf)=.5; cfg.source.prepost_PLV_targets(1,xf)=.2;
cfg.source.sig_PLI_targets(1,xf)=-.5; cfg.source.prepost_PLI_targets(1,xf)=.2;
cfg.source.sig_PLV_targets(2,xf)=.5; cfg.source.prepost_PLV_targets(2,xf)=.1;
cfg.source.sig_PLI_targets(2,xf)=.5; cfg.source.prepost_PLI_targets(2,xf)=.1;
cfg.source.sig_PLV_targets(3,xf)=.5; cfg.source.prepost_PLV_targets(3,xf)=.1;
cfg.source.sig_PLI_targets(3,xf)=.5; cfg.source.prepost_PLI_targets(3,xf)=.1;


% Phase-Amplitude Couplings - all set to zero
cfg.source.phase_amp_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; % sig_contrasts = cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
cfg.source.sig_phase_amp_freq_idx=zeros(size(cfg.source.phase_amp_contrasts,1),size(cfg.source.sig_freqs,2)); %(sig_contrasts x Nfreqs);  % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
cfg.source.prepost_phase_amp_freq_idx=zeros(size(cfg.source.phase_amp_contrasts,1),size(cfg.source.sig_freqs,2));  %(PLV_contrasts x Nfreqs);     % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
cfg.source.sig_phase_amp_depth_perc=zeros(size(cfg.source.phase_amp_contrasts,1),size(cfg.source.sig_freqs,2)); % amplitude-modulation depth as a percentage of signal's amplitude (sig_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
cfg.source.prepost_phase_amp_depth_perc=zeros(size(cfg.source.phase_amp_contrasts,1),size(cfg.source.sig_freqs,2)); % depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
cfg.source.sig_phase_amp_depth_perc_range=zeros(size(cfg.source.phase_amp_contrasts,1),size(cfg.source.sig_freqs,2)); % +/- range of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: sig_phase_amp_depth_perc +/- sig_phase_amp_depth_perc_range must be within [0 100]
cfg.source.prepost_phase_amp_depth_perc_range=zeros(size(cfg.source.phase_amp_contrasts,1),size(cfg.source.sig_freqs,2)); % +/- range devitaion of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: prepost_phase_amp_depth_perc +/- prepost_phase_amp_depth_perc_range must be within [0 100]
%% Run simulation 
run_multiple=0; % (0) run simulation only once and replicate simulation num_iter (1) run simulation multiple times = num_iter
num_iter=1;

% replicate same data to get many trials so that noise phase can be added in other frequency bands and cancel out across more trials for TFR analyses.
[sig_final,sig_wav,prepost_wav,noise_wav,cfg,prepost_win,sig_win] = SimSignals(cfg);
sfg=repmat(sig_final,[1 1 1 num_iter]);
sig_phase=repmat(cfg.sig_phase,[1 1 1 num_iter]);
prepost_phase=repmat(cfg.prepost_phase,[1 1 1 num_iter]);
sig_final=reshape(sfg,[size(sfg,1) size(sfg,2) size(sfg,3)*size(sfg,4)]); sig_final_org=sig_final; % reshape back to single iteration run

% OR
if run_multiple==1
    % Run multiple iterations to yield more trials and so that noise phase can be added in other frequency bands and cancel out across more trials for TFR analyses
    for tt=1:num_iter
        t1=tic;
        fprintf('Iteration #%.f\n',tt);
        [sig_final,sig_wav,prepost_wav,noise_wav,cfg,prepost_win,sig_win] = SimSignals(cfg);
        p_time(tt)=toc(t1);
        sfg(:,:,:,tt)=sig_final;
        cfg_tt(tt)=cfg;
        sig_phase(:,:,:,tt)=cfg_tt(tt).sig_phase; prepost_phase(:,:,:,tt)=cfg_tt(tt).prepost_phase;
    end
    dims=size(sig_phase); sg=permute(reshape(permute(sig_phase,[1 3 2 4]),[dims(1) dims(3) dims(2)*dims(4)]),[1 3 2]);
    dims=size(prepost_phase); pg=permute(reshape(permute(prepost_phase,[1 3 2 4]),[dims(1) dims(3) dims(2)*dims(4)]),[1 3 2]);
end


%% Adding gain factor to convert sig_final to microVolts using base interval.
bs=floor(cfg.study.base_int*cfg.study.srate)-(round(cfg.study.lat_sim(1)*cfg.study.srate)); base_samps=bs(1):bs(2);
max_erp=max(max(rms(erp_final(base_samps,:,:))));
max_sig=max(max(rms(sig_final_org(base_samps,:,:))));
sig_gain=max_erp/max_sig;
sig_final=sig_final_org*sig_gain;
% v=2; sig_final(:,v,:)=sig_final(:,v,:)*-1;
figure(4); clf;
for v=1:3 
    subplot(3,1,v); cla; hold on; plot(lat,squeeze(sig_final(:,v,:)),'b'); plot(lat,squeeze(erp_final(:,v,:)),'color',[1 1 1]*.6);  plot(lat,squeeze(nanmean(sig_final(:,v,:),3)),'r','linewidth',3); plot(lat,squeeze(nanmean(erp_final(:,v,:),3)),'k','linewidth',3); 
    axis([-0.4 1 -100 100]);
end
%% plot ERP by trial
% sig_final
v=1;
c_axis=[-30 30]; act_samps=1306:1408;
[~,xi]=sort(range(squeeze(sig_final(act_samps,v,:))));
xs=sig_final;
figure(1); clf; hold on; surf(cfg.study.lat_sim,1:size(xs,3),squeeze(xs(:,v,xi))'); view(0,90); shading interp; axis tight; caxis([-1 1]); colormap(jet); colorbar; axis([-0.5 1 1 size(xs,3)]);
plot3(lat,(squeeze(nanmean(xs(:,v,:),3))/c_axis(2)*(size(xs,3)/2))+(size(xs,3)/2),ones(size(sig_final,1),1)*c_axis(2),'color',[1 1 1]*0,'linewidth',3); colormap('default')
caxis(c_axis); title('Simulated Data');
% erp_final
[~,xi]=sort(range(squeeze(erp_final(act_samps,v,:))));
xs=erp_final;
figure(2); clf; hold on; surf(cfg.study.lat_sim,1:size(xs,3),squeeze(xs(:,v,xi))'); view(0,90); shading interp; axis tight; caxis([-1 1]); colormap(jet); colorbar; axis([-0.5 1 1 size(xs,3)]);
plot3(lat,(squeeze(nanmean(xs(:,v,:),3))/c_axis(2)*(size(xs,3)/2))+(cfg.study.num_trials/2),ones(size(sig_final,1),1)*c_axis(2),'color',[1 1 1]*0,'linewidth',3); colormap('default')
caxis(c_axis); title('Real ERP Data');

%% %%%%%%%%%%%%%%%%%%% Time-Frequency Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adding notched-noise after simulating with no noise to get clearer TFR results
noise_gain=2.5;
fprintf('Adding notched white-noise. This might take some time ...\n');
freqs=cfg.study.noise_freqs;
num_reps=ceil((4/freqs(1))/(range(lat))); % making sure that at least 4 cycles of lowest filter freq will be within the latency interval
if mod(num_reps,2)==0 % even number ad one rep to get odd number so that middle will be saddled by even number of samples
    num_reps=num_reps+1;
end
y=rand(cfg.study.num_samps*num_reps,num_chans,cfg.study.num_trials)-.5;

nyq=cfg.study.srate*.5;
clear fn
for ff1=1:num_freqs
    f1=squeeze(nanmean(unique(cfg.source.sig_freqs(:,ff1,:)))); fstop=[f1*.95 f1*1.05];
    [b1,a1] = butter(3,fstop/nyq,'stop');
    fn(ff1)=  dfilt.df2t(b1,a1);
end
hcas=dfilt.cascade(fn);
%     fvtool(hcas,'Fs',cfg.study.srate);
y2=filter(hcas,y);

xsamps=round(size(y2,1)/2)-round(.5*cfg.study.num_samps);   % getting data in center of filtered data x to avoid windowing effects
if xsamps(1)==0; xsamps=1; end
xsamps=xsamps:xsamps+cfg.study.num_samps; xsamps=xsamps(1:cfg.study.num_samps);
y2=y2/max(max(max(abs(y2))));
noise_wav2=y2(xsamps,:,:)*(noise_gain);
sig_final=(sig_final_org+noise_wav2)*sig_gain;
%% TFR & PLV/PLI parameters
[num_chans,num_freqs,num_minmaxfr]=size(cfg.source.sig_freqs);
lat=cfg.study.lat_sim;
min_max_freq=cfg.study.plot_freq_int;
%% calculate wavelets (total & induced) under signal final
clear wt wt_ind wt_evk;
sig_ind=bsxfun(@minus,sig_final,nanmean(sig_final,3)); % induced by subtracting mean across trials (i.e., evoked response)
wt=nan(85,2560,3,size(sig_final,3));
wt_ind=wt; wt_evk=nan(85,2560,3);
fprintf('Calculating wavelets ...\n')
for v=1:num_chans
    fprintf('Source = %.f\n',v);
    %% Wavelets
    wt_param=[3 60];
    for t=1:size(sig_final,3)
        [wt(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_final(:,v,t)),'morse',cfg.study.srate,'WaveletParameters',wt_param); % total power
        [wt_ind(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_ind(:,v,t)),'morse',cfg.study.srate,'WaveletParameters',wt_param);    %induced power
    end
    [wt_evk(:,:,v),F,coi_wt]=cwt(squeeze(nanmean(sig_final(:,v,:),3)),'morse',cfg.study.srate,'WaveletParameters',wt_param);    % evoked power
end
F2=flipud(F); wt2=flipud(wt); wt2_ind=flipud(wt_ind); wt2_evk=flipud(wt_evk);
ss=find(cfg.study.lat_sim<=cfg.study.base_int(1)); bs(1)=ss(end);
ss=find(cfg.study.lat_sim<=cfg.study.base_int(2)); bs(2)=ss(end);
base_samps=bs(1):bs(2);
wt3=abs(wt2); % converting to real
wt3_ind=abs(wt2_ind); % converting to real
wt3_evk=abs(wt2_evk); % converting to real
% dividing by baseline then multiply 100 to get percent then baseline
wt_based=bsxfun(@rdivide,wt3,nanmean(wt3(:,base_samps,:,:),2))*100;   % percentage
wt_ind_based=bsxfun(@rdivide,wt3_ind,nanmean(wt3_ind(:,base_samps,:,:),2))*100;   % percentage
wt_evk_based=bsxfun(@rdivide,wt3_evk,nanmean(nanmean(wt3(:,base_samps,:),2),3))*100;   % percentage
% baselining
wt_based=bsxfun(@minus,wt_based,nanmean(wt_based(:,base_samps,:,:),2));   % percentage
wt_ind_based=bsxfun(@minus,wt_ind_based,nanmean(wt_ind_based(:,base_samps,:,:),2));   % percentage
wt_evk_based=bsxfun(@minus,wt_evk_based,nanmean(wt_evk_based(:,base_samps,:),2));   % percentage
avg_wt=nanmean(wt_based,4);
avg_wt_ind=nanmean(wt_ind_based,4);
avg_wt_evk=nanmean(wt_evk_based,4);
%% PLV/PLI calculations based on wavelets
sf=find(F2<=min_max_freq(1)); if isempty(sf); sf=1;end
ef=find(F2<=min_max_freq(2)); if isempty(ef); ef=length(F2);end
f_samps=sf(end):ef(end);
phase_data=angle(wt2(f_samps,:,:,:));

F_plv=F2(f_samps);
coi_wt2=coi_wt; coi_wt2(coi_wt>max(F2(f_samps)))=nan; coi_wt2(coi_wt<min(F2(f_samps)))=nan;
clear plv_data pli_data;
chan_contrasts=nchoosek(1:size(sig_final,2),2); surg_flag=0; num_resamps=1;
clear plv_data pli_data dpli_data;
for f=1:size(phase_data,1)
    [PLV]=calc_PLV_ath(squeeze(phase_data(f,:,:,:)),chan_contrasts,surg_flag,num_resamps);
    PLI_win=range(cfg.study.lat_sim)/50; PLI_win_overlap=PLI_win/2;
    [PLI]=calc_PLI_ath(squeeze(phase_data(f,:,:,:)),cfg.study.srate,cfg.study.lat_sim,PLI_win,PLI_win_overlap,chan_contrasts,surg_flag,num_resamps);
    plv_data(f,:,:)=PLV.PLV; pli_data(f,:,:)=PLI.PLI; dpli_data(f,:,:)=PLI.dPLI;
end
pli_lat=PLI.lat;
plv_based=bsxfun(@minus,plv_data,nanmean(plv_data(:,:,base_samps),3));
ss=find(pli_lat<=0);
ss=find(pli_lat<=cfg.study.base_int(1)); bs(1)=ss(end);
ss=find(pli_lat<=cfg.study.base_int(2)); bs(2)=ss(end);
base_samps_pli=bs(1):bs(2);
pli_based=bsxfun(@minus,pli_data,nanmean(pli_data(:,:,base_samps_pli),3));
dpli_based=bsxfun(@minus,dpli_data,nanmean(dpli_data(:,:,base_samps_pli),3));

%% %%%%%%%%%%%%%%%%%%% Plotting Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializing plotting parameters
% num_iter=1500;  % number of iterations to find PLVs
sig_name={'Signal 1' 'Signal 2' 'Signal 3'};
sig_name=erp_chan_names;
ps=find(cfg.study.lat_sim<=cfg.study.plot_time_int(1)); ps=ps(end); pe=find(cfg.study.lat_sim<=cfg.study.plot_time_int(2)); pe=pe(end); plot_samps=ps:pe;

num_clmns=num_chans; num_rows=num_freqs;
min_max=[-100 100]; % time domain scale as percent of baseline    %[-abs(max(max(max(sig_final)))) abs(max(max(max(sig_final))))];
min_max2=[-60 60]; % wavelet color axis scale as percent of baseline
min_max3=ceil([-max(max(max(abs(sig_final(plot_samps,:,:))))) max(max(max(abs(sig_final(plot_samps,:,:)))))])*1.1; % wavelet color axis scale as percent of baseline
plv_caxis=[-.5 .5]; pli_caxis=[-.5 .5]; dpli_caxis=pli_caxis/2; %[-0.25 0.25];
freq_axis=[5 30];
mrk_clr=[0 .5 1; 0 .6 0; 1 0 0];
plv_clr=[.7 0 .9; 1 0 1; 1 .6 0];
xtik=[-.4:.2:1.2];
f_size=10; % font size for axis & title
f_size2=8; % font size for legend
ln_style={'-' '-' '-'};
t0=find(cfg.study.lat_sim<=0); t0=t0(end);
%% figure(997): Signal final waves
figure(997); set(gcf,'color','w'); clf;
ax=subplot_axes(4,num_clmns,.06,.05,0,0,0);
for v=1:num_chans
    %% Time-domain waves
    axes(ax(v)); cla;  hold on; axis on;
    p1=plot(cfg.study.lat_sim,squeeze(sig_final(:,v,:)),'color',[1 1 1]*.6);
    p2=plot(cfg.study.lat_sim,squeeze(nanmean(sig_final(:,v,:),3)),'color',mrk_clr(v,:),'linewidth',2);
    plot([0 0],[min_max3],'k--');
    axis([cfg.study.plot_time_int min_max3]);  set(gca,'XTick',xtik);
    title(sprintf('Final %s',sig_name{v}),'Color',mrk_clr(v,:)); set(gca,'Fontsize',f_size); box on;
    legend([p1(1),p2],{'Trials','Average'},'Location','NorthWest','FontSize',f_size2)
    if v==1; ylabel('Amplitude (%)'); end
end
%% Power wavelets
for v=1:num_chans
    axes(ax(v+3)); cla;  hold on; axis on;
    surf(cfg.study.lat_sim,F2,squeeze(avg_wt(:,:,v))); view(0,90); shading interp; colormap(jet); axis tight;
    plot3(cfg.study.lat_sim,coi_wt,ones(size(coi_wt)),'color',[1 1 1]*.7,'linewidth',2)
    plot3([0 0],[min_max_freq],[min_max2(2) min_max2(2)],'k--');
    axis([cfg.study.plot_time_int freq_axis]); caxis(min_max2); set(gca,'XTick',xtik,'Fontsize',f_size);
    title(sprintf('Total Power %s',sig_name{v}),'Color',mrk_clr(v,:));
    if v==1; ylabel('Frequency (Hz)'); end
    axes(ax(v+6)); cla;  hold on; axis on;
    surf(cfg.study.lat_sim,F2,squeeze(avg_wt_evk(:,:,v))); view(0,90); shading interp; colormap(jet); axis tight;
    plot3(cfg.study.lat_sim,coi_wt,ones(size(coi_wt)),'color',[1 1 1]*.7,'linewidth',2)
    plot3([0 0],[min_max_freq],[min_max2(2) min_max2(2)],'k--');
    axis([cfg.study.plot_time_int freq_axis]); caxis(min_max2); set(gca,'XTick',xtik,'Fontsize',f_size);
    title(sprintf('Evoked Power: %s',sig_name{v}),'Color',mrk_clr(v,:));
    if v==1; ylabel('Frequency (Hz)'); end
    axes(ax(v+9)); cla;  hold on; axis on;
    surf(cfg.study.lat_sim,F2,squeeze(avg_wt_ind(:,:,v))); view(0,90); shading interp; colormap(jet); axis tight;
    plot3(cfg.study.lat_sim,coi_wt,ones(size(coi_wt)),'color',[1 1 1]*.7,'linewidth',2)
    plot3([0 0],[min_max_freq],[min_max2(2) min_max2(2)],'k--');
    axis([cfg.study.plot_time_int freq_axis]); caxis(min_max2); set(gca,'XTick',xtik,'Fontsize',f_size);
    title(sprintf('Induced Power: %s',sig_name{v}),'Color',mrk_clr(v,:));
    if v==1; ylabel('Frequency (Hz)'); end
    xlabel('Time (sec');
end
ax1=axes('Position',[.84 ax(6).Position(2) .1 ax(12).Position(4)]); axis off; hc=colorbar('peer',ax1,'Location','EastOutside');  ax1.Position(3)=.1; ylabel(hc,'Power (% baseline)'); caxis(min_max2); hc.Label.Position=[2 0 0];
ax2=axes('Position',[.84 ax(9).Position(2) .1 ax(12).Position(4)]); axis off; hc=colorbar('peer',ax2,'Location','EastOutside');  ax2.Position(3)=.1; ylabel(hc,'Power (% baseline)'); caxis(min_max2); hc.Label.Position=[2 0 0];
ax3=axes('Position',[.84 ax(12).Position(2) .1 ax(12).Position(4)]); axis off; hc=colorbar('peer',ax3,'Location','EastOutside');  ax3.Position(3)=.1; ylabel(hc,'Power (% baseline)'); caxis(min_max2); hc.Label.Position=[2 0 0];
%% figure(998): PLV & PLI plots
figure(998); clf; set(gcf,'color','w');
ax=subplot_axes(3,num_clmns,.06,.05,0,0,0);
for vx=1:length(chan_contrasts)
    axes(ax(vx)); cla;  hold on; axis on;
    surf(cfg.study.lat_sim,F_plv,squeeze(plv_based(:,vx,:))); view(0,90); shading interp; colormap(jet);
    %         surf(cfg.study.lat_sim,F_plv,squeeze(plv_data(:,vx,:))); view(0,90); shading interp; colormap(jet);
    plot3(cfg.study.lat_sim,coi_wt2,ones(size(coi_wt2)),'color',[1 1 1]*.7,'linewidth',2);
    plot3([0 0],[min_max_freq],[1 1],'k--');
    title(sprintf('PLV %.f vs %.f',chan_contrasts(vx,:)),'Color',plv_clr(vx,:));
    axis([cfg.study.plot_time_int freq_axis]); caxis(plv_caxis); set(gca,'XTick',xtik,'Fontsize',f_size);
    if vx==1; ylabel('Freq (Hz)'); end
    
    axes(ax(vx+3)); cla;  hold on; axis on;
    surf(pli_lat,F_plv,squeeze(pli_based(:,vx,:))); view(0,90); shading interp; colormap(jet);
    %         surf(pli_lat,F_plv,squeeze(pli_data(:,vx,:))); view(0,90); shading interp; colormap(jet);
    plot3(cfg.study.lat_sim,coi_wt2,ones(size(coi_wt2)),'color',[1 1 1]*.7,'linewidth',2);
    plot3([0 0],[min_max_freq],[1 1],'k--');
    title(sprintf('PLI %.f vs %.f',chan_contrasts(vx,:)),'Color',plv_clr(vx,:));
    axis([cfg.study.plot_time_int freq_axis]); caxis(pli_caxis); set(gca,'XTick',xtik,'Fontsize',f_size);
    if vx==1; ylabel('Freq (Hz)'); end
    
    axes(ax(vx+6)); cla;  hold on; axis on;
    surf(pli_lat,F_plv,squeeze(dpli_based(:,vx,:))); view(0,90); shading interp; colormap(jet);
    %          surf(pli_lat,F_plv,squeeze(dpli_data(:,vx,:))-0.5); view(0,90); shading interp; colormap(jet);
    plot3(cfg.study.lat_sim,coi_wt2,ones(size(coi_wt2)),'color',[1 1 1]*.7,'linewidth',2);
    %         surf(cfg.study.lat_sim,Fcoh,squeeze(nanmean(wcoh,3))); view(0,90); shading interp; colormap(jet); axis tight;
    plot3([0 0],[min_max_freq],[1 1],'k--');
    title(sprintf('dPLI %.f vs %.f',chan_contrasts(vx,:)),'Color',plv_clr(vx,:));
    axis([cfg.study.plot_time_int freq_axis]); caxis(dpli_caxis);  set(gca,'XTick',xtik,'Fontsize',f_size);
    xlabel('Time (sec)');
    if vx==1; ylabel('Freq (Hz)'); end
    
    
end
ax1=axes('Position',[.84 ax(3).Position(2) .1 ax(3).Position(4)]); axis off; hc=colorbar('peer',ax1,'Location','EastOutside');  ax1.Position(3)=.1; ylabel(hc,'PLV'); caxis(plv_caxis); hc.Label.Position=[2 0 0];
ax2=axes('Position',[.84 ax(6).Position(2) .1 ax(6).Position(4)]); axis off; hc=colorbar('peer',ax2,'Location','EastOutside');  ax2.Position(3)=.1; ylabel(hc,'PLI'); caxis(pli_caxis); hc.Label.Position=[2 0 0];
ax3=axes('Position',[.84 ax(9).Position(2) .1 ax(9).Position(4)]); axis off; hc=colorbar('peer',ax3,'Location','EastOutside');  ax3.Position(3)=.1; ylabel(hc,'dPLI'); caxis(dpli_caxis); hc.Label.Position=[2 0 0];

