function [sig_final,sig_wav,prepost_wav,noise_wav,cfg,prepost_win,sig_win] = SimSignals(cfg)
%function [sig_final,sig_wav,prepost_wav,noise_wav,cfg] = sim_PLV_3signals (cfg);
%
% Function generates signals composed of a Nfreqs sinusoids whose amplitudes
% are being enhanced (ERS) or suppressed (ERD) after a certain position using a hanning
% windowing defined sig_dur and an ampltidue define by sig_amp_perc. Each
% sinusoid will be composed of combined sinusoids:
%       prepost_sig = sinusoid(Nfreqs) for from the pre/post-signal intervals (inverse ahnning window from (100-sig_amp_perc) to zero amplitude and PLV within the sig_interval)
%       sig = sinusoid(Nfreqs) for the signal interval (windowed from zero amp and PLV in the pre/post intervals to sig_amp_perc amplitude within the sig_interval)
%
% Noise of different types can be added to the overall epocch for whitening.
%
% NOTE: Requires parameters for 3 signals with N sinusoids --> You can ONLY simulate 3 signals; no more, no less.
%
% Required inputs:
%   Study parameters
%     cfg.study.
%       srate=cfg.study.srate;
%       lat_sim = cfg.study.lat_sim;    % latency values for study
%       num_samps=cfg.study.num_samps;
%       num_trials=cfg.study.num_trials;    %(best results = 160); Note: must be an integer value of cfg.study.max_perm_plv.
%       max_perm=cfg.study.max_perm_plv;    % maximum number of permutations to search for PLV_trials (for best results set to 8; must be 8<cfg.study.max_perm_plv<10 or memory will fail on most computers); Note: must be able to equally divide  cfg.study.num_trials by cfg.study.max_perm_plv .
%       noise_flag=cfg.study.noise_flag;
%                       (1) broad-band (white) noise to be added across epochs with ampltidue defined by 'prestim_perc'
%                       (2) narrow-band noise power between 'noise_freqs' to be added across epochs with ampltidue defined by 'prestim_perc'
%                       (3) broad-band (white) noise with notches at signal frequencies to be added across epochs with ampltidue defined by 'prestim_perc'
%                       (4) Brown noise to be added across epochs with ampltidue defined by 'prestim_perc'
%       plot_sim_flag=cfg.study.plot_sim_flag;  % plots evoked, trial, and PLV results. Note: This can take time because of filtering.
%       noise_amp_perc = cfg.study.noise_amp_perc; % percent of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
%       noise_freqs = cfg.study.noise_freqs; % [1 100] (or [1:100] for brain noise) Frequency (Hz) of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
%       plv_thresh=cfg.study.plv_thresh;   % stoppping criterion when search for best PLV matched to PLV_targs and prestim_PLV (e.g., 0.05).
%       mvar_flag = cfg.study.mvar_flag; % Not Implmented yet - Placeholder for future versions: signals will be simulated using multi-variate autoregressive (mVAR) 
%
%   Signal Parameters
%     cfg.source.   
%       sig_freqs = cfg.source.sig_freqs(3signals x Nfreqs x [minfreq maxfreq]); % sinusoids will be randomly assigned betwee [minfreq maxfreq] across num_trials;
%       sig_amp_perc = cfg.source.sig_amp_perc(3signals x Nfreqs); % defines amount of amplitude within the signal duration (sig_dur) for the peak of the Hanning window across sig_dur.
%                    = sig_amp_perc = 100 --> ampltidue of signal goes from 0 pre-signal to 100% in sig_dur back to 0 post-signal interval.
%                    = sig_amp_perc = 60 --> ampltidue of signal goes from 0 pre-signal to 100% in sig_dur back to 0% post-signal interval.
%       cfg.source.sig_amp_perc_std = [5 5; 5 5; 5 5];   % 3sigs x Nfreqs % standard devitaion as a percent of amplitude within the signal duration (sig_dur) for the peak of the Hanning window across sig_dur.
%
%       sig_evoked_perc=cfg.source.sig_evoked_perc(3signals x Nfreqs); % based on signal phases so that the evoked (averaged amplitude) signal
%                    will be the sig_evoked_perc of the sig_amp_perc e.g., if sig_amp_perc=0.8 and sig_evoked_per=0.5, then evoked amplitude will be as 0.4.
%       sig_durs=cfg.source.sig_durs(3signals x Nfreqs);    % duration of Hanning windowing of signals amplitude defined by sig_amp_perc .
%       sig_start=cfg.source.sig_start(3signals x Nfreqs);   % start (onset) time of the signals relative to 0=sample(1).
%       sig_PLV_targets=cfg.source.sig_PLV_targets(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
%       sig_PLI_targets=cfg.source.sig_PLI_targets(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
%       sig_win_type = cfg.source.sig_win_type(3signals x Nfreqs) % type of windowing function (1)='Hann' (2)='Gauss' (3)='Triang'  (4)='Blackman' ;
%       sig_win_rise_time(3signals x Nfreqs)  = cfg.source.sig_win_rise_time % rise_time of windowing function. If less than 1/2 sig_duration then window function will have a plateau of the difference in duration between the rise/fall time and the signal duration;
%       sig_phase_lag=cfg.source.sig_phase_lag(3signals x Nfreqs);  % phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval --> (sig_phase_lag/360)*2*pi);    cos(sig_phase_lag) = correlation of signal relative to zero-phase onset
%                   Warning: large phase lag differences among sources can obliterate PLI and PLI_targets will not be found. The larger the PLI targets, the smaller sig_phase_lag differences are needed among sources..
%
%   Pre/Post Parameters
%     cfg.source.   
%       prepost_amp_perc = cfg.source.sig_amp_perc(3signals x Nfreqs); % defines amplitude for pre/post signal interval that declines to zero within the sig_dur interval using an inverse Hanning window.
%                        = 100 --> ampltidue of signal goes from 100% pre-signal to 0% in sig_dur back to 100% post-signal interval.
%                        = 60 --> ampltidue of signal goes from 60% pre-signal to 0% in sig_dur back to 60% post-signal interval.
%               Note: ERS/ERD =  sig_amp_perc - prepost_amp_perc; Thus, (+)values = event-related amplitude supplement(ERAS); (-)values = event-related amplitude depression (ERAD.
%       prepost_evoked_perc=cfg.source.prepost_evoked_perc(3signals x Nfreqs);    % based on signal phases so that the evoked (averaged amplitude) signal
%                    will be the sig_evoked_perc of the sig_amp_perc e.g., if sig_amp_perc=0.8 and sig_evoked_per=0.5, then evoked amplitude will be as 0.4.
%       prepost_PLV_targets =cfg.source.prepost_PLV(PLV_contrasts x Nfreqs); % PLVcontrast(sig1-sig2, sig1-sig3, sig2-sig3) x Nfreqs --> usually set this to 0 for randomized phase in pre/post intervals for non-event lock phase
%           NOTE: ERS and ERD for each PLV contrast is based on the relative PLVs in sig_PLV and prepost_PLV as defined by:
%                   deltaPLV(sig_interval)=sig_PLV-prepost_PLV;  where positive=event-relate synchronization (ERS), negative=event-related desynchronization(ERD).
%       prepost_PLI_targets =cfg.source.prepost_PLI(PLV_contrasts x Nfreqs); % PLVcontrast(sig1-sig2, sig1-sig3, sig2-sig3) x Nfreqs --> usually set this to 0 for randomized phase in pre/post intervals for non-event lock phase
%           NOTE: ERS and ERD for each PLI contrast is based on the relative PLIs in sig_PLI and prepost_PLI as defined by:
%                   deltaPLI(sig_interval)=sig_PLI-prepost_PLI;  where positive=ERS, negative=ERD.
%       prepost_phase_lag=(cfg.source.prepost_phase_lag/360)*2*pi;   % phase-lag (radians; relative to sample(1)) of each 3 signals within the pre/post-signal intervals --> (prepost_phase_lag/360)*2*pi);    cos(prepost_phase_lag) = correlation of signal relative to zero-phase onset
%                   Warning: large phase lag differences among sources can obliterate PLI and PLI_targets will not be found. The larger the PLI targets, the smaller prepost_phase_lag differences are needed among sources..
%
%
% Output:
%   sig_final = simulated final waveforms, summed across frequency and across sig_wav and noise_wav;  [samples x channels x num_trials]
%   sig_wav = simulated signal waves only for each frequency; [samples x channels x num_trials x freqs]
%   prepost_wav = simulated noise waves only for each frequency; [samples x channels x num_trials x freqs]
%   prepost_win = = windowing envelope used to window the prepost_wav;  [samples x channels x num_trials x freqs]
%   sig_win = windowing envelope used to window the sig_wav; [samples x channels x num_trials x freqs]
%   cfg = structure of configuration inputs plus:
%       cfg.signal inputs plus
%             .sig_PLV_trials_est = Found PLVs within the plv_thresh range of the PLV_targets.
%             .sig_PLV_evoked_est = Found PLVs within the plv_thresh range of the sig_evoked_perc.
%             .sig_PLI_trials_est = Found PLVs within the plv_thresh range of the PLIV_targets.
%             .sig_dPLI_trials_est = Found PLVs within the plv_thresh range of the "directed" PLI_targets.
%             .prepost_PLV_trials_est
%             .prepost_PLV_evoked_est
%             .prepost_PLI_trials_est
%             .prepost_dPLI_trials_est
%       .study
%       .sig_phase = starting pahse for sig_wav
%       .prepost_phase = starting pahse for prepost_wav
%
% written by Anthony Hedman (UBC) on July 19, 2018;
%   updated April 16, 2019 by Anthony Herdman
%       - includes PLV and PLI target estimates
%
%  You can do cross frequency phase-phase coupling PLV;
%   e.g., sig_freqs(1,:,:)=[4 4; 12 12]; sig_freqs(2,:,:)=[10 10; 22 22].
%   This will simulate 4-10Hz and 12-22 Hz phase-to-phase couplings
%
% %% Tutorial: input variables for testing
% clear all; close all;
%
%
%
% %% Run script
% [sig_final,sig_wav,prepost_wav,noise_wav,cfg,prepost_win,sig_win] = sim_PLV_3signals(cfg);

sig_final=[]; sig_wav=[]; prepost_wav=[]; noise_wav=[]; prepost_win=[]; sig_win=[]; 

%% Initializing in-program variables
% num_iter=1500;  % number of iterations to find PLVs
pos=round(cfg.source.sig_start*cfg.study.srate)-(cfg.study.lat_sim(1)*cfg.study.srate); % 3sigs x Nfreqs
freq_rand=0:1/cfg.study.num_trials:1-(1/cfg.study.num_trials); % setting up equally distributed  variable for phase randomization based on # epochs when there are a range of frequncies --> minfr~=maxfr
[num_chans,num_freqs,num_minmaxfr]=size(cfg.source.sig_freqs);
% [num_contrasts,~]=size(sig_PLV_targets);
lat=0:1/cfg.study.srate:cfg.study.num_samps/cfg.study.srate; lat=lat(1:cfg.study.num_samps)';
wavefr=nan(num_chans,cfg.study.num_trials);
sig_wav=nan(cfg.study.num_samps,num_chans,cfg.study.num_trials,num_freqs);
prepost_wav=sig_wav;
sig_win=zeros(cfg.study.num_samps,num_chans,cfg.study.num_trials,num_freqs);
prepost_win=sig_win;
noise_wav=zeros(cfg.study.num_samps,num_chans,cfg.study.num_trials);



%% checking validity of inputs
%
if sum(sum(cfg.source.sig_phase_amp_depth_perc+cfg.source.sig_phase_amp_depth_perc_range>100))>0
    warndlg(sprintf('Phase-Amplitude Coupling parameters exceed 100%\nPlese adjust and try again'),'WARNING!')
    return
elseif sum(sum(cfg.source.sig_phase_amp_depth_perc+cfg.source.sig_phase_amp_depth_perc_range<0))>0
    warndlg(sprintf('Phase-Amplitude Coupling parameters < 0%\nPlese adjust and try again'),'WARNING!')
    return
end

%% running simulation
if num_chans>3
    fprintf('Please simulate ONLY up to 3 signals [currently size(maxfr,2)=%.f but should be <=3 ]\n',num_chans)
elseif num_minmaxfr~=2
    fprintf('Make sure that there sig_freqs has (3signals x Nfreqs x [minfreq maxfreq])\n');
    return;
else
    
    %% Simulating Signal
    %          save('temp_phaselag','sig_phase','prepost_phase');
    %              load('temp_phaselag','sig_phase','prepost_phase');   % for temporary testing
    for f=1:size(cfg.source.sig_freqs,2)
        %% Finding PLV phases for each 3sigs
        fprintf('\n %%%% Finding Signal Phases %%%%\n');
        %          [sig_phase(:,:,f),cfg.source.sig_PLV_trials_est(:,f), cfg.source.sig_PLV_evoked_est(:,f)]=bl_find_PLV_phases(sig_PLV_targets(:,f),sig_evoked_perc(:,f),num_trials,plv_thresh,num_iter);
%        try
       [sig_phase(:,:,f),cfg.source.sig_PLV_trials_est(:,f), cfg.source.sig_PLV_evoked_est(:,f), cfg.source.sig_PLI_trials_est(:,f),cfg.source.sig_dPLI_trials_est(:,f)]=bl_find_PLV_phases_v2(cfg.source.sig_PLV_targets(:,f),cfg.source.sig_evoked_perc(:,f)/100,cfg.source.sig_PLI_targets(:,f),cfg.source.sig_phase_lag(:,f),cfg.study.num_trials,cfg.study.plv_thresh,cfg.study.max_perm_plv);
        fprintf('\n %%%% Finding PrePost Phases %%%%\n');
        %           [prepost_phase(:,:,f),cfg.source.prepost_PLV_trials_est(:,f), cfg.source.prepost_PLV_evoked_est(:,f)]=bl_find_PLV_phases(prepost_PLV_targets(:,f),prepost_evoked_perc(:,f),num_trials,plv_thresh,num_iter);
        [prepost_phase(:,:,f),cfg.source.prepost_PLV_trials_est(:,f), cfg.source.prepost_PLV_evoked_est(:,f), cfg.source.prepost_PLI_trials_est(:,f),cfg.source.prepost_dPLI_trials_est(:,f)]=bl_find_PLV_phases_v2(cfg.source.prepost_PLV_targets(:,f),cfg.source.prepost_evoked_perc(:,f)/100,cfg.source.prepost_PLI_targets(:,f),cfg.source.prepost_phase_lag(:,f),cfg.study.num_trials,cfg.study.plv_thresh,cfg.study.max_perm_plv);
        cfg.sig_phase=sig_phase; cfg.prepost_phase=prepost_phase;
%        catch
%            return
%        end
        for v=1:num_chans
            clear hwin hwin3 dur_samps win_samps han_win1
            %% Signal Interval Sine Waves with windowing
            wavefr(v,:) = (freq_rand * (cfg.source.sig_freqs(v,f,2)-cfg.source.sig_freqs(v,f,1))) + cfg.source.sig_freqs(v,f,1);
            fprintf('Source #%.f   Freqs = %.1f - %.1f\n',v, cfg.source.sig_freqs(v,f,:))
            % iterating across frequency range but phases will be eqaul for these groups.
                        freqs2=cfg.source.sig_freqs(v,f,1):cfg.source.sig_freqs(v,f,2); % in 1-Hz freq steps
            %
            for ff=1:length(freqs2)
%                 wavefr(v,t) = freqs2(ff)
                for t=1:cfg.study.num_trials
                    wavefr(v,t) = freqs2(ff);   
                    %                for fx=1:length(freqs2)
                    %                    wavefr=freqs2(fx);
                    %                    % note: using same wave frequency (wavefr) as that for the signal interval in order to simulate ERS & ERD.
                    %                    sig_wav(:,v,t,f,fx) = sin(lat*2*pi*wavefr + (squeeze(sig_phase(v,t,f))));
                    %                    % note: using same wave frequency (wavefr) as that for the signal interval in order to simulate ERS & ERD.
                    %                    prepost_wav(:,v,t,f,fx) = sin(lat*2*pi*wavefr + (squeeze(prepost_phase(v,t,f))));
                    
                    sig_wav(:,v,t,f,ff) = sin(lat*2*pi*wavefr(v,t) + (squeeze(sig_phase(v,t,f))));
                    prepost_wav(:,v,t,f,ff) = sin(lat*2*pi*wavefr(v,t) + (squeeze(prepost_phase(v,t,f))));
                    
                    %                end
                end
            end
            
            sig_wav=nansum(sig_wav,5);
            prepost_wav=nansum(prepost_wav,5);
            
            % Signal windowing
            dur_samps=ceil(cfg.source.sig_durs(v,f)*cfg.study.srate);
            win_samps=2*ceil(cfg.source.sig_win_rise_time(v,f)*cfg.study.srate);
            
            %             switch cfg.source.sig_win_type{v,f}
            %                 case 'Hann'
            %                     hwin3=hanning( win_samps);
            %                 case 'Triang'
            %                     hwin3=triang( win_samps);
            %                 case 'Gauss'
            %                     hwin3=gausswin( win_samps);
            %                 case 'Blackman'
            %                     hwin3=blackman( win_samps);
            %             end
            if cfg.source.sig_win_type(v,f)==1      % Hanning
                hwin3=hanning(win_samps);
            elseif cfg.source.sig_win_type(v,f)==2  %'Gauss'
                hwin3=gausswin(win_samps);
            elseif cfg.source.sig_win_type(v,f)==3  % 'Triang'
                hwin3=triang(win_samps);
            elseif cfg.source.sig_win_type(v,f)==4  %'Blackman'
                hwin3=blackman(win_samps);
            end
            
            if win_samps<=dur_samps
                hwin=ones(dur_samps,1);
                mid_pt=floor(length(hwin3)/2);
                rhwin=hwin3(1:mid_pt);   % rise window
                fhwin=flipud(rhwin); % fall window
                hwin(1:mid_pt)=rhwin;
                hwin(end-length(fhwin)+1:end)=fhwin;
            else
                hwin=hwin3;
            end
            % scaling hwin from 0-1
            hwin=hwin-min(hwin); hwin=hwin/max(hwin);
            
            han_win1=hwin*cfg.source.sig_amp_perc(v,f)/100; % multipying hanning window by sig_amp_perc
            
            % Prepost windowing
            hwin2=hwin;
            %                 hwin2=2*(-(1./(1+hwin2))+1); % this is a speciallized inverted version of hwin
            han_win2=-(hwin2*(cfg.source.prepost_amp_perc(v,f)/100))+(cfg.source.prepost_amp_perc(v,f)/100);  % inverting hanning window and adding prepost_amp_perc
            %                 clf; hold on; plot(han_win1,'r','linewidth',2);  plot(han_win2,'b','linewidth',2);
            %                 plot((han_win1+han_win2),'k','linewidth',2);
            
            %                 hwin4=han_win1+han_win2; % find the slope
            sig_win(pos(v,f):pos(v,f)+size(han_win1,1)-1,v,:,f)=repmat(han_win1,[1 cfg.study.num_trials]); % placing hanning window into epoch at pos
            prepost_win(:,v,:,f)=cfg.source.prepost_amp_perc(v,f)/100; % adding prepost_amp_perc to entire epoch
            prepost_win(pos(v,f):pos(v,f)+size(han_win2,1)-1,v,:,f)=repmat(han_win2,[1 cfg.study.num_trials]);    % placing inverted hanning window into epoch at pos
            
        end
    end
    sig_wav = sig_wav./max(max(max(abs(sig_wav))));
    prepost_wav = prepost_wav./max(max(max(abs(prepost_wav))));
    
    %% Noise --> Whitening the signals by adding noise to sig_final based on designated type 'noise_flag'.
    if cfg.study.noise_flag==1   % adding broad-band white-noise to ssp defined by noise_amp_perc
        for v=1:num_chans
            noise_wav(:,v,:)=((rand(cfg.study.num_samps,1,cfg.study.num_trials)-0.5)*2)*(cfg.study.noise_amp_perc/100);
        end
    elseif   cfg.study.noise_flag==2    % Narrow-band white noise for prestimulus noise
        fprintf('Filtering white-noise. This might take some time ...\n');
        if cfg.study.noise_freqs(1)<=0
            %             f_type='low'; f_method='fir';
            freqs=cfg.study.noise_freqs(2);
        elseif cfg.study.noise_freqs(2)<=0
            %             f_type='high'; f_method='fir';
            freqs=cfg.study.noise_freqs(1);
        else
            %             f_type='bandpass'; f_method='fir';
            freqs=cfg.study.noise_freqs;
        end
        
        num_reps=ceil((4/freqs(1))/(range(lat))); % making sure that at least 4 cycles of lowest filter freq will be within the latency interval
        if mod(num_reps,2)==0 % even number ad one rep to get odd number so that middle will be saddled by even number of samples
            num_reps=num_reps+1;
        end
        y=rand(cfg.study.num_samps*num_reps,num_chans,cfg.study.num_trials)-.5;
        nyq=cfg.study.srate/2; % nyquist frequency
        
        fprintf('Adding narrow-Band Noise ([%.2f %.2f] Hz) = %.f percent\n',cfg.study.noise_freqs,cfg.study.noise_amp_perc);
        
        b=fir1(nyq,cfg.study.noise_freqs/nyq); fn =  dfilt.df2t(b); %fvtool(fn,'Fs',cfg.study.srate);
        y2=filter(fn,y);
        
        xsamps=round(size(y2,1)/2)-round(.5*num_samps);   % getting data in center of filtered data x to avoid windowing effects
        if xsamps(1)==0; xsamps=1; end
        xsamps=xsamps:xsamps+cfg.study.num_samps; xsamps=xsamps(1:cfg.study.num_samps);
        noise_wav=y2(xsamps,:,:)*2*(cfg.study.noise_amp_perc/100);
    elseif cfg.study.noise_flag==3            % adding broad-band (white) noise with notches at signal frequencies to be added across epochs with ampltidue defined by 'noise_amp_perc'
        fprintf('Adding notched white-noise. This might take some time ...\n');
        if cfg.study.noise_freqs(1)<=0
            %             f_type='low'; f_method='fir';
            freqs=cfg.study.noise_freqs(2);
        elseif cfg.study.noise_freqs(2)<=0
            %             f_type='high'; f_method='fir';
            freqs=cfg.study.noise_freqs(1);
        else
            %             f_type='bandpass'; f_method='fir';
            freqs=cfg.study.noise_freqs;
        end
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
        noise_wav=y2(xsamps,:,:)*2*(cfg.study.noise_amp_perc/100);
    elseif cfg.study.noise_flag==4            % adding pink noise
        nyq=cfg.study.srate/2;
        fprintf('Adding Brownian noise. This might take some time ...\n');
        % trpiling num_samps to deal with filtering edge effects - only
        % selecting middle num_samps after filtering
        for vv=1:num_chans
            cn = dsp.ColoredNoise('Color','brown','SamplesPerFrame',3*cfg.study.num_samps,'NumChannels',cfg.study.num_trials,'OutputDataType','double');
            y = cn();
            b=fir1(nyq,cfg.study.noise_freqs/nyq); fn =  dfilt.df2t(b); %fvtool(fn,'Fs',cfg.study.srate);
            y=filter(fn,y);
            y=y/max(max(abs(y)));
            
            xsamps=round(size(y,1)/2)-round(.5*cfg.study.num_samps);   % getting data in center of filtered data x to avoid windowing effects
            if xsamps(1)==0; xsamps=1; end
            xsamps=xsamps:xsamps+cfg.study.num_samps; xsamps=xsamps(1:cfg.study.num_samps);
            noise_wav(:,vv,:)=y(xsamps,:,:)*(cfg.study.noise_amp_perc/100);
        end
    else
        fprintf('No noise added.\n');
        noise_wav=nan(size(sig_wav));
    end
    
    %% Phase-Amplitude Coupling (see Canolty & Knight 2010 TICS).
    %     sg2=sig_wav;     pg2=prepost_wav;
    sig_win_ap=ones(size(sig_win));
    prepost_win_ap=ones(size(prepost_win));
    rn1=-1+1/size(sig_wav,3):2/size(sig_wav,3):1-1/size(sig_wav,3); % setting up equally distributed randomly assigned windowing amplitude range for phase-amplitude couplings.
    
    for vx=1:size(cfg.source.phase_amp_contrasts,1) % vx = signal contrasts for phase-amp coupling
        va=cfg.source.phase_amp_contrasts(vx,1); % signal index that will have ampltidue modulated
        vp=cfg.source.phase_amp_contrasts(vx,2); % signal index for phase that modulate v1 ampltidue
        for f=1:size(cfg.source.sig_freqs,2)
            
            if cfg.source.sig_phase_amp_freq_idx(vx,f)>0 % apply phase-ampltidue coupling
                rn=rn1(randperm(size(sig_wav,3))); % randomizing window amplitudes across trials
                sig_win_ap(:,va,:,f)=(sig_wav(:,vp,:,cfg.source.sig_phase_amp_freq_idx(vx,f))./max(max(squeeze(sig_wav(:,vp,:,cfg.source.sig_phase_amp_freq_idx(vx,f)))))); % renormalizing back to -1 to 1 for applying phase_amp windowing
                sig_win_ap(:,va,:,f)=squeeze((sig_win_ap(:,va,:,f)+1)/2); % reset range 0 to 1 so that troughs =0 not -1
                rn_range=repmat(rn*(cfg.source.sig_phase_amp_depth_perc_range(vx,f)/100),size(sig_wav,1),1);
                rn_perc=(cfg.source.sig_phase_amp_depth_perc(vx,f)/100)+rn_range;
                rn_diff=1-rn_perc;
                sig_win_ap(:,va,:,f)=(squeeze(sig_win_ap(:,va,:,f)).*rn_perc);
                sig_win_ap(:,va,:,f)=squeeze(sig_win_ap(:,va,:,f))+rn_diff;
            end
            if cfg.source.prepost_phase_amp_freq_idx(vx,f)>0 % apply phase-ampltidue coupling
                rn=rn1(randperm(size(prepost_wav,3))); % randomizing window amplitudes across trials
                prepost_win_ap(:,va,:,f)=(prepost_wav(:,vp,:,cfg.source.prepost_phase_amp_freq_idx(vx,f))./max(max(squeeze(prepost_wav(:,vp,:,cfg.source.prepost_phase_amp_freq_idx(vx,f)))))); % renormalizing back to -1 to 1 for applying phase_amp windowing
                prepost_win_ap(:,va,:,f)=squeeze((prepost_win_ap(:,va,:,f)+1)/2); % reset range 0 to 1 so that troughs =0 not -1
                rn_range=repmat(rn*(cfg.source.prepost_phase_amp_depth_perc_range(vx,f)/100),size(prepost_wav,1),1);
                rn_perc=(cfg.source.prepost_phase_amp_depth_perc(vx,f)/100)+rn_range;
                rn_diff=1-rn_perc;
                prepost_win_ap(:,va,:,f)=(squeeze(prepost_win_ap(:,va,:,f)).*rn_perc);
                prepost_win_ap(:,va,:,f)=squeeze(prepost_win_ap(:,va,:,f))+rn_diff;
            end
        end
    end
    
    %% adding amplitude variation to sig_win and prepost_win based on sig_amp_std prepost_amp_std standard deviation
    sig_rn=permute(repmat(1+(cfg.source.sig_amp_perc_std(v,f)/100)*randn(size(sig_wav,2),size(sig_wav,3),size(sig_wav,4)),[1 1 1 size(sig_wav,1)]),[4 1 2 3]);
    prepost_rn=permute(repmat(1+(cfg.source.prepost_amp_perc_std(v,f)/100)*randn(size(prepost_wav,2),size(prepost_wav,3),size(prepost_wav,4)),[1 1 1 size(prepost_wav,1)]),[4 1 2 3]);
    
    sig_win=sig_win.*sig_rn;
    prepost_win=prepost_win.*prepost_rn;
    clear sig_rn prepost_rn;
    
    %% multiplying phase-amp window with amp window to yield final windowing function
    sig_win=sig_win.*sig_win_ap;
    prepost_win=prepost_win.*prepost_win_ap;
    clear sig_win_ap prepost_win_ap sig_rn prepost_rn;
    
    %% windowing sig_wav and prepost_wav
    sig_wav=sig_wav.*sig_win;
    prepost_wav=prepost_wav.*prepost_win;
    
    
    %% Finalizing Signals by adding everything together
    sig_final=nansum(cat(4,sig_wav/3,prepost_wav/3,noise_wav/3),4)*3; %/sqrt(3);   % sqrt(3) converts it back to -1 to 1 scale
    %     sig_final=nansum(cat(4,sig_wav,prepost_wav,noise_wav),4); %/sqrt(3);   % sqrt(3) converts it back to -1 to 1 scale
end

%% Plotting --> Confirming
if cfg.study.plot_sim_flag==1
    %% %%%%%%%%%%%%%%%%%%% Time-Frequency Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TFR & PLV/PLI parameters
    [num_chans,num_freqs,num_minmaxfr]=size(cfg.source.sig_freqs);
    lat=cfg.study.lat_sim;
    min_max_freq=cfg.study.plot_freq_int;
    %% calculate wavelets (total & induced) under signal final
    clear wt wt_ind wt_evk;
    sig_ind=bsxfun(@minus,sig_final,nanmean(sig_final,3)); % induced by subtracting mean across trials (i.e., evoked response)
    fprintf('Calculating wavelets ...\n')
    for v=1:num_chans
        %% Wavelets - Total Power
%         wt_param=[3 30]; %[3 60];
        TB = 30; % The larger the time-bandwidth parameter, the more spread out the wavelet is in time and narrower the wavelet is in frequency.
        for t=1:size(sig_final,3)
%             [wt(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_final(:,v,t)),'morse',cfg.study.srate,'WaveletParameters',wt_param); % total power
%             [wt_ind(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_ind(:,v,t)),'morse',cfg.study.srate,'WaveletParameters',wt_param);    %induced power
            [wt(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_final(:,v,t)),'morse',cfg.study.srate,'TimeBandwidth',TB); % total power
            [wt_ind(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_ind(:,v,t)),'morse',cfg.study.srate,'TimeBandwidth',TB);    %induced power
        end
%         [wt_evk(:,:,v),F,coi_wt]=cwt(squeeze(nanmean(sig_final(:,v,:),3)),'morse',cfg.study.srate,'WaveletParameters',wt_param);    % evoked power
        [wt_evk(:,:,v),F,coi_wt]=cwt(squeeze(nanmean(sig_final(:,v,:),3)),'morse',cfg.study.srate,'TimeBandwidth',TB);    % evoked power
    end
    F2=flipud(F); wt2=flipud(wt); wt2_ind=flipud(wt_ind); wt2_evk=flipud(wt_evk);
    ss=find(cfg.study.lat_sim<=cfg.study.base_int(1)); bs(1)=ss(end);
    ss=find(cfg.study.lat_sim<=cfg.study.base_int(2)); bs(2)=ss(end);
    base_samps=bs(1):bs(2);
    wt3=abs(wt2); % converting to real 
    wt3_ind=abs(wt2_ind); % converting to real
    wt3_evk=abs(wt2_evk); % converting to real
    %     wt_based=20*log10(bsxfun(@rdivide,wt3,nanmean(wt3(:,base_samps,:),2))); % dB
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
    pos=round(cfg.source.sig_start*cfg.study.srate)-(cfg.study.lat_sim(1)*cfg.study.srate); % 3sigs x Nfreqs
    num_clmns=num_chans; num_rows=num_freqs;
    min_max=[-100 100]; % time domain scale as percent of baseline    %[-abs(max(max(max(sig_final)))) abs(max(max(max(sig_final))))]*100;
    min_max2=[-100 100]; % wavelet color axis scale as percent of baseline
    min_max3=[-max(max(max(abs(sig_final)))) max(max(max(abs(sig_final))))]*110; % wavelet color axis scale as percent of baseline
    plv_caxis=[-.5 .5]; pli_caxis=[-.5 .5]; dpli_caxis=pli_caxis/2; %[-0.25 0.25];
    %     plv_caxis=[0 1]; pli_caxis=[0 1]; dpli_caxis=[-.5 .5]; %[-0.25 0.25];
    
    mrk_clr=[0 .5 1; 0 .6 0; 1 0 0];
    plv_clr=[.7 0 .9; 1 0 1; 1 .6 0];
    xtik=[-.4:.2:1.2];
    f_size=10; % font size for axis & title
    f_size2=8; % font size for legend
    ln_style={'-' '-' '-'};
    t0=find(cfg.study.lat_sim<=0); t0=t0(end);
    
    %% figure(995): Signal & Prepost example waves separate plots
    figure(995); clf; set(gcf,'color','w');
    ax=subplot_axes(num_rows+1,num_clmns,.06,.06,0,0,0);
    v=1; tx=[1 2]; % trial# for example waves
%     min_max4=[-max(abs(squeeze(nansum(sig_win(:,v,tx,:),4)+nansum(prepost_win(:,v,tx,:),4)))) max(abs(squeeze(nansum(sig_win(:,v,tx,:),4)+nansum(prepost_win(:,v,tx,:),4))))]*110;
    min_max4=[-max(max(abs(squeeze(nansum(sig_win(:,v,tx,:),4)+nansum(prepost_win(:,v,tx,:),4))))) max(max(abs(squeeze(nansum(sig_win(:,v,tx,:),4)+nansum(prepost_win(:,v,tx,:),4)))))]*110;
    a=0;
    for f=1:num_freqs
        a=a+1;
        axes(ax(a)); cla; hold on; axis on;
        plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
        p1=plot(cfg.study.lat_sim,squeeze(sig_wav(:,v,tx,f))*100,'color',mrk_clr(v,:),'linewidth',1); %(2).LineWidth=2;
        p2=plot(cfg.study.lat_sim,squeeze(sig_win(:,v,tx,f))*100,'-','color',mrk_clr(v,:),'linewidth',2);
        axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Source %.f Signal(%.1f to %.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:)); box on;
        a=a+1;
        axes(ax(a)); cla; hold on; axis on; box on;
        plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
        p1=plot(cfg.study.lat_sim,squeeze(prepost_wav(:,v,tx,f))*100,'color',mrk_clr(v,:),'linewidth',1); %p1(2).LineWidth=2;
        p2=plot(cfg.study.lat_sim,squeeze(prepost_win(:,v,tx,f))*100,'-','color',mrk_clr(v,:),'linewidth',2);
        axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Source %.f Prepost (%.1f to %.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:)); box on;
        
        a=a+1;
        axes(ax(a)); cla; hold on; axis on; box on;
        plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
        p1=plot(cfg.study.lat_sim,squeeze(sig_wav(:,v,tx,f)+prepost_wav(:,v,tx,f))*100,'color',mrk_clr(v,:),'linewidth',1); %p1(2).LineWidth=2;
        p2=plot(cfg.study.lat_sim,squeeze(sig_win(:,v,tx,f)+prepost_win(:,v,tx,f))*100,'-','color',mrk_clr(v,:),'linewidth',2);
        axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Source %.f (%.1f to %.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:)); box on;
    end
    %% final signal
    a=a+1;
    axes(ax(a)); cla; hold on; axis on;
    plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
    p1=plot(cfg.study.lat_sim,squeeze(nansum(sig_wav(:,v,tx,:),4))*100,'color',mrk_clr(v,:),'linewidth',1); %p1(2).LineWidth=2;
    p2=plot(cfg.study.lat_sim,squeeze(nansum(sig_win(:,v,tx,:),4))*100,'-','color',mrk_clr(v,:),'linewidth',2);
    axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
    title(sprintf('Source %.f Signal Sum',v),'Color',mrk_clr(v,:)); box on;
    a=a+1;
    axes(ax(a)); cla; hold on; axis on; box on;
    plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
    p1=plot(cfg.study.lat_sim,squeeze(nansum(prepost_wav(:,v,tx,:),4))*100,'color',mrk_clr(v,:),'linewidth',1); %p1(2).LineWidth=2;
    p2=plot(cfg.study.lat_sim,squeeze(nansum(prepost_win(:,v,tx,:),4))*100,'-','color',mrk_clr(v,:),'linewidth',2);
    axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
    title(sprintf('Source %.f PrePost Sum',v),'Color',mrk_clr(v,:)); box on;
    
    a=a+1;
    axes(ax(a)); cla; hold on; axis on; box on;
    plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
    p1=plot(cfg.study.lat_sim,squeeze(nansum(sig_wav(:,v,tx,:),4)+nansum(prepost_wav(:,v,tx,:),4))*100,'color',mrk_clr(v,:),'linewidth',1); %p1(2).LineWidth=2;
    p2=plot(cfg.study.lat_sim,squeeze(nansum(sig_win(:,v,tx,:),4)+nansum(prepost_win(:,v,tx,:),4))*100,'-','color',mrk_clr(v,:),'linewidth',2);
    axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
    title(sprintf('Source %.f Sum',v),'Color',mrk_clr(v,:)); box on;
    
    %% figure(996): Signal & Prepost waves overlaid
    figure(996); clf; set(gcf,'color','w');
    ax=subplot_axes(num_rows,num_clmns,.06,.05,0,0,0);
    a=0;
    for f=1:num_freqs
        for v=1:num_chans
            a=a+1;
            axes(ax(a)); cla; hold on; axis on;
            p1=plot(cfg.study.lat_sim,squeeze(prepost_wav(:,v,:,f))*100,'color',[1 1 1]*.6);
            p2=plot(cfg.study.lat_sim,squeeze(sig_wav(:,v,:,f))*100,'color',mrk_clr(v,:));
            p3=plot(cfg.study.lat_sim,squeeze(nanmean(prepost_wav(:,v,:,f),3))*100,'color',[1 1 1]*.4,'linewidth',2);
            p4=plot(cfg.study.lat_sim,squeeze(nanmean(sig_wav(:,v,:,f),3))*100,'color',mrk_clr(v,:)*.75,'linewidth',2);
            plot([0 0],[min_max],'k--');
            axis([cfg.study.plot_time_int min_max]);  set(gca,'XTick',xtik,'Fontsize',f_size);
            title(sprintf('Source %.f (%.1f-%.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:)); box on;
            if v==1; ylabel('Amplitude (%)'); end
            if f==num_freqs; xlabel('Time (sec'); end
            
            if f==1
                legend([p2(1), p4(1), p1(1),  p3(1) ], {'Signal','Avg Signal','Prepost','Avg Prepost'},'Location','NorthWest','Fontsize',f_size2)
            end
        end
    end
    
    %% figure(997): Signal final waves
    figure(997); set(gcf,'color','w'); clf;
    ax=subplot_axes(4,num_clmns,.06,.05,0,0,0);
    for v=1:num_chans
        %% Time-domain waves
        axes(ax(v)); cla;  hold on; axis on;
        p1=plot(cfg.study.lat_sim,squeeze(sig_final(:,v,:))*100,'color',[1 1 1]*.6);
        p2=plot(cfg.study.lat_sim,squeeze(nanmean(sig_final(:,v,:),3))*100,'color',mrk_clr(v,:),'linewidth',2);
        plot([0 0],[min_max3],'k--');
        axis([cfg.study.plot_time_int min_max3]);  set(gca,'XTick',xtik);
        title(sprintf('Source %.f',v),'Color',mrk_clr(v,:)); set(gca,'Fontsize',f_size); box on;
        legend([p1(1),p2],{'Trials','Average'},'Location','NorthWest','FontSize',f_size2)
        if v==1; ylabel('Amplitude (%)'); end
    end
    %% Power wavelets
    for v=1:num_chans
        axes(ax(v+3)); cla;  hold on; axis on;
        surf(cfg.study.lat_sim,F2,squeeze(avg_wt(:,:,v))); view(0,90); shading interp; colormap(jet); axis tight;
        plot3(cfg.study.lat_sim,coi_wt,ones(size(coi_wt)),'color',[1 1 1]*.7,'linewidth',2)
        plot3([0 0],[min_max_freq],[min_max2(2) min_max2(2)],'k--');
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(min_max2); set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Total Power: Source %.f',v),'Color',mrk_clr(v,:));
        if v==1; ylabel('Frequency (Hz)'); end
        axes(ax(v+6)); cla;  hold on; axis on;
        surf(cfg.study.lat_sim,F2,squeeze(avg_wt_evk(:,:,v))); view(0,90); shading interp; colormap(jet); axis tight;
        plot3(cfg.study.lat_sim,coi_wt,ones(size(coi_wt)),'color',[1 1 1]*.7,'linewidth',2)
        plot3([0 0],[min_max_freq],[min_max2(2) min_max2(2)],'k--');
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(min_max2); set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Evoked Power: Source %.f',v),'Color',mrk_clr(v,:));
        if v==1; ylabel('Frequency (Hz)'); end
        axes(ax(v+9)); cla;  hold on; axis on;
        surf(cfg.study.lat_sim,F2,squeeze(avg_wt_ind(:,:,v))); view(0,90); shading interp; colormap(jet); axis tight;
        plot3(cfg.study.lat_sim,coi_wt,ones(size(coi_wt)),'color',[1 1 1]*.7,'linewidth',2)
        plot3([0 0],[min_max_freq],[min_max2(2) min_max2(2)],'k--');
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(min_max2); set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Induced Power: Source %.f',v),'Color',mrk_clr(v,:));
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
        title(sprintf('PLV Source %.f vs %.f',chan_contrasts(vx,:)),'Color',plv_clr(vx,:));
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(plv_caxis); set(gca,'XTick',xtik,'Fontsize',f_size);
        if vx==1; ylabel('Freq (Hz)'); end
        
        axes(ax(vx+3)); cla;  hold on; axis on;
        surf(pli_lat,F_plv,squeeze(pli_based(:,vx,:))); view(0,90); shading interp; colormap(jet);
        %         surf(pli_lat,F_plv,squeeze(pli_data(:,vx,:))); view(0,90); shading interp; colormap(jet);
        plot3(cfg.study.lat_sim,coi_wt2,ones(size(coi_wt2)),'color',[1 1 1]*.7,'linewidth',2);
        plot3([0 0],[min_max_freq],[1 1],'k--');
        title(sprintf('PLI Source %.f vs %.f',chan_contrasts(vx,:)),'Color',plv_clr(vx,:));
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(pli_caxis); set(gca,'XTick',xtik,'Fontsize',f_size);
        if vx==1; ylabel('Freq (Hz)'); end
        
        axes(ax(vx+6)); cla;  hold on; axis on;
        surf(pli_lat,F_plv,squeeze(dpli_based(:,vx,:))); view(0,90); shading interp; colormap(jet);
        %          surf(pli_lat,F_plv,squeeze(dpli_data(:,vx,:))-0.5); view(0,90); shading interp; colormap(jet);
        plot3(cfg.study.lat_sim,coi_wt2,ones(size(coi_wt2)),'color',[1 1 1]*.7,'linewidth',2);
        %         surf(cfg.study.lat_sim,Fcoh,squeeze(nanmean(wcoh,3))); view(0,90); shading interp; colormap(jet); axis tight;
        plot3([0 0],[min_max_freq],[1 1],'k--');
        title(sprintf('dPLI Source %.f vs %.f',chan_contrasts(vx,:)),'Color',plv_clr(vx,:));
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(dpli_caxis);  set(gca,'XTick',xtik,'Fontsize',f_size);
        xlabel('Time (sec)');
        if vx==1; ylabel('Freq (Hz)'); end
    end
    ax1=axes('Position',[.84 ax(3).Position(2) .1 ax(3).Position(4)]); axis off; hc=colorbar('peer',ax1,'Location','EastOutside');  ax1.Position(3)=.1; ylabel(hc,'PLV'); caxis(plv_caxis); hc.Label.Position=[2 0 0];
    ax2=axes('Position',[.84 ax(6).Position(2) .1 ax(6).Position(4)]); axis off; hc=colorbar('peer',ax2,'Location','EastOutside');  ax2.Position(3)=.1; ylabel(hc,'PLI'); caxis(pli_caxis); hc.Label.Position=[2 0 0];
    ax3=axes('Position',[.84 ax(9).Position(2) .1 ax(9).Position(4)]); axis off; hc=colorbar('peer',ax3,'Location','EastOutside');  ax3.Position(3)=.1; ylabel(hc,'dPLI'); caxis(dpli_caxis); hc.Label.Position=[2 0 0];
    
    %% figure(1000) & figure(1001): Polar plot of phases
    bin_wdth=20;
    % get histcounts to set maximum polarhistogram values.
    for f=1:num_freqs
        for v=1:num_chans
            ph=histcounts(squeeze(cfg.prepost_phase(v,:,f)),bin_wdth);
            p1_max(v,f,1)=max(ph);
            ph=histcounts(squeeze(cfg.sig_phase(v,:,f)),bin_wdth);
            p1_max(v,f,2)=max(ph);
            
            pd1=cfg.prepost_phase(chan_contrasts(v,1),:,f)-cfg.prepost_phase(chan_contrasts(v,2),:,f);
            pd2=cfg.sig_phase(chan_contrasts(v,1),:,f)-cfg.sig_phase(chan_contrasts(v,2),:,f);
            ph=histcounts(pd1,bin_wdth);
            pd_max(v,f,1)=max(ph);
            ph=histcounts(pd2,bin_wdth);
            pd_max(v,f,2)=max(ph);
        end
    end
    bin_axis=[0 360 0 max(max(max(p1_max)))+1]; bin_axis2=[0 360 0 max(max(max(pd_max)))+1];
    ab=0; ac=0;
    figure(1000); clf; set(gcf,'color','w');
    ax1=subplot_axes(num_rows,6,0.05,0.05,0,0.05,0); ax1_pos={ax1.Position}; clf;
    figure(1001); clf; set(gcf,'color','w');
    ax2=subplot_axes(num_rows,6,0.05,0.05,0,0,0);  ax2_pos={ax2.Position}; clf;
    
    for f=1:num_freqs
        % Signal phases
        figure(1000);
        for v=1:num_chans
            ab=ab+1;
            axes('Position',ax1_pos{ab}); cla;
            %         polarscatter(squeeze(cfg.prepost_phase(v,:,f)),ones(size(squeeze(cfg.prepost_phase(v,:,f)))),'markerfacecolor',mrk_clr(v,:),'markeredgecolor',mrk_clr(v,:)*.5);
            polarhistogram(squeeze(cfg.prepost_phase(v,:,f)),bin_wdth,'facecolor',mrk_clr(v,:),'edgecolor',mrk_clr(v,:)*.5)
            title(sprintf('PrePost Phases\nSource %.f\n(%.1f-%.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:),'Fontsize',f_size);
            axis(bin_axis);
            axes('Position',ax1_pos{ab+3}); cla;
            %         polarscatter(squeeze(cfg.sig_phase(v,:,f)),ones(size(squeeze(cfg.sig_phase(v,:,f)))),'markerfacecolor',mrk_clr(v,:),'markeredgecolor',mrk_clr(v,:)*.5);
            polarhistogram(squeeze(cfg.sig_phase(v,:,f)),bin_wdth,'facecolor',mrk_clr(v,:),'edgecolor',mrk_clr(v,:)*.5)
            title(sprintf('Signal Phases\nSource %.f\n(%.1f-%.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:),'Fontsize',f_size);
            axis(bin_axis);
        end
        ab=ab+3;
        % Signal phase difference for PLV
        figure(1001);
        for vx=1:num_chans
            ac=ac+1;
            axes('Position',ax2_pos{ac}); cla;
            pd=cfg.prepost_phase(chan_contrasts(vx,1),:,f)-cfg.prepost_phase(chan_contrasts(vx,2),:,f);
            polarhistogram(pd,bin_wdth,'facecolor',plv_clr(vx,:),'edgecolor',plv_clr(vx,:)*.5)
            %         polarscatter(pd,ones(size(pd)),'markerfacecolor',plv_clr(vx,:),'markeredgecolor',plv_clr(vx,:)*.5);
            title(sprintf('PrePost Phase Diff\nSource %.f-%.f \n(%.1f-%.1f Hz)',chan_contrasts(vx,:),squeeze(cfg.source.sig_freqs(v,f,:))),'Color',plv_clr(vx,:),'Fontsize',f_size);
            plv12=abs((sum(exp(1i*(squeeze(pd))),2))')./size(cfg.sig_phase,2);
            my_sine = round(sin(pd)*(10^6))/(10^6); sign_test = sign(my_sine); pli12 = squeeze(abs(mean(squeeze(nanmean(sign_test, 2)), 3)));
            Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; dpli12 = 2*(squeeze(mean(mean(Y, 2), 4))-.5);
            text((250/360)*2*pi,max(max(max(pd_max)))*2.1,sprintf('PLV = %.2f\nPLI = %.2f\ndPLI = %.2f',plv12,pli12,dpli12),'Color',plv_clr(vx,:));
            axis(bin_axis2);
            
            axes('Position',ax2_pos{ac+3}); cla;
            pd=cfg.sig_phase(chan_contrasts(vx,1),:,f)-cfg.sig_phase(chan_contrasts(vx,2),:,f);
            polarhistogram(pd,bin_wdth,'facecolor',plv_clr(vx,:),'edgecolor',plv_clr(vx,:)*.5)
            %       polarscatter(pd,ones(size(pd)),'markerfacecolor',plv_clr(vx,:),'markeredgecolor',plv_clr(vx,:)*.5);
            title(sprintf('Signal Phase Diff\n Source %.f-%.f \n(%.1f-%.1f Hz)',chan_contrasts(vx,:),squeeze(cfg.source.sig_freqs(v,f,:))),'Color',plv_clr(vx,:),'Fontsize',f_size);
            plv12=abs((sum(exp(1i*(squeeze(pd))),2))')./size(cfg.sig_phase,2);
            my_sine = round(sin(pd)*(10^6))/(10^6); sign_test = sign(my_sine); pli12 = squeeze(abs(mean(squeeze(nanmean(sign_test, 2)), 3)));
            Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; dpli12 = 2*(squeeze(mean(mean(Y, 2), 4))-.5);
            text((250/360)*2*pi,max(max(max(pd_max)))*2.1,sprintf('PLV = %.2f\nPLI = %.2f\ndPLI = %.2f',plv12,pli12,dpli12),'Color',plv_clr(vx,:));
            axis(bin_axis2);
        end
        ac=ac+3;
        
    end
end





