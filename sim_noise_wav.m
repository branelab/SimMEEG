function noise_wav = sim_noise_wav(cfg,num_chans)


%% Noise --> Whitening the signals by adding noise to sig_final based on designated type 'noise_flag'.
if cfg.study.synthetic_noise_flag==1   % adding broad-band white-noise to ssp defined by noise_amp_perc
    noise_wav = nan(cfg.study.num_samps,num_chans,cfg.study.num_trials);
    for v=1:num_chans
        noise_wav(:,v,:)=((rand(cfg.study.num_samps,1,cfg.study.num_trials)-0.5)*2);
    end
elseif   cfg.study.synthetic_noise_flag==2    % Narrow-band white noise for prestimulus noise
    fprintf('Filtering white-noise. This might take some time ...\n');
    if cfg.study.synthetic_noise_freqs(1)<=0
        %             f_type='low'; f_method='fir';
        freqs=cfg.study.synthetic_noise_freqs(2);
    elseif cfg.study.synthetic_noise_freqs(2)<=0
        %             f_type='high'; f_method='fir';
        freqs=cfg.study.synthetic_noise_freqs(1);
    else
        %             f_type='bandpass'; f_method='fir';
        freqs=cfg.study.synthetic_noise_freqs;
    end
    
    num_reps=ceil((4/freqs(1))/(range(cfg.study.lat_sim))); % making sure that at least 4 cycles of lowest filter freq will be within the latency interval
    if mod(num_reps,2)==0 % even number ad one rep to get odd number so that middle will be saddled by even number of samples
        num_reps=num_reps+1;
    end
    y=rand(cfg.study.num_samps*num_reps,num_chans,cfg.study.num_trials)-.5;
    nyq=cfg.study.srate/2; % nyquist frequency
    
%     fprintf('Adding narrow-Band Noise ([%.2f %.2f] Hz) = %.f percent\n',cfg.study.synthetic_noise_freqs,cfg.study.synthetic_noise_amp_perc);
    fprintf('Adding narrow-Band Noise ([%.2f %.2f] Hz) = %.f percent\n',cfg.study.synthetic_noise_freqs);
    
    b=fir1(nyq,cfg.study.synthetic_noise_freqs/nyq); fn =  dfilt.df2t(b); %fvtool(fn,'Fs',cfg.study.srate);
    y2=filter(fn,y);
    xsamps=round(size(y2,1)/2)-round(.5*cfg.study.num_samps);   % getting data in center of filtered data x to avoid windowing effects
    if xsamps(1)==0; xsamps=1; end
    xsamps=xsamps:xsamps+cfg.study.num_samps; xsamps=xsamps(1:cfg.study.num_samps);
    noise_wav=y2(xsamps,:,:);
    noise_wav=noise_wav./max(max(max(abs(noise_wav)))); %scaling -1 to 1
    
elseif cfg.study.synthetic_noise_flag==3            % adding broad-band (white) noise with notches at signal frequencies to be added across epochs with ampltidue defined by 'noise_amp_perc'
    fprintf('Adding notched white-noise. This might take some time ...\n');
    if cfg.study.synthetic_noise_freqs(1)<=0
        %             f_type='low'; f_method='fir';
        freqs=cfg.study.synthetic_noise_freqs(2);
    elseif cfg.study.synthetic_noise_freqs(2)<=0
        %             f_type='high'; f_method='fir';
        freqs=cfg.study.synthetic_noise_freqs(1);
    else
        %             f_type='bandpass'; f_method='fir';
        freqs=cfg.study.synthetic_noise_freqs;
    end
    num_reps=ceil((4/freqs(1))/(range(cfg.study.lat_sim))); % making sure that at least 4 cycles of lowest filter freq will be within the latency interval
    if mod(num_reps,2)==0 % even number ad one rep to get odd number so that middle will be saddled by even number of samples
        num_reps=num_reps+1;
    end
    y=rand(cfg.study.num_samps*num_reps,num_chans,cfg.study.num_trials)-.5;
    
    nyq=cfg.study.srate*.5;
    [~,num_freqs,~]=size(cfg.source.sig_freqs);
    
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
    noise_wav=y2(xsamps,:,:);
    noise_wav=noise_wav./max(max(max(abs(noise_wav)))); %scaling -1 to 1
elseif cfg.study.synthetic_noise_flag==4            % adding pink noise
    nyq=cfg.study.srate/2;
    fprintf('Generating Colored noise. This might take some time ...\n');
    % trpiling num_samps to deal with filtering edge effects - only
    % selecting middle num_samps after filtering
    for vv=1:num_chans
        %             cn = dsp.ColoredNoise('Color','Pink','SamplesPerFrame',3*cfg.study.num_samps,'NumChannels',cfg.study.num_trials,'OutputDataType','double');
        % padding noise with 2 trial at beginning and 2 trials at end to reduce filter edge effects
        cn = dsp.ColoredNoise(cfg.study.synthetic_pink_noise_slope,'SamplesPerFrame',(4*cfg.study.num_samps)+(cfg.study.num_trials*cfg.study.num_samps),'NumChannels',1,'OutputDataType','double');
        y = cn();
        b=fir1(nyq,cfg.study.synthetic_noise_freqs/nyq); fn =  dfilt.df2t(b); %fvtool(fn,'Fs',cfg.study.srate);
        y=filter(fn,y);
        y=y/max(max(abs(y)));
        
        % selecting noise after first 2 trials
        xn = y(2*cfg.study.num_samps:end); xn = xn(1:cfg.study.num_trials*cfg.study.num_samps);
        xnt = reshape(xn,cfg.study.num_samps,cfg.study.num_trials);
        noise_wav(:,vv,:)=xnt;
        noise_wav=noise_wav./max(max(max(abs(noise_wav)))); %scaling -1 to 1
    end
else
    fprintf('No noise added.\n');
    noise_wav=nan(size(sig_wav));
end

fprintf('Finished Simulating Noise Waves\n'); 
