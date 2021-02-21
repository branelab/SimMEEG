function [data_filt]=filter_data(data,srate,freq,filt_type,bp,Norder)
%function [data_out]=filter_data(data,srate,freq,filt_type,bp,Norder);
%
%   data = [samples x channels x trials]
%   srate = sampling rate
%   freq = frequencies for low, high, or bandpass filter (e.g., [30] for 'low' pass filter; [1 30] for 'bandpass' filter)
%   filt_type = type of filter {'fir' 'butter'}   Note: 'fir' filtering
%       should only be done on long duration trials, for epoched data use
%       the 'butter' filter
%   bp = type of filtering {'low' 'high' 'bandpass'}
%   Norder = filter order
%


% WARNING --> need to restore matlab filter directories because Field Trip's 'external' filtfilt.m do not work with my function!.
%set_paths_BraneLab;

%Norder=2;

if length(freq)==2 & freq(2)>srate/2;
    errordlg(sprintf('High Pass Frequency is GREATER than Nyquist Frequency (sample_rate/2)!\nAdjusting High Pass to Nyquist','Filtering Parameter Error'));
    freq(2)=floor(srate/2);
elseif length(freq)==1 & freq(1)>srate/2;
    errordlg(sprintf('Low Pass Frequency is GREATER than Nyquist Frequency (sample_rate/2)!\nAdjusting Low Pass to Nyquist','Filtering Parameter Error'));
    freq(1)=floor(srate/2);
end
    
[num_samp,num_chans,num_trials]=size(data);
cat_flag=0;

if strcmp(filt_type,'fir')==1
Norder=[];

    nyq = srate*0.5;  % Nyquist frequency
    MINFREQ = 0.1/nyq;
    minfac         = 3; %3;    % this many (lo)cutoff-freq cycles in filter (i.e, need to have at least this number of cycles of lowest freq in an epoch to resolve that frequency).
    min_filtorder  = 15;   % minimum filter length
    trans          = .25; % 0.25; % fractional width of transition zones
    locutoff=freq(1);
    hicutoff=freq(2);
    if exist('Norder')==0 | isempty(Norder)==1; Norder = minfac*fix(srate/freq(1));end
    if Norder < min_filtorder; Norder = min_filtorder; end

    if strcmp(bp,'bandpass')==1
        %fprintf(1,'%.2f\n',Norder)
        if size(data,1)<minfac*Norder; 
            lowest_cf=((3*minfac*srate)/size(data,1))+(0.05*((3*minfac*srate)/size(data,1)));
            fprintf(1,'WARNING! Number of samples %.f must be %.f times the Norder = %.f --> lowest centre frequency = %.2f\n',size(data,1),minfac,Norder,lowest_cf);
            data_filt=[]; cat_flag=1;
            num_loops=ceil( (minfac*Norder)/size(data,1))-1;
            fprintf('Concatenating data %.f times in order to filter data.\n',num_loops);  
            data2=data; for lp=1:num_loops;   data2 = cat(1,data2,data);  end     % concatenating data to allow for filtering. 
            data=data2; clear data2
        %    return;
        end
        if hicutoff>nyq; fprintf(1,'High Cut-off Freq (i.e., low-pass freq) %.2f Hz must be less than Nyquist Freq %.2f Hz\n',hicutoff,nyq);end
%        fprintf('Performing %d-point FIR bandpass filtering. Please wait ...\n',Norder);
        f=[MINFREQ (1-trans)*locutoff/nyq locutoff/nyq hicutoff/nyq (1+trans)*hicutoff/nyq 1];
        m=[0       0                      1            1            0                      0];
    elseif strcmp(bp,'low')==1
%        fprintf('Performing %d-point FIR highpass filtering. Please wait ...\n',Norder);
        f=[MINFREQ locutoff*(1-trans)/nyq locutoff/nyq 1];
        m=[   0             0                   1      1];
    elseif strcmp(bp,'high')==1
%        fprintf('Performing %d-point FIR lowpass filtering. Please wait ...\n',Norder);
        f=[MINFREQ hicutoff/nyq hicutoff*(1+trans)/nyq 1];
        m=[     1           1              0                 0];
    end

    filtwts = firls(Norder,f,m); % get FIR filter coefficients
%    keyboard;
    data_filt=zeros(size(data));

    if num_chans>num_trials;
        for t =1:num_trials
            data_filt(:,:,t) = filtfilt(filtwts,1,squeeze(data(:,:,t)));
        end
    else
        for v =1:num_chans
            data_filt(:,v,:) = filtfilt(filtwts,1,squeeze(data(:,v,:)));
        end
    end
    %   hold off;plot(squeeze(data(:,:,1)),'k','linewidth',2);hold on;plot(squeeze(data_filt(:,:,1)),'r','linewidth',2);

end

if strcmp(filt_type,'butter')==1;
    %if size(data,1)<3*Norder; Norder=size(data,1)/4;end
    Norder=2;
    if strcmp(bp,'bandpass')==1
        lowfreq=freq(2);
        highfreq=freq(1);
        if isempty(Norder);
            [Norder, Wn] = buttord(highfreq/(srate/2), lowfreq/(srate/2), 3, 12); 
        end
        %Norder=2;
        [b,a]=butter(Norder,[highfreq/(srate/2) lowfreq/(srate/2)],'bandpass');
        data_filt=filtfilt(b,a,data);
        
%         keyboard;
%         h1=dfilt.df2(b,a);
%         hfvt=fvtool(h1,'FrequencyScale','log');
    end
    if strcmp(bp,'low')==1
         lowfreq=freq;
         %% lowpass filtering
         Norder=2;
        [b,a]=butter(Norder,lowfreq/(srate/2),'low');
        data_filt=filtfilt(b,a,data);
    end
    if strcmp(bp,'high')==1
        highfreq=freq;
        %%% high pass filtering
         Norder=2;
        [b,a]=butter(Norder,highfreq/(srate/2),'high');
        data_filt=filtfilt(b,a,data);
    end
end

if cat_flag==1; 
    data_filt=data_filt(1:num_samp,1:num_chans,1:num_trials);
end

%hold off;plot(squeeze(data(:,:,1)),'k','linewidth',2);hold on;plot(squeeze(data_filt(:,:,1)),'r','linewidth',2);




