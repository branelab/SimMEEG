function [filt_design,filt_data]=filter_data_new(filt_design,data,srate,freqs,filt_type,filt_method,PassbandRipple,StopbandAttenuation,FilterOrder,plot_filter_flag)
%function [filt_design,filt_data]=filter_data_new(filt_design,data,srate,freqs,filt_type,filt_method,PassbandRipple,StopbandAttenuation,FilterOrder,plot_filter_flag)
%   Create a filter and filter data (optional) using Minimum-Order or Known-Order Filter designs using designfilt.m .
%      (see "designfilt.m" at https://www.mathworks.com/help/signal/ref/designfilt.html#bt612s9-1)
%
%   filt_design = filter design and parameters from desigfilt.m; leave empty filt_design=[] to calculate the filter design; else it will be applpited to the data using filt_data=filtfilt(filt_design,data).
%   data = [samples x channels x trials]; % if empty data=[], then only filt_design is returned for applying it to other using "filt_data=filtfilt(filt_design,data);".
%   srate = sampling rate
%   freqs = frequencies for  [lpStartFreq lpStopFreq hpStartFreq hpStopFreq];
%           lowpass filter at 20-Hz with 4-Hz stop-band transiton = [nan nan 20 24].
%           highpass filter at 10-Hz with 2-Hz stop-band transiton = [8 10 nan nan].
%           bandpass filter at 10- to 20-Hz with 2- and 4-Hz stop-band transitons = [8 10 20 24].
%           bandstop filter at 10- to 20-Hz with 2- and 4-Hz stop-band transitons = [10 12 18 20].
%   filt_type = 'lowpassfir' | 'lowpassiir' | 'highpassfir' | 'highpassiir' | 'bandpassfir' | 'bandpassiir' | 'bandstopfir' | 'bandstopiir' | 'differentiatorfir' | 'hilbertfir' | 'arbmagfir'
%               (see designfilt.m for more information);
%   filt_method = if using IIR filter then you can specificy filter method
%           for IIR-filters =  'butter';    'cheby1';    'cheby2';    'ellip';
%           for FIR-filters (minimum-order design)=  'equiripple'; 'kaiserwin';  NOTE: Kaiserwin yields flatter passbands than equiripple.
%           for FIR-filters (Known-order design)=  'equiripple'; 'maxflat'; 'window'; 'ls';
%   PassbandRipple = amount of ripple in passband (decibels)
%   StopbandAttenuation = amount of attenuation in stopband (decibels)
%   FilterOrder = filter order - if empty or not provided then default is to use the PassbandRipple and StopbandAttenuation to determine minimum fitler order.
%                         - if provided then: (see "designfilt.m" at https://www.mathworks.com/help/signal/ref/designfilt.html#bt612s9-1)
%                         IIR   'butter' =  PassbandRipple and StopbandAttenuation are NOT used to design the filter
%                               'cheby1' =  PassbandRipple is only used to design the filter
%                               'cheby2' =  StopbandAttenuation is only used to design the filter
%                               'ellip'  =  PassbandRipple and StopbandAttenuation are used to design the filter
%                         FIR   'equiripple'  =  use PassbandRipple 'PassbandWeight' [0-1] (higher weight=greater ripple attenuation in passband)
%                                               use StopbandAttenuation as 'StopbandWeight' [0-1] (higher weight=greater attenuation in stopband)
%                               'ls'          =  use PassbandRipple 'PassbandWeight' [0-1] (higher weight=greater ripple attenuation in passband)
%                                               use StopbandAttenuation as 'StopbandWeight' [0-1] (higher weight=greater attenuation in stopband)
%                               'window'     =  uses only Hamming window
%   plot_filter_flag = plot fitler design using fvtool.m
%
% OUTPUT
%   filt_design = filter design and parameters from desigfilt.m - can use this to filter other data using "filt_data=filtfilt(filt_design,data);"
%   filt_data = filtered data [samples x channels x trials]
%
%   Created by Anthony Herdman on June 30, 2018
%
%
if nargin==2 % for filtering data only
    srate=[];freqs=[];filt_type=[];filt_method=[];PassbandRipple=[];StopbandAttenuation=[];FilterOrder=[];plot_filter_flag=0;
elseif nargin<9; FilterOrder=[]; plot_filter_flag=0;
elseif nargin<10; plot_filter_flag=0;
end
try
    
    if isempty(filt_design)  % design a filter --> "filt_design"
%         fprintf('Generating a filter design\n');
        switch filt_type
            %% IIR filter designs
            case 'lowpassiir' % low-pass IIR filtering
                if isempty(FilterOrder)
                    filt_design = designfilt('lowpassiir','DesignMethod',filt_method,'StopbandFrequency',freqs(4),'PassbandFrequency',freqs(3),...
                        'PassbandRipple',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,'MatchExactly','passband','SampleRate',srate);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'butter')
                    filt_method='butter';
                    filt_design = designfilt('lowpassiir','DesignMethod',filt_method,'HalfPowerFrequency',freqs(4),'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'cheby1')
                    filt_method='cheby1';
                    filt_design = designfilt('lowpassiir','DesignMethod',filt_method,'PassbandFrequency',freqs(3),...
                        'PassbandRipple',PassbandRipple,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'cheby2')
                    filt_method='cheby2';
                    filt_design = designfilt('lowpassiir','DesignMethod',filt_method,'StopbandFrequency',freqs(4),...
                        'StopbandAttenuation',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'ellip')
                    filt_method='ellip';
                    filt_design = designfilt('lowpassiir','DesignMethod',filt_method,'PassbandFrequency',freqs(3),...
                        'PassbandRipple',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                end
            case 'highpassiir'     % high-pass IIR filtering
                if isempty(FilterOrder)
                    filt_design = designfilt('highpassiir','DesignMethod',filt_method,'StopbandFrequency',freqs(1),'PassbandFrequency',freqs(2),...
                        'PassbandRipple',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,'MatchExactly','passband','SampleRate',srate);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'butter')
                    filt_method='butter';
                    filt_design = designfilt('highpassiir','DesignMethod',filt_method,'HalfPowerFrequency',freqs(1),'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'cheby1')
                    filt_method='cheby1';
                    filt_design = designfilt('highpassiir','DesignMethod',filt_method,'PassbandFrequency',freqs(2),...
                        'PassbandRipple',PassbandRipple,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'cheby2')
                    filt_method='cheby2';
                    filt_design = designfilt('highpassiir','DesignMethod',filt_method,'StopbandFrequency',freqs(1),...
                        'StopbandAttenuation',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'ellip')
                    filt_method='ellip';
                    filt_design = designfilt('highpassiir','DesignMethod',filt_method,'PassbandFrequency',freqs(2),...
                        'PassbandRipple',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                end
            case 'bandpassiir' % band-pass IIR filtering
                if mod(FilterOrder,2)==1; FilterOrder=FilterOrder+1; end % changing odd to even FilterOrder
                if isempty(FilterOrder)
                    filt_design = designfilt('bandpassiir', 'StopbandFrequency1',freqs(1),'PassbandFrequency1',freqs(2), ...
                        'PassbandFrequency2',freqs(3),'StopbandFrequency2',freqs(4),'StopbandAttenuation1',StopbandAttenuation, ...
                        'PassbandRipple',PassbandRipple,'StopbandAttenuation2',StopbandAttenuation, 'DesignMethod',filt_method, ...
                        'MatchExactly','passband', 'SampleRate',srate);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'butter')
                    filt_method='butter';
                    filt_design = designfilt('bandpassiir','DesignMethod',filt_method,'HalfPowerFrequency1',freqs(1),'HalfPowerFrequency2',freqs(4),'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'cheby1')
                    filt_method='cheby1';
                    filt_design = designfilt('bandpassiir','DesignMethod',filt_method,'PassbandFrequency1',freqs(2),'PassbandFrequency2',freqs(3),...
                        'PassbandRipple',PassbandRipple,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'cheby2')
                    filt_method='cheby2';
                    filt_design = designfilt('bandpassiir','DesignMethod',filt_method,'StopbandFrequency1',freqs(1),'StopbandFrequency2',freqs(4),...
                        'StopbandAttenuation',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'ellip')
                    filt_method='ellip';
                    filt_design = designfilt('bandpassiir','DesignMethod',filt_method,'PassbandFrequency1',freqs(2),'PassbandFrequency2',freqs(3),...
                        'PassbandRipple',PassbandRipple,'StopbandAttenuation1',StopbandAttenuation,'StopbandAttenuation2',StopbandAttenuation,...
                        'SampleRate',srate,'FilterOrder',FilterOrder);
                end
            case 'bandstopiir' % band-pass IIR filtering
                if mod(FilterOrder,2)==1; FilterOrder=FilterOrder+1; end % changing odd to even FilterOrder
                if isempty(FilterOrder)
                    filt_design = designfilt('bandstopiir', 'DesignMethod',filt_method,...
                        'PassbandFrequency1',freqs(1), 'PassbandFrequency2',freqs(4),...
                        'StopbandFrequency1',freqs(2),'StopbandFrequency2',freqs(3),...
                        'PassbandRipple1',PassbandRipple, ...
                        'PassbandRipple2',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,  ...
                        'MatchExactly','passband', 'SampleRate',srate);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'butter')
                    filt_method='butter';
                    filt_design = designfilt('bandstopiir','DesignMethod',filt_method,'HalfPowerFrequency1',freqs(1),'HalfPowerFrequency2',freqs(4),'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'cheby1')
                    filt_method='cheby1';
                    filt_design = designfilt('bandstopiir','DesignMethod',filt_method,'PassbandFrequency1',freqs(2),'PassbandFrequency2',freqs(3),...
                        'PassbandRipple',PassbandRipple,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'cheby2')
                    filt_method='cheby2';
                    filt_design = designfilt('bandstopiir','DesignMethod',filt_method,'StopbandFrequency1',freqs(1),'StopbandFrequency2',freqs(4),...
                        'StopbandAttenuation',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'ellip')
                    filt_method='ellip';
                    filt_design = designfilt('bandstopiir','DesignMethod',filt_method,'PassbandFrequency1',freqs(2),'PassbandFrequency2',freqs(3),...
                        'PassbandRipple',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,...
                        'SampleRate',srate,'FilterOrder',FilterOrder);
                end
                
                
                %% FIR filter designs
            case 'lowpassfir' % low-pass IIR filtering
%                 fprintf('Generating a FIR filter can take a few moments. Please be patient ...\n');
                % Minimum-order Designs
                if isempty(FilterOrder) && strcmp(filt_method,'equiripple')
                    filt_method='equiripple';
                    filt_design = designfilt('lowpassfir','DesignMethod',filt_method,'StopbandFrequency',freqs(4),'PassbandFrequency',freqs(3),'PassbandRipple',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,'SampleRate',srate);
                elseif isempty(FilterOrder) && strcmp(filt_method,'kaiserwin')
                    filt_method='kaiserwin';
                    filt_design = designfilt('lowpassfir','DesignMethod',filt_method,'StopbandFrequency',freqs(4),'PassbandFrequency',freqs(3),'PassbandRipple',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,'ScalePassband',true,'SampleRate',srate);
                    % Known-order Designs
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'equiripple')
                    filt_method='equiripple';
                    filt_design = designfilt('lowpassfir','DesignMethod',filt_method,'StopbandFrequency',freqs(4),'PassbandFrequency',freqs(3),'PassbandWeight',PassbandRipple,'StopbandWeight',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'ls')
                    filt_method='ls';
                    filt_design = designfilt('lowpassfir','DesignMethod',filt_method,'StopbandFrequency',freqs(4),'PassbandFrequency',freqs(3),'PassbandWeight',PassbandRipple,'StopbandWeight',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'window')
                    filt_method='window';
                    if mod(FitlerOrder,2)==1; FilterOrder=FilterOrder+1; end % changing odd to even FilterOrder
                    filt_design = designfilt('lowpassfir','DesignMethod',filt_method,'CutoffFrequency',freqs(4),'Window','hamming','ScalePassband',true,'SampleRate',srate,'FilterOrder',FilterOrder);
                end
            case 'highpassfir'     % high-pass IIR filtering
%                 fprintf('Generating a FIR filter can take a few moments. Please be patient ...\n');
                % Minimum-order Designs
                if isempty(FilterOrder) && strcmp(filt_method,'equiripple')
                    filt_method='equiripple';
                    filt_design = designfilt('highpassfir','DesignMethod',filt_method,'StopbandFrequency',freqs(1),'PassbandFrequency',freqs(2),'PassbandRipple',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,'SampleRate',srate);
                elseif isempty(FilterOrder) && strcmp(filt_method,'kaiserwin')
                    filt_method='kaiserwin';
                    filt_design = designfilt('highpassfir','DesignMethod',filt_method,'StopbandFrequency',freqs(1),'PassbandFrequency',freqs(2),'PassbandRipple',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,'ScalePassband',true,'SampleRate',srate);
                    % Known-order Designs
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'equiripple')
                    filt_method='equiripple';
                    filt_design = designfilt('highpassfir','DesignMethod',filt_method,'StopbandFrequency',freqs(1),'PassbandFrequency',freqs(2),'PassbandWeight',PassbandRipple,'StopbandWeight',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'ls')
                    filt_method='ls';
                    filt_design = designfilt('highpassfir','DesignMethod',filt_method,'StopbandFrequency',freqs(1),'PassbandFrequency',freqs(2),'PassbandWeight',PassbandRipple,'StopbandWeight',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'window')
                    filt_method='window';
                    if mod(FitlerOrder,2)==1; FilterOrder=FilterOrder+1; end % changing odd to even FilterOrder
                    filt_design = designfilt('highpassfir','DesignMethod',filt_method,'CutoffFrequency',freqs(1),'Window','hamming','ScalePassband',true,'SampleRate',srate,'FilterOrder',FilterOrder);
                end
            case 'bandpassfir' % band-pass IIR filtering
%                 fprintf('Generating a FIR filter can take a few moments. Please be patient ...\n');
                if mod(FilterOrder,2)==1; FilterOrder=FilterOrder+1; end % changing odd to even FilterOrder
                % Minimum-order Designs
                if isempty(FilterOrder) && strcmp(filt_method,'equiripple')
                    filt_method='equiripple';
                    filt_design = designfilt('bandpassfir','DesignMethod',filt_method,'StopbandFrequency1',freqs(1),'PassbandFrequency1',freqs(2),'StopbandFrequency2',freqs(4),'PassbandFrequency2',freqs(3),'PassbandRipple',PassbandRipple,'StopbandAttenuation1',StopbandAttenuation,'StopbandAttenuation2',StopbandAttenuation,'SampleRate',srate);
                elseif isempty(FilterOrder) && strcmp(filt_method,'kaiserwin')
                    filt_method='kaiserwin';
                    filt_design = designfilt('bandpassfir','DesignMethod',filt_method,'StopbandFrequency1',freqs(1),'PassbandFrequency1',freqs(2),'StopbandFrequency2',freqs(4),'PassbandFrequency2',freqs(3),'PassbandRipple',PassbandRipple,'StopbandAttenuation1',StopbandAttenuation,'StopbandAttenuation2',StopbandAttenuation,'SampleRate',srate);
                    % Known-order Designs
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'equiripple')
                    filt_method='equiripple';
                    filt_design = designfilt('bandpassfir','DesignMethod',filt_method,'StopbandFrequency1',freqs(1),'PassbandFrequency1',freqs(2),'StopbandFrequency2',freqs(4),'PassbandFrequency2',freqs(3),'PassbandWeight',PassbandRipple,'StopbandWeight1',StopbandAttenuation,'StopbandWeight2',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'ls')
                    filt_method='ls';
                    filt_design = designfilt('bandpassfir','DesignMethod',filt_method,'StopbandFrequency1',freqs(1),'PassbandFrequency1',freqs(2),'StopbandFrequency2',freqs(4),'PassbandFrequency2',freqs(3),'PassbandWeight',PassbandRipple,'StopbandWeight1',StopbandAttenuation,'StopbandWeight2',StopbandAttenuation,'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'window')
                    filt_method='window';
                    if mod(FitlerOrder,2)==1; FilterOrder=FilterOrder+1; end % changing odd to even FilterOrder
                    filt_design = designfilt('bandpassfir','DesignMethod',filt_method,'CutoffFrequency1',freqs(1),'CutoffFrequency2',freqs(4),'Window','hamming','ScalePassband',true,'SampleRate',srate,'FilterOrder',FilterOrder);
                end
            case 'bandstopfir' % band-pass IIR filtering
%                 fprintf('Generating a FIR filter can take a few moments. Please be patient ...\n');
                % Minimum-order Designs
                if isempty(FilterOrder) && strcmp(filt_method,'equiripple')
                    filt_method='equiripple';
                    filt_design = designfilt('bandstopfir','DesignMethod',filt_method,...
                        'PassbandFrequency1',freqs(1), 'PassbandFrequency2',freqs(4),'StopbandFrequency1',freqs(2),'StopbandFrequency2',freqs(3),...
                        'PassbandRipple1',PassbandRipple, ...
                        'PassbandRipple2',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,  ...
                        'SampleRate',srate);
                elseif isempty(FilterOrder) && strcmp(filt_method,'kaiserwin')
                    filt_method='kaiserwin';
                    filt_design = designfilt('bandstopfir','DesignMethod',filt_method,...
                        'PassbandFrequency1',freqs(1), 'PassbandFrequency2',freqs(4),'StopbandFrequency1',freqs(2),'StopbandFrequency2',freqs(3),...
                        'PassbandRipple1',PassbandRipple, ...
                        'PassbandRipple2',PassbandRipple,'StopbandAttenuation',StopbandAttenuation,  ...
                        'SampleRate',srate);
                    % Known-order Designs
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'equiripple')
                    filt_method='equiripple';
                    filt_design = designfilt('bandstopfir','DesignMethod',filt_method,...
                        'PassbandFrequency1',freqs(1), 'PassbandFrequency2',freqs(4),'StopbandFrequency1',freqs(2),'StopbandFrequency2',freqs(3),...
                        'PassbandWeight1',PassbandRipple,'PassbandWeight2',PassbandRipple,'StopbandWeight',StopbandAttenuation,...
                        'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'ls')
                    filt_method='ls';
                    filt_design = designfilt('bandstopfir','DesignMethod',filt_method,...
                        'PassbandFrequency1',freqs(1), 'PassbandFrequency2',freqs(4),'StopbandFrequency1',freqs(2),'StopbandFrequency2',freqs(3),...
                        'PassbandWeight1',PassbandRipple,'PassbandWeight2',PassbandRipple,'StopbandWeight',StopbandAttenuation,...
                        'SampleRate',srate,'FilterOrder',FilterOrder);
                elseif ~isempty(FilterOrder) && strcmp(filt_method,'window')
                    filt_method='window';
                    filt_design = designfilt('bandstopfir','DesignMethod',filt_method,'CutoffFrequency1',freqs(2),'CutoffFrequency2',freqs(3),'Window','hamming','ScalePassband',true,'SampleRate',srate,'FilterOrder',FilterOrder);
                end
        end
        filt_data=[];
    elseif ~isempty(filt_design) && isempty(data)   % do nothing because "data" is empty
%         fprintf('Warning! Variable "data" is empty, thus can not apply "filt_design"\n') ;
        filt_data=[];
    else
%         fprintf('Using provided "filt_design" to filter the "data"\n') ;
    end
    
    if ~isempty(filt_design) && ~isempty(data)  % filter the "data" using "filt_design"
%         fprintf('Applying "%s %s %s" fitler to the "data"\n',filt_design.FrequencyResponse, filt_design.ImpulseResponse,filt_design.DesignMethod);
        if ~isempty(data)   % if not empty then apply filt_design to data
            try
                filt_data=filtfilt(filt_design,double(data));
            catch
%                 fprintf('WARNING! data length is too small for filter.\nPadding with zeros at beginning and end of data before filtering\n');
                min_size=3*size(filt_design.Coefficients,2);
                num_reps=ceil(min_size/size(data,1))/2;
%                  fprintf('WARNING! data length is too small for filter.\nPadding with %.f zeros at beginning and end of data before filtering\n',size(data,1)*num_reps);
               y=cat(1,zeros(size(data,1)*num_reps,size(data,2),size(data,3)),data,zeros(size(data,1)*num_reps,size(data,2),size(data,3)));
                y=filtfilt(filt_design,double(y));
                filt_data=y( (1:size(data,1))+(size(data,1)*num_reps),:,:);
            end
        end
        
    end
    
    % plotting "filt_design"
    if plot_filter_flag==1
        fvtool(filt_design,'Analysis','freq');
        legend({sprintf('Magnitude (%s%s %s)',filt_design.FrequencyResponse, filt_design.ImpulseResponse,filt_design.DesignMethod) sprintf('Phase (%s%s %s)',filt_design.FrequencyResponse, filt_design.ImpulseResponse,filt_design.DesignMethod)});
        if strcmp(filt_design.FrequencyResponse,'lowpass')
            axis([0 filt_design.StopbandFrequency*2 -filt_design.StopbandAttenuation*1.5 10]);
        elseif strcmp(filt_design.FrequencyResponse,'highpass')
            axis([filt_design.StopbandFrequency*.5 filt_design.PassbandFrequency*2 -filt_design.StopbandAttenuation*1.5 10]);
        elseif strcmp(filt_design.FrequencyResponse,'bandpass')
            if filt_design.StopbandAttenuation1<=1 && filt_design.StopbandAttenuation2<=1
                axis([filt_design.StopbandFrequency1*.5 filt_design.StopbandFrequency2*1.5 -60 10]);
            elseif filt_design.StopbandAttenuation1>1 || filt_design.StopbandAttenuation2>1
                axis([filt_design.StopbandFrequency1*.5 filt_design.StopbandFrequency2*1.5 -filt_design.StopbandAttenuation1*1.5 10]);
            end
        end
    end
    
catch
    fprintf('WARNING! Filter properties are not set properly.\nPlease make sure the "filt_type" and "filt_method" align\n');
end

return


%% code for testing program
%
% xdata=bl_data.data(3,:)';
% s=13110;
% xdata=bsxfun(@minus,xdata,nanmean(xdata(s:s+srate,:),1));
% fdata=filter_data(xdata,bl_data.srate,freqs(2:3),'butter','bandpass',[]);
% filt_method={'butter'    'cheby1'    'cheby2'    'ellip'};
% filt_method={'equiripple' 'kaiserwin'};
% clear filt_data
% for f=1:length(filt_method)
%     [filt_data(:,f)]=filter_data_new(xdata(:,1),srate,freqs,'bandpassfir',filt_method{f},1,60,[],1);
%
% figure(100+f); clf; hold on;
% plot(squeeze(xdata(:,1)),'k'); plot(squeeze(filt_data(:,f)),'r','linewidth',2); plot(squeeze(fdata),'g','linewidth',2);axis([s s+srate -15 15]);
% end
% figure(999); clf; hold on;
% plot(squeeze(xdata(:,1)),'k'); plot(squeeze(filt_data(:,:)),'linewidth',2); plot(squeeze(fdata),'g','linewidth',2);axis([s s+srate -15 15]);
% legend({'raw' filt_method{:} 'old'});


