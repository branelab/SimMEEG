function [PLI_data]=calc_PLI_ath(phase_data,srate,lat,PLI_win,PLI_overlap,chan_contrasts,surg_flag,num_resamps)
% function [PLI_data]=calc_PLI_ath(phase_data,srate,lat,PLI_win,PLI_overlap,chan_contrasts,surg_flag,num_resamps);
%
% INPUT:
%   phase_data = phase data [samples x channels x trials];
%   srate = sample rate
%   lat = latency for phase_data samples
%   PLI_win = (sec) duration of consecutive intervals to calculate PLI then average
%   PLI_overlap = (sec) duration of overlap between consecutive intervals; must not exceed PLI_win.
%   chan_contrasts = channel contrasts to calculate PLI
%   surg_flag = flag to calculate surrogate data by reshuffling one channel's trial order.
%   num_resamps = number of ersamplings to perform when calculating the surrogate data.
%
% OUTPUT:
%   PLI_data. 
%       PLI = [0, 1] Phase-lag index
%       dPLI = 0.5 when no directionality.
%           signal 1 "leads" signal2 = dPLI < 0.5 .
%           signal 1 "lags" signal2  = dPLI > 0.5 .
%       PLI_surg = PLI surrogate data
%       PLI_lat = latency (sec) of the PLI samples
%
%   Written by A. Herdman July 9, 2018
%

% fprintf(1,'Calculating PLI. Please wait ...\n');
[~,~,num_trials]=size(phase_data);   % Technically, the second output is length(chan_contrasts)

% create phase_data = [window_num, window_samp, num_chan, num_trials].
PLI_samp=round(PLI_win*srate); % number of samples
overlap_samp=round(PLI_overlap*srate); % number of samples

%% Setting up PLI intervals for sign test and averaging
clear PLI_samps PLI_samps2
x1=[1:PLI_samp-overlap_samp:size(phase_data,1)];
for t=1:length(x1)
    PLI_samps(t,:)=x1(t):x1(t)+PLI_samp;
end
dims=size(PLI_samps);
PLI_samps2=reshape(PLI_samps,[numel(PLI_samps) 1]);
pad_samps=max(PLI_samps2)-size(phase_data,1);

if pad_samps>0
    pad_time=pad_samps/srate;
%     fprintf('Warning! Unequal number of PLI windows within data.\nLast PLI interval is padded with %.1f ms (%.f samples) of zeros\n',pad_time*1000,pad_samps);
    pd=padarray(phase_data,[pad_samps 0 0],'post');
else
    pd=phase_data;
end
clear phase_data
PLI_lat=lat(x1)+(PLI_win/2); % PLI_lat is at the beginning of lat
pd=reshape(pd(PLI_samps2,:,:),[dims(1) dims(2) size(pd,2) size(pd,3)]);

%% Calculate PLI: across trials, what is the likelihood that the sign
% is biased toward either positive or negative?
[num_samps,num_win,num_chans,num_trials]=size(pd);
num_cont=size(chan_contrasts,1);

% Memory check
num_numbers=num_cont*num_samps*num_trials*3*num_win;  % 8bytse per number
perc_remain=15;
[mem_flag]=memory_check(num_numbers,perc_remain);
% [xmem,sys_mem]=memory;  free_mem=sys_mem.PhysicalMemory.Available; free_mem=free_mem*.75; % available memory in bytes minus 500 MBytes
% if free_mem<0;
%     fprintf('ERROR! Not enough free memory. Please free up memory at least %.f MegaBytes. \nProgram Terminated.\n',free_mem/1e6);
%     return;
% end
% needed_mem=num_cont*num_samps*num_trials*8*3*num_win;  % 8bytse per number
% mem_diff=(needed_mem-free_mem);
if mem_flag==0 %mem_diff>=0
%     allow_cont=mem_diff/(num_samps*num_trials*8);
%     fprintf('ERROR! Not enough free memory.\nPlease free up memory at least %.f MegaBytes. Reduce PLI_overlap Or  Reduce number of contrasts to %.f\nProgram Terminated.\n',mem_diff/1e6,allow_cont);
    return;
else
%     fprintf('Memory Check! Needed memory = %.f MB of %.f MB Free memory (-25%% overhead).\n',needed_mem/1e6,free_mem/1e6);
%     fprintf(1,'Calculating PLI data. Please wait ...\n');
%     PLI_window=single(nan(num_samps,num_cont,num_trials));
%     dPLI_window=PLI_window;
%    for t=1:num_trials  % need to iterate through to reduce memory overload
%         my_sine = round(sin(squeeze( pd(:, :, chan_contrasts(:, 1), t) - pd(:, :, chan_contrasts(:, 2), t) ))*(10^6))/(10^6); % For MATLAB, 0 can be just an approximate.  Make approximates real zeros.
%         sign_test = sign(my_sine);
%         PLI_window(:,:,t) = squeeze(mean(sign_test, 2));
%         %% Calculate dPLI
%         Y = zeros(size(my_sine));
%         Y(my_sine > 0) = 1;
%         Y(my_sine == 0) = .5;
%         dPLI_window(:,:,t) = mean(Y, 2);
%     end
%     PLI = squeeze(abs(mean(PLI_window, 3))); clear PLI_window;
%     dPLI = squeeze(mean(dPLI_window, 3)); clear dPLI_window;
%         
    % vectorized version
    pd_diff=squeeze( pd(:, :, chan_contrasts(:, 1), :) - pd(:, :, chan_contrasts(:, 2), :) );
    my_sine = round(sin(pd_diff)*(10^6))/(10^6); % For MATLAB, 0 can be just an approximate.  Make approximates real zeros.
    sign_test = sign(my_sine);
    PLI = squeeze(abs(mean(squeeze(nanmean(sign_test, 2)), 3))); pd_diff=[]; sign_test=[]; 
    % dPLI
    Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; my_sine=[];
    dPLI = squeeze(mean(mean(Y, 2), 4)); Y=[];
    
    %% Surrogates
    if surg_flag ==1
        PLI_surg=nan(num_resamps,num_samps,num_cont); dPLI_surg=PLI_surg;
        for k = 1:num_resamps
%             fprintf(1,'Calculating %.f of %.f Surrogate PLI data. Please wait ...\n',k,num_resamps);
            % vectorized version
            rt1=randperm(num_trials); rt2=randperm(num_trials);
            pd_diff=squeeze( pd(:, :, chan_contrasts(:, 1), rt1) - pd(:, :, chan_contrasts(:, 2), rt2) );
            my_sine = round(sin(pd_diff)*(10^6))/(10^6); % For MATLAB, 0 can be just an approximate.  Make approximates real zeros.
            sign_test = sign(my_sine);
            PLI_surg(k,:,:) = squeeze(abs(mean(squeeze(nanmean(sign_test, 2)), 3))); pd_diff=[]; sign_test=[];
            % dPLI
            Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; my_sine=[];
            dPLI_surg(k,:,:) = squeeze(mean(mean(Y, 2), 4)); Y=[];
        end
    else
        PLI_surg=[];  dPLI_surg=[];
    end
end
PLI_data.PLI=PLI';
PLI_data.dPLI=dPLI';
PLI_data.PLI_surg=permute(PLI_surg,[1 3 2]);
PLI_data.dPLI_surg=permute(dPLI_surg,[1 3 2]);
PLI_data.lat=PLI_lat;



