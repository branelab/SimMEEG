function [PLI_data]=sm_calc_PLI_ath(phase_data,srate,lat,PLI_win,PLI_overlap,chan_contrasts,surg_flag,num_resamps)
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
%       PLI_surg_mean = mean of PLI surrogate data
%       PLI_surg_std = standard devitiaton of PLI surrogate data
%       PLI_lat = latency (sec) of the PLI samples
%
%   Written by A. Herdman July 9, 2018
%   updated Dec 28, 2020 - to save memory the surrogate outputs are now mean and stdev only, num_resamps are not in output.
%

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
perc_remain=5;
num_numbers=num_cont*num_samps*num_trials*3*num_win;  % 8bytse per number
[mem_flag,free_mem]=memory_check(num_numbers,perc_remain);

%% PLI
if mem_flag==0 %mem_diff>=0
    blk2_size=floor(free_mem/(num_samps*num_trials*8));
    n_chans=[1:blk2_size:size(chan_contrasts,1)];
    fprintf('Calculating PLI in Blocks\n');
    for nl=1:length(n_chans)
        %     fprintf('Calculating PLI in blocks of %.f samples. Block = %.f of %.f\n',num_samps,nl,length(n_chans));
        if nl==length(n_chans)
            pd_diff=squeeze( pd(:, :, chan_contrasts(n_chans(nl):size(chan_contrasts,1), 1), :) - pd(:, :, chan_contrasts(n_chans(nl):size(chan_contrasts,1), 2), :) );
            my_sine = round(sin(pd_diff)*(10^6))/(10^6); % For MATLAB, 0 can be just an approximate.  Make approximates real zeros.
            sign_test = sign(my_sine);
            if length(n_chans(nl):size(chan_contrasts,1))==1    % only 1 contrast
                PLI(:,n_chans(nl):size(chan_contrasts,1)) = squeeze(abs(mean(nanmean(sign_test, 2), 3))); pd_diff=[]; sign_test=[];
                % dPLI
                Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; my_sine=[];
                dPLI(:,n_chans(nl):size(chan_contrasts,1)) = squeeze(mean(mean(Y, 2), 3)); Y=[];
            else
                PLI(:,n_chans(nl):size(chan_contrasts,1)) = squeeze(abs(mean(nanmean(sign_test, 2), 4))); pd_diff=[]; sign_test=[];
                % dPLI
                Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; my_sine=[];
                dPLI(:,n_chans(nl):size(chan_contrasts,1)) = squeeze(mean(mean(Y, 2), 4)); Y=[];
            end
        else
            pd_diff=squeeze( pd(:, :, chan_contrasts(n_chans(nl):n_chans(nl+1)-1, 1), :) - pd(:, :, chan_contrasts(n_chans(nl):n_chans(nl+1)-1, 2), :) );
            my_sine = round(sin(pd_diff)*(10^6))/(10^6); % For MATLAB, 0 can be just an approximate.  Make approximates real zeros.
            sign_test = sign(my_sine);
            PLI(:,n_chans(nl):n_chans(nl+1)-1) = squeeze(abs(mean(nanmean(sign_test, 2), 4))); pd_diff=[]; sign_test=[];
            % dPLI
            Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; my_sine=[];
            dPLI(:,n_chans(nl):n_chans(nl+1)-1) = squeeze(mean(mean(Y, 2), 4)); Y=[];
        end
    end
    
    
    
else
    
    % vectorized version
    pd_diff=squeeze( pd(:, :, chan_contrasts(:, 1), :) - pd(:, :, chan_contrasts(:, 2), :) );
    my_sine = round(sin(pd_diff)*(10^6))/(10^6); % For MATLAB, 0 can be just an approximate.  Make approximates real zeros.
    sign_test = sign(my_sine);
    PLI = squeeze(abs(mean(squeeze(nanmean(sign_test, 2)), 3))); pd_diff=[]; sign_test=[];
    % dPLI
    Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; my_sine=[];
    dPLI = squeeze(mean(mean(Y, 2), 4)); Y=[];
end

%% Surrogates
if surg_flag ==1
    % Memory check
    %     num_numbers=num_cont*num_samps*num_trials*3*num_win  + (num_cont*num_samps*num_trials*3*num_win*num_resamps);  % 8bytse per number
    num_numbers=num_cont*num_samps*num_trials*3*num_win;  % 8bytse per number
    [mem_flag,free_mem]=memory_check(num_numbers,perc_remain);
    if mem_flag==0
        blk2_size=floor(free_mem/(num_samps*num_trials*8));
        n_chans=[1:blk2_size:size(chan_contrasts,1)];
        fprintf('Calculating PLI Surrogates in %.f Blocks\n',length(n_chans));
        for k = 1:num_resamps
            for nl=1:length(n_chans)
                %             fprintf('Calculating PLI_surg in blocks of %.f samples. Block = %.f of %.f. For Surrgoates %.f of %.f.\n',num_samps,nl,length(n_chans),k,num_resamps);
                if nl==length(n_chans)
                    rt1=randperm(num_trials); rt2=randperm(num_trials);
                    pd_diff=squeeze( pd(:, :, chan_contrasts(n_chans(nl):size(chan_contrasts,1), 1), rt1) - pd(:, :, chan_contrasts(n_chans(nl):size(chan_contrasts,1), 2), rt2) );
                    my_sine = round(sin(pd_diff)*(10^6))/(10^6); % For MATLAB, 0 can be just an approximate.  Make approximates real zeros.
                    sign_test = sign(my_sine);
                    if length(n_chans(nl):size(chan_contrasts,1))==1    % only 1 contrast
                        PLI_surg(k,:,n_chans(nl):size(chan_contrasts,1)) = squeeze(abs(mean(nanmean(sign_test, 2), 3))); pd_diff=[]; sign_test=[];
                        % dPLI
                        Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; my_sine=[];
                        dPLI_surg(k,:,n_chans(nl):size(chan_contrasts,1)) = squeeze(mean(mean(Y, 2), 3)); Y=[];
                    else
                        PLI_surg(k,:,n_chans(nl):size(chan_contrasts,1)) = squeeze(abs(mean(nanmean(sign_test, 2), 4))); pd_diff=[]; sign_test=[];
                        % dPLI
                        Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; my_sine=[];
                        dPLI_surg(k,:,n_chans(nl):size(chan_contrasts,1)) = squeeze(mean(mean(Y, 2), 4)); Y=[];
                    end
                else
                    rt1=randperm(num_trials); rt2=randperm(num_trials);
                    pd_diff=squeeze( pd(:, :, chan_contrasts(n_chans(nl):n_chans(nl+1)-1, 1), rt1) - pd(:, :, chan_contrasts(n_chans(nl):n_chans(nl+1)-1, 2), rt2) );
                    my_sine = round(sin(pd_diff)*(10^6))/(10^6); % For MATLAB, 0 can be just an approximate.  Make approximates real zeros.
                    sign_test = sign(my_sine);
                    PLI_surg(k,:,n_chans(nl):n_chans(nl+1)-1) = squeeze(abs(mean(nanmean(sign_test, 2), 4))); pd_diff=[]; sign_test=[];
                    % dPLI
                    Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; my_sine=[];
                    dPLI_surg(k,:,n_chans(nl):n_chans(nl+1)-1) = squeeze(mean(mean(Y, 2), 4)); Y=[];
                end
            end
        end
        
    else
        PLI_surg=nan(num_resamps,num_samps,num_cont); dPLI_surg=PLI_surg;
        for k = 1:num_resamps
            %             fprintf(1,'Calculating %.f of %.f Surrogate PLI data. Please wait ...\n',k,num_resamps);
            % vectorized version
            rt1=randperm(num_trials); rt2=randperm(num_trials);
            pd_diff=squeeze( pd(:, :, chan_contrasts(:, 1), rt1) - pd(:, :, chan_contrasts(:, 2), rt2) );
            my_sine = round(sin(pd_diff)*(10^6))/(10^6); % For MATLAB, 0 can be just an approximate.  Make approximates real zeros.
            sign_test = sign(my_sine);
            PLI_surg(k,:,:) = squeeze(abs(nanmean(nanmean(sign_test, 2), 4))); pd_diff=[]; sign_test=[];
            % dPLI
            Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; my_sine=[];
            dPLI_surg(k,:,:) = squeeze(mean(mean(Y, 2), 4)); Y=[];
        end
    end
else
    PLI_surg = [];
    dPLI_surg = [];
end



%% OUTPUT
PLI_data.PLI=PLI';
PLI_data.dPLI=dPLI';
PLI_data.PLI_surg=permute(PLI_surg,[1 3 2]);
PLI_data.dPLI_surg=permute(dPLI_surg,[1 3 2]);
PLI_data.lat=PLI_lat;



