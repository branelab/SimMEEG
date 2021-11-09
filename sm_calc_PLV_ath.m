function [PLV_data]=sm_calc_PLV_ath(phase_data,chan_contrasts,surg_flag,num_resamps)
%function [PLV_data]=calc_PLV(phase_data,chan_contrasts,surg_flag,num_resamps);
%
% INPUT:
%   phase_data = instantaneous phases of the filtered-epoched data from hilbert transform [samples x channels x trial]
%   chan_contrasts = channel indices for contrasting to calculate inter-channel PLV
%   surg_flag = (1) calculate surrogate data (0) no surrogate data - default = 0
%   num_resamps = number of times to resample the PLV surrogate data to create a null distribution
%
% OUTPUT:
%   PLV_data.
%       PLV = [contrasts x samps]
%       PLV_surg = mean of PLV surrogates [contrasts x samps]
%
%   Written by A. Herdman July 9, 2018


if nargin <3; surg_flag=0;end
PLV_surg=[];PLV=[];

%%% Calculating PLV data
[num_samp,num_chan,num_trials]=size(phase_data);
PLV=single(zeros(size(chan_contrasts,1),num_samp));

%%% Calculating PLV data
[~,sys_mem]=memory;  free_mem=sys_mem.PhysicalMemory.Available; free_mem=free_mem*.95; % if exceeding 95% of available memory in bytes, then calulating blocks of PLV calculations across contrasts
if free_mem<0
    fprintf('ERROR! Not enough free memory. Please free up memory at least %.f MegaBytes. \nProgram Terminated.\n',free_mem/1e6);
    return;
end
if num_resamps>0
    needed_mem=size(chan_contrasts,1)*num_samp*num_trials*8 + (size(chan_contrasts,1)*num_samp*num_trials*8*num_resamps);  % 8bytse per number
else
    needed_mem=size(chan_contrasts,1)*num_samp*num_trials*8;  % 8bytse per number
end

%% Calculating PLV
if size(chan_contrasts,1)==1
    PLV=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(:,1),:)-phase_data(:,chan_contrasts(:,2),:)))),2))')./num_trials;
elseif size(phase_data,1)==1   % only 1 sample
    PLV=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(:,1),:)-phase_data(:,chan_contrasts(:,2),:)))),2))')./num_trials;
elseif needed_mem>free_mem
    blk2_size=floor(free_mem/(num_samp*num_trials*8));
    n_chans=[1:blk2_size:size(chan_contrasts,1)];
        fprintf('Calculating PLV in Blocks\n');

    for nl=1:length(n_chans)
        %         try
        fprintf('Calculating PLV in blocks of %.f samples. Block = %.f of %.f\n',num_samp,nl,length(n_chans));
        if nl==length(n_chans)
            PLV(n_chans(nl):size(chan_contrasts,1),:)=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(n_chans(nl):size(chan_contrasts,1),1),:)-phase_data(:,chan_contrasts(n_chans(nl):size(chan_contrasts,1),2),:)))),3))')./num_trials;
        else
            PLV(n_chans(nl):n_chans(nl+1)-1,:)=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(n_chans(nl):n_chans(nl+1)-1,1),:)-phase_data(:,chan_contrasts(n_chans(nl):n_chans(nl+1)-1,2),:)))),3))')./num_trials;
        end
        %         catch ME
        %             fprintf('ERROR! %s\n',ME);
        %         end
    end
else
    PLV=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(:,1),:)-phase_data(:,chan_contrasts(:,2),:)))),3))')./num_trials;
end



%%% Calculating surrogate PLV data
if surg_flag ==1
    PLV_surg=single(zeros(num_resamps,size(chan_contrasts,1),num_samp));
    needed_mem=2*size(chan_contrasts,1)*num_samp*num_trials*8;  % 8bytse per number
    [~,sys_mem]=memory;  free_mem=sys_mem.PhysicalMemory.Available; free_mem=free_mem*.95; % available memory in bytes minus 500 MBytes
    if needed_mem>free_mem;    fprintf('Calculating PLI in Blocks\n'); end

    for np=1:num_resamps
        if size(chan_contrasts,1)==1 || size(phase_data,1)==1
            PLV_surg(np,:,:)=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(:,1),:)-phase_data(:,chan_contrasts(:,2),randperm(size(phase_data,3)))))),2))')./num_trials;
        elseif needed_mem>free_mem
            blk2_size=floor(free_mem/(num_samp*num_trials*8));
            n_chans=[1:blk2_size:size(chan_contrasts,1)];
            
            for nl=1:length(n_chans)
                %                 try
%                 fprintf('Calculating PLV_surg in blocks of %.f samples. Block = %.f of %.f\n',num_samp,nl,length(n_chans));
                if nl==length(n_chans)
                    PLV_surg(np,n_chans(nl):size(chan_contrasts,1),:)=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(n_chans(nl):size(chan_contrasts,1),1),:)-phase_data(:,chan_contrasts(n_chans(nl):size(chan_contrasts,1),2),randperm(size(phase_data,3)))))),3))')./num_trials;
                else
                    PLV_surg(np,n_chans(nl):n_chans(nl+1)-1,:)=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(n_chans(nl):n_chans(nl+1)-1,1),:)-phase_data(:,chan_contrasts(n_chans(nl):n_chans(nl+1)-1,2),randperm(size(phase_data,3)))))),3))')./num_trials;
                end
                %                 catch ME
                %                     fprintf('ERROR! %s\n',ME);
                %                 end
            end
        else
            PLV_surg(np,:,:)=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(:,1),:)-phase_data(:,chan_contrasts(:,2),randperm(size(phase_data,3)))))),3))')./num_trials;
        end
        
        %toc
    end
else
    PLV_surg=[];
end
PLV_data.PLV=PLV;
PLV_data.PLV_surg = PLV_surg;

