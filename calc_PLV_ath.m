function [PLV_data]=calc_PLV_ath(phase_data,chan_contrasts,surg_flag,num_resamps)
%function [PLV_data]=calc_PLV(phase_data,chan_contrasts,surg_flag,num_resamps);
%
%   phase_data = instantaneous phases of the filtered-epoched data from hilbert transform [samples x channels x trial]
%   chan_contrasts = channel indices for contrasting to calculate inter-channel PLV
%   surg_flag = (1) calculate surrogate data (0) no surrogate data - default = 0
%   num_resamps = number of times to resample the PLV surrogate data to create a null distribution
%
%keyboard;
if nargin <3; surg_flag=0;end
PLV_surg=[];PLV=[];
%keyboard;
%%% Calculating PLV data
% fprintf(1,'Calculating PLV. Please wait ...\n');
[num_samp,num_chan,num_trials]=size(phase_data);
%phase_data=unwrap(phase_data);
PLV=single(zeros(size(chan_contrasts,1),num_samp));
%phase_data=single(phase_data);
%phase_data=codistributed(single(phase_data));
%tic; parfor s=1:size(phase_data,1); PLV(s,:)=abs((sum(exp(1i*(squeeze(phase_data(s,chan_contrasts(:,1),:)-phase_data(s,chan_contrasts(:,2),:)))),2))')./num_trials; end; toc
%tic; for s=1:size(phase_data,1); PLV(s,:)=abs((sum(exp(1i*(squeeze(phase_data(s,chan_contrasts(:,1),:)-phase_data(s,chan_contrasts(:,2),:)))),2))')./num_trials; end; toc

%%% Calculating PLV data
[~,sys_mem]=memory;  free_mem=sys_mem.PhysicalMemory.Available; free_mem=free_mem*.35; % available memory in bytes minus 500 MBytes
if free_mem<0;
    fprintf('ERROR! Not enough free memory. Please free up memory at least %.f MegaBytes. \nProgram Terminated.\n',free_mem/1e6);
    return;
end
needed_mem=size(chan_contrasts,1)*num_samp*num_trials*8;  % 8bytse per number

if size(chan_contrasts,1)==1;
%    fprintf(1,'Calculating plv for 1 contrast\n');
    PLV=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(:,1),:)-phase_data(:,chan_contrasts(:,2),:)))),2))')./num_trials;
elseif size(phase_data,1)==1;   % only 1 sample
%    fprintf(1,'Calculating plv for 1 sample\n');
    PLV=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(:,1),:)-phase_data(:,chan_contrasts(:,2),:)))),2))')./num_trials;
%elseif size(chan_contrasts,1)*num_samp*num_trials>8e8     % 800 million is an arbitrary limit to reduce memory consumption
%    blk_size=floor(6e8/(size(chan_contrasts,1)*num_trials)); % finding sample bock size that will not exceed the 600 million threshold
%     if blk_size<1;
elseif needed_mem>free_mem                            
%        fprintf('"chan_contrasts" too large to calculate PLV! subdividing into blocks of chan_contrasts.\n');
        %blk2_size=floor(8e8/(num_samp*num_trials));
        blk2_size=floor(free_mem/(num_samp*num_trials*8));
        n_chans=[1:blk2_size:size(chan_contrasts,1)];
        for nl=1:length(n_chans);
            try
%                fprintf('Calculating PLV in blocks of %.f samples. Block = %.f of %.f\n',num_samp,nl,length(n_chans));
                if nl==length(n_chans)
                    PLV(n_chans(nl):size(chan_contrasts,1),:)=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(n_chans(nl):size(chan_contrasts,1),1),:)-phase_data(:,chan_contrasts(n_chans(nl):size(chan_contrasts,1),2),:)))),3))')./num_trials;
                else
                    PLV(n_chans(nl):n_chans(nl+1)-1,:)=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(n_chans(nl):n_chans(nl+1)-1,1),:)-phase_data(:,chan_contrasts(n_chans(nl):n_chans(nl+1)-1,2),:)))),3))')./num_trials;
                end
            catch ME
                keyboard
            end
        end
%     else
%         n_samps=[1:blk_size:size(phase_data,1)];
%         for nl=1:length(n_samps);   % looping through blks
%             try
%                 fprintf('Calculating PLV in blocks of %.f samples. Block = %.f of %.f\n',blk_size,nl,length(n_samps));
%                 if nl==length(n_samps)
%                     PLV(:,n_samps(nl):size(phase_data,1))=abs((sum(exp(1i*(squeeze(phase_data(n_samps(nl):size(phase_data,1),chan_contrasts(:,1),:)-phase_data(n_samps(nl):size(phase_data,1),chan_contrasts(:,2),:)))),3))')./num_trials;
%                 else
%                     PLV(:,n_samps(nl):n_samps(nl+1)-1)=abs((sum(exp(1i*(squeeze(phase_data(n_samps(nl):n_samps(nl+1)-1,chan_contrasts(:,1),:)-phase_data(n_samps(nl):n_samps(nl+1)-1,chan_contrasts(:,2),:)))),3))')./num_trials;
%                 end
%             catch ME
%                 keyboard
%             end
%         end
%     end
else
    PLV=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(:,1),:)-phase_data(:,chan_contrasts(:,2),:)))),3))')./num_trials;
end



%%% Calculating surrogate PLV data
if surg_flag ==1;
    %    keyboard;
    %PLV_surg=zeros(num_resamps,num_samp,num_chan);
    PLV_surg=single(zeros(num_resamps,size(chan_contrasts,1),num_samp));
%    fprintf(1,'Calculating %.f Surrogate PLV data. Please wait ...\n',num_resamps);
    %PLV_surg=single(zeros(num_resamps,size(chan_contrasts,1),num_samp));
    needed_mem=2*size(chan_contrasts,1)*num_samp*num_trials*8;  % 8bytse per number
    [~,sys_mem]=memory;  free_mem=sys_mem.PhysicalMemory.Available; free_mem=free_mem*.35; % available memory in bytes minus 500 MBytes

    for np=1:num_resamps
%         fprintf(1,'Surrogate # %.f\n',np);
        %tic
        if size(chan_contrasts,1)==1 || size(phase_data,1)==1;
            PLV_surg(np,:,:)=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(:,1),:)-phase_data(:,chan_contrasts(:,2),randperm(size(phase_data,3)))))),2))')./num_trials;
            
%         elseif size(chan_contrasts,1)*num_samp*num_trials>8e8     % 800 million is an arbitrary limit to reduce memory consumption
%             fprintf('"chan_contrasts" too large to calculate PLV_surg! subdividing into blocks of chan_contrasts.\n');
%             blk2_size=floor(8e8/(num_samp*num_trials));
        elseif needed_mem>free_mem
%             fprintf('"chan_contrasts" too large to calculate PLV_surg! subdividing into blocks of chan_contrasts.\n');
            %blk2_size=floor(8e8/(num_samp*num_trials));
            blk2_size=floor(free_mem/(num_samp*num_trials*8));
            n_chans=[1:blk2_size:size(chan_contrasts,1)];
            for nl=1:length(n_chans);
                try
%                     fprintf('Calculating PLV_surg in blocks of %.f samples. Block = %.f of %.f\n',num_samp,nl,length(n_chans));
                    if nl==length(n_chans)
                        PLV_surg(np,n_chans(nl):size(chan_contrasts,1),:)=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(n_chans(nl):size(chan_contrasts,1),1),:)-phase_data(:,chan_contrasts(n_chans(nl):size(chan_contrasts,1),2),randperm(size(phase_data,3)))))),3))')./num_trials;
                    else
                        PLV_surg(np,n_chans(nl):n_chans(nl+1)-1,:)=abs((sum(exp(1i*(squeeze(phase_data(:,chan_contrasts(n_chans(nl):n_chans(nl+1)-1,1),:)-phase_data(:,chan_contrasts(n_chans(nl):n_chans(nl+1)-1,2),randperm(size(phase_data,3)))))),3))')./num_trials;
                    end
                catch ME
                    keyboard
                end
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
PLV_data.PLV_surg=PLV_surg;
