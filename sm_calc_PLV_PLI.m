function sm_calc_PLV_PLI(varargin)
% This program used matlab's "cwt.m" --> future versions will incorporate "ft_freqanalysis.m" because it has a lot more functionality and evidence for use with M/EEG
%
%

global h

%% no peak sources found
if all(isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.Hits))
    h.inv_soln(h.current_inv_soln).plv_seed_idx = h.inv_soln(h.current_inv_soln).classifier_metrics.Hits;
    hm = warndlg(sprintf('No Peak Sources Found\nChange Image Scale --> Select Seed and Comparions Locations'),'No Sources Found '); %WinOnTop(hm);
    hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
    return
end

%% Initialize varibles
plv_surg_based_mean = []; 
plv_surg_based_std = []; 
pli_surg_based_mean = []; 
pli_surg_based_std = []; 
dpli_surg_based_mean = []; 
dpli_surg_based_std = []; 


%%
calc_flag=0;
if isfield(h.inv_soln(h.current_inv_soln),'TFR_results') 
    if isempty(h.inv_soln(h.current_inv_soln).TFR_results)
        calc_flag=1;
    else
        answ = questdlg(sprintf('Peak Source Connectivity Exists\n\nWould you like to recalculate?'),'Peak Connectivity Exists','Yes','No','No');
%      answ = 'Yes';
     switch answ
            case 'Yes'
                calc_flag=1;
            case 'No'
                calc_flag=0;
        end
    end
else
    calc_flag=1;
end


if isempty(h.inv_soln(h.current_inv_soln).plv_seed_idx)   % no seed peaks found to conduct PLV analyses
    calc_flag=0;
    h.inv_soln(h.current_inv_soln).TFR_results=[];
    fprintf('No seed voxels selected or found. To conduct PLV/PLI analysis\n'); 
end

if calc_flag==1
    tic;
    
    
    
    %% %%%%% gathering parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% checking memory is suficient for # FC contrasts
    % num_contrasts = size(h.inv_soln(h.current_inv_soln).plv_contrasts,1);
    % num_freqs=60; num_samps = length(h.sim_data.cfg.study.lat_sim);
    % num_numbers = num_contrasts*num_freqs*num_samps*8;
    %
    % perc_remain = 10;
    % [mem_flag]=memory_check(num_numbers,perc_remain);
    %
    % if mem_flag == 0   % PLV calculations will exceed memory but can still proceed because "calc_PLV" will block the plv_contrasts into separate blocks per calculation depending on available memory -- it will just take time
    %     answ = questdlg(sprintf('Number of PLV contrasts is %.f\n\nThis will take a long time for PLV to compute\nand may exceed computer''s memory capacity\n\nWould you like to continue?\n',num_contrasts),'Save SimMEEG Dataset?','Yes','No','No');
    % else
    %     answ = 'Yes';
    % end
    
    %% setting visibiltiy of waitfor panel
    if h.monte_carlo_flag == 1
        h.waitfor_txt.String = sprintf('Time-Frequency Analyses of Peak Data\n\n Calculating ...'); drawnow;
    else
        h.waitfor_panel.Visible='on';
        h.waitfor_txt.String = sprintf('Time-Frequency Analyses of Peak Data\n\n Calculating ...'); drawnow;
    end
    
    hm1 = msgbox(sprintf('Running Time-Frequency Analyses on selected\nSeed and Comparison Locations\n\nThis can take time.\nPlease be patient.\n\nCalculating ...'));
    hm1.Position(3:4)=[300 130]; htext = findobj(hm1, 'Type', 'Text'); htext.FontSize = 11; htext.HorizontalAlignment = 'left'; drawnow; % setting fontsize to being readable
    TB = str2num(h.edit_wavelet_TB.String); % wavelet parameter --> The larger the time-bandwidth parameter, the more spread out the wavelet is in time and narrower the wavelet is in frequency.
    
    % getting config params for simulated data
    cfg=h.sim_data.cfg;
    
    %% checking if freqs/octave are within 4 to 48
    if str2num(h.edit_inv_wvlt_num_voices.String)<4
        h.edit_inv_wvlt_num_voices.String = 4;
        warndlg(sprintf('Number of frequencies/octave must be between 4 and 48\n Now set to 4'),'Fres/Ocatve too low');
    elseif str2num(h.edit_inv_wvlt_num_voices.String)>48
        h.edit_inv_wvlt_num_voices.String = 48;
        warndlg(sprintf('Number of frequencies/octave must be between 4 and 48\n Now set to 48'),'Fres/Ocatve too high');
    end
    
    
    %% %%%% create source waveforms (swf) depending on invSoln type
    % Issue: How to deal with vector-beamformers? --> currently set to using "max" power. This is because it is computationally explosive to include all 3 dipole orientations for
    % connectivity analyses using vector beamformers.
    %
    
    if size(h.inv_soln(h.current_inv_soln).plv_contrasts,1)>1
        
        seed_idx = h.inv_soln(h.current_inv_soln).plv_seed_idx;     % seed dipole indices of current lead fields
        comp_idx  = h.inv_soln(h.current_inv_soln).plv_comp_idx;    % comparison dipole indices of current lead fields
        plv_idx = h.inv_soln(h.current_inv_soln).plv_contrast_idx;  % comparison indices for contrasts of "seed_swf"  and "comp_swf" below; clmn 1 = seed_idx;  clmn 2 = comp_idx
        seed_val_idx = find(~isnan(seed_idx));
        comp_val_idx = find(~isnan(comp_idx));
        switch h.inv_soln(h.current_inv_soln).Type
            case {'Dipole', 'SPA','LCMV (FT)', 'SAM (FT)','sLORETA (FT)','sMCMV','bRAPBeam','TrapMUSIC'}    % BRANE Lab beamformers
                seed_swf = nan(size(h.sim_data.sens_final,1), length(seed_idx),size(h.sim_data.sens_final,3)); 
                comp_swf = nan(size(h.sim_data.sens_final,1), length(comp_idx),size(h.sim_data.sens_final,3)); 
                for t=1:size(h.sim_data.sens_final,3)
                    seed_swf(:,seed_val_idx,t) = [h.inv_soln(h.current_inv_soln).soln.wts(:,seed_idx(seed_val_idx))'*squeeze(h.sim_data.sens_final(:,:,t))']';
                    seed_swf(:,:,t) = bsxfun(@minus, seed_swf(:,:,t), nanmean(seed_swf(h.sim_data.cfg.study.base_samps,:,t)));
                    comp_swf(:,comp_val_idx,t) = [h.inv_soln(h.current_inv_soln).soln.wts(:,comp_idx(comp_val_idx))'*squeeze(h.sim_data.sens_final(:,:,t))']';
                    comp_swf(:,:,t) = bsxfun(@minus, comp_swf(:,:,t), nanmean(comp_swf(h.sim_data.cfg.study.base_samps,:,t)));
                end
            case {'SIA','MIA'}    % BRANE Lab beamformers - have to use nulled_wts or plv_locs that were not in wts matrix will have nans
                seed_swf = nan(size(h.sim_data.sens_final,1), length(seed_idx),size(h.sim_data.sens_final,3));
                comp_swf = nan(size(h.sim_data.sens_final,1), length(comp_idx),size(h.sim_data.sens_final,3));
                for t=1:size(h.sim_data.sens_final,3)
                    %                 % % using nulled wts calculated from Ninv Noise Covarianace may get better supression of noise because it doesn't include possible active source as in Rinv
                    seed_swf(:,seed_val_idx,t) = [h.inv_soln(h.current_inv_soln).soln.nulled_wts(:,seed_idx(seed_val_idx))'*squeeze(h.sim_data.sens_final(:,:,t))']';
                    seed_swf(:,seed_val_idx,t) = bsxfun(@minus, seed_swf(:,seed_val_idx,t), nanmean(seed_swf(h.sim_data.cfg.study.base_samps,seed_val_idx,t)));
                    comp_swf(:,comp_val_idx,t) = [h.inv_soln(h.current_inv_soln).soln.nulled_wts(:,comp_idx(comp_val_idx))'*squeeze(h.sim_data.sens_final(:,:,t))']';
                    comp_swf(:,comp_val_idx,t) = bsxfun(@minus, comp_swf(:,comp_val_idx,t), nanmean(comp_swf(h.sim_data.cfg.study.base_samps,comp_val_idx,t)));
                end
            case {'eLORETA (FT)' 'MNE (FT)' 'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}    % vector inverse solutions
                % picking orientation with maximal response in active interval to generate a source waveform
                clear seed_swf comp_swf;
                % re-order wts in case different dimensions --> wts dims should be [chans x voxels x ori]
                dims1 = size(h.inv_soln(h.current_inv_soln).soln.wts);
                dims2 = [length(h.inv_soln(h.current_inv_soln).leadfield.label), size(h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,1) 3];
                if ~isequal(dims1,dims2) % reordering wts
                    dims=[]; for a=1:3; dims(a) = find(dims1(a)==dims2); end
                    wts = permute( h.inv_soln(h.current_inv_soln).soln.wts,dims);
                else
                    wts = h.inv_soln(h.current_inv_soln).soln.wts;
                end
                
                % source waveforms for 3-vector dipoles per voxel
                swf=[];
                for ox=1:size(wts,3)
                    seed_swf(:,ox,seed_val_idx) = (squeeze(wts(:,seed_idx(seed_val_idx),ox))'*squeeze(nanmean(h.sim_data.sens_final,3))')';
                    seed_swf(:,ox,seed_val_idx) = bsxfun(@minus, seed_swf(:,ox,seed_val_idx), nanmean(seed_swf(h.sim_data.cfg.study.base_samps,ox,seed_val_idx)));
                    comp_swf(:,ox,comp_val_idx) = (squeeze(wts(:,comp_idx(comp_val_idx),ox))'*squeeze(nanmean(h.sim_data.sens_final,3))')';
                    comp_swf(:,ox,comp_val_idx) = bsxfun(@minus, comp_swf(:,ox,comp_val_idx), nanmean(comp_swf(h.sim_data.cfg.study.base_samps,ox,comp_val_idx)));
                end
                % For vector inv_soln --> selecting dipole orientation (X,Y,or Z) that has max RMS across active interval time samples 'act_samp'
                try
                    [~,seed_max_ori] = nanmax(squeeze(rms(seed_swf(h.sim_data.cfg.study.bl_bmf.act_samps,:,:),1)));
                    [~,comp_max_ori] = nanmax(squeeze(rms(comp_swf(h.sim_data.cfg.study.bl_bmf.act_samps,:,:),1)));
                catch
                    h.sim_data.cfg.study.bl_bmf = h.cfg.study.bl_bmf;
                    [~,seed_max_ori] = nanmax(squeeze(rms(seed_swf(h.sim_data.cfg.study.bl_bmf.act_samps,:,:),1)));
                    [~,comp_max_ori] = nanmax(squeeze(rms(comp_swf(h.sim_data.cfg.study.bl_bmf.act_samps,:,:),1)));
                end
                
                seed_swf = nan(size(h.sim_data.sens_final,1), length(seed_idx),size(h.sim_data.sens_final,3));
                comp_swf = nan(size(h.sim_data.sens_final,1), length(comp_idx),size(h.sim_data.sens_final,3));
                for v=1:length(seed_idx)
                    if ~isnan(seed_idx(v))
                        for t=1:size(h.sim_data.sens_final,3)
                            seed_swf(:,v,t) = [squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,seed_idx(v),seed_max_ori(v)))'*squeeze(h.sim_data.sens_final(:,:,t))']';
                        end
                    end
                end
                
                for v=1:length(comp_idx)
                    if ~isnan(comp_idx(v))
                        for t=1:size(h.sim_data.sens_final,3)
                            comp_swf(:,v,t) = [squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,comp_idx(v),comp_max_ori(v)))'*squeeze(h.sim_data.sens_final(:,:,t))']';
                        end
                    end
                end
                
                seed_swf = bsxfun(@minus, seed_swf, nanmean(seed_swf(h.sim_data.cfg.study.base_samps,:,:)));
                comp_swf = bsxfun(@minus, comp_swf, nanmean(comp_swf(h.sim_data.cfg.study.base_samps,:,:)));
       end
        
        %% %%%%%%%%%%%%%%%%%%% Time-Frequency Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TFR & PLV/PLI parameters - for standard wavelet analyses using matlab ---> alternative would be to use "ft_freqanalysis.m" but will leave this to later because this will take some time to test conversion of data structures and parameters
        
        % get setting from "Analysis" subpanel on "Source Modeling" Tab
        flimits = str2num(h.edit_inv_plv_freq_bands.String);
        min_max_freq=flimits;
        
        %% calculate wavelets (total & induced) under signal final
        %% "seed_swf" wavelets
        mat_size = cwt(squeeze(seed_swf(:,seed_val_idx(1),t)),'morse',cfg.study.srate,'TimeBandwidth',TB,'FrequencyLimits',flimits,'VoicesPerOctave',str2num(h.edit_inv_wvlt_num_voices.String)); % total power
        dims =size(mat_size);
        seed_wt = nan(dims(1),dims(2),length(seed_idx));
        seed_wt_ind = seed_wt; seed_wt_evk = seed_wt; 
        seed_ind = bsxfun(@minus,seed_swf,nanmean(seed_swf,3)); % induced by subtracting mean across trials (i.e., evoked response)
        fprintf('Calculating wavelets for seed source waveforms...\n')
        hw = waitbar(0,'Calculating Wavelets for each Seed Location');
        for v=1:size(seed_swf,2)    % num seed sources
            waitbar(v/size(seed_swf,2),hw);
            if ~isnan(seed_idx(v))
                %% Wavelets - Total Power
                for t=1:size(seed_swf,3)
                    [seed_wt(:,:,v,t),F,coi_seed_wt]=cwt(squeeze(seed_swf(:,v,t)),'morse',cfg.study.srate,'TimeBandwidth',TB,'FrequencyLimits',flimits,'VoicesPerOctave',str2num(h.edit_inv_wvlt_num_voices.String)); % total power
                    [seed_wt_ind(:,:,v,t),F,coi_seed_wt]=cwt(squeeze(seed_ind(:,v,t)),'morse',cfg.study.srate,'TimeBandwidth',TB,'FrequencyLimits',flimits,'VoicesPerOctave',str2num(h.edit_inv_wvlt_num_voices.String));    %induced power
                end
                [seed_wt_evk(:,:,v),F,coi_seed_wt]=cwt(squeeze(nanmean(seed_swf(:,v,:),3)),'morse',cfg.study.srate,'TimeBandwidth',TB,'FrequencyLimits',flimits,'VoicesPerOctave',str2num(h.edit_inv_wvlt_num_voices.String));    % evoked power
            end
        end
        delete(hw);
        %% cwt has outputs freq indices in descending order
        F2=flipud(F); seed_wt2=single(flipud(seed_wt)); seed_wt2_ind=single(flipud(seed_wt_ind)); seed_wt2_evk=single(flipud(seed_wt_evk));
         clear seed_wt seed_wt_ind seed_wt_evk;
         % Decibel
        % all
        db_seed_wt = 10*bsxfun(@minus,log10(abs(seed_wt2)),log10(nanmean(abs(seed_wt2(:,cfg.study.base_samps,:,:)),2)));
        seed_wt_based = bsxfun(@minus,db_seed_wt,nanmean(db_seed_wt(:,cfg.study.base_samps,:,:),2));   % decibel baselined
        % induced
        db_seed_wt = 10*bsxfun(@minus,log10(abs(seed_wt2_ind)),log10(nanmean(abs(seed_wt2_ind(:,cfg.study.base_samps,:,:)),2)));
        seed_wt_ind_based = bsxfun(@minus,db_seed_wt,nanmean(db_seed_wt(:,cfg.study.base_samps,:,:),2));   % decibel baselined
        % evoked
        db_seed_wt = 10*bsxfun(@minus,log10(abs(seed_wt2_evk)),log10(nanmean(abs(seed_wt2_evk(:,cfg.study.base_samps,:,:)),2)));
        seed_wt_evk_based = bsxfun(@minus,db_seed_wt,nanmean(db_seed_wt(:,cfg.study.base_samps,:,:),2));   % decibel baselined
              
        avg_seed_wt=nanmean(seed_wt_based,4);
        avg_seed_wt_ind=nanmean(seed_wt_ind_based,4);
        avg_seed_wt_evk=nanmean(seed_wt_evk_based,4);
        
        
        %% "comp_swf" wavelets
        clear comp_wt comp_wt_ind comp_wt_evk;
        comp_ind = bsxfun(@minus,comp_swf,nanmean(comp_swf,3)); % induced by subtracting mean across trials (i.e., evoked response)
        fprintf('Calculating wavelets for comparison source waveforms...\n')
        hw = waitbar(0,'Calculating Wavelets for each Comparison Location');
        for v=1:size(comp_swf,2)    % num comp sources
            waitbar(v/size(comp_swf,2),hw);
            if ~isnan(comp_idx(v))
                %% Wavelets - Total Power
            for t=1:size(comp_swf,3)
                [comp_wt(:,:,v,t),F,coi_comp_wt]=cwt(squeeze(comp_swf(:,v,t)),'morse',cfg.study.srate,'TimeBandwidth',TB,'FrequencyLimits',flimits,'VoicesPerOctave',str2num(h.edit_inv_wvlt_num_voices.String)); % total power
                [comp_wt_ind(:,:,v,t),F,coi_comp_wt]=cwt(squeeze(comp_ind(:,v,t)),'morse',cfg.study.srate,'TimeBandwidth',TB,'FrequencyLimits',flimits,'VoicesPerOctave',str2num(h.edit_inv_wvlt_num_voices.String));    %induced power
            end
            [comp_wt_evk(:,:,v),F,coi_comp_wt]=cwt(squeeze(nanmean(comp_swf(:,v,:),3)),'morse',cfg.study.srate,'TimeBandwidth',TB,'FrequencyLimits',flimits,'VoicesPerOctave',str2num(h.edit_inv_wvlt_num_voices.String));    % evoked power
            end
        end
        delete(hw);
        
        %% cwt has outputs freq indices in descending order
        F2=flipud(F); 
        comp_wt=single(comp_wt(fliplr(1:size(comp_wt,1)),:,:,:)); 
        comp_wt_ind=single(comp_wt_ind(fliplr(1:size(comp_wt_ind,1)),:,:,:)); 
        comp_wt_evk=single(comp_wt_evk(fliplr(1:size(comp_wt_evk,1)),:,:)); 
        
%         clear comp_wt comp_wt_ind comp_wt_evk;
        
        % Decibel
        % all
        db_comp_wt = 10*bsxfun(@minus,log10(abs(comp_wt)),log10(nanmean(abs(comp_wt(:,cfg.study.base_samps,:,:)),2)));
        comp_wt_based = bsxfun(@minus,db_comp_wt,nanmean(db_comp_wt(:,cfg.study.base_samps,:,:),2));   % decibel baselined
        % induced
        db_comp_wt = 10*bsxfun(@minus,log10(abs(comp_wt_ind)),log10(nanmean(abs(comp_wt_ind(:,cfg.study.base_samps,:,:)),2)));
        comp_wt_ind_based = bsxfun(@minus,db_comp_wt,nanmean(db_comp_wt(:,cfg.study.base_samps,:,:),2));   % decibel baselined
        % evoked
        db_comp_wt = 10*bsxfun(@minus,log10(abs(comp_wt_evk)),log10(nanmean(abs(comp_wt_evk(:,cfg.study.base_samps,:,:)),2)));
        comp_wt_evk_based = bsxfun(@minus,db_comp_wt,nanmean(db_comp_wt(:,cfg.study.base_samps,:,:),2));   % decibel baselined
        
        avg_comp_wt=nanmean(comp_wt_based,4);
        avg_comp_wt_ind=nanmean(comp_wt_ind_based,4);
        avg_comp_wt_evk=nanmean(comp_wt_evk_based,4);
    end
    
    %% PLV/PLI calculations based on wavelets
    fprintf('Calculating PLV ...\n')
    
    seed_phase_data=angle(seed_wt2);
    comp_phase_data=angle(comp_wt);
%     plv_phase = cat(3,seed_phase_data,comp_phase_data); % concatenating phases to work with calc_PLV_ath.m
    plv_phase = comp_phase_data; 
    plv_phase(:,:,1:size(seed_phase_data,3),:) = seed_phase_data; % concatenating phases to work with calc_PLV_ath.m
    
    chan_contrasts = plv_idx;
%     chan_contrasts(:,2) = chan_contrasts(:,2)+length(seed_idx);   % making matrix compatible with plv_phase for calc_PLV_ath.m
    
    F_plv=F2;
    coi_wt2=coi_seed_wt; coi_wt2(coi_seed_wt>max(F2))=nan; coi_wt2(coi_seed_wt<min(F2))=nan;
    
    
    num_resamps = str2num(h.edit_inv_plv_surrogates.String);
    if num_resamps<=0; surg_flag=0; else; surg_flag=1; end
    
    clear plv_data pli_data dpli_data;
    fprintf('Calculating Connectivity Analysis for:\n'); fprintf('%s\n',h.listbox_inv_FC_type.String{h.listbox_inv_FC_type.Value}); fprintf('Please wait ...\n');
    hw = waitbar(0,'Running connectivitiy analysis across frequencies');
    for f=1:size(plv_phase,1)
        waitbar(f/size(plv_phase,1),hw)
        
        if any( strcmp(h.listbox_inv_FC_type.String(h.listbox_inv_FC_type.Value),'PLV') )
            [PLV]=sm_calc_PLV_ath(squeeze(plv_phase(f,:,:,:)),chan_contrasts,surg_flag,num_resamps);
            plv_data(f,:,:)=PLV.PLV;
            plv_surg=PLV.PLV_surg;
        else
            plv_data = []; plv_surg = [];
        end
        if any( strcmp(h.listbox_inv_FC_type.String(h.listbox_inv_FC_type.Value),'PLI') )
            PLI_win=range(cfg.study.lat_sim)/50; PLI_win_overlap=PLI_win/2;
            [PLI]=sm_calc_PLI_ath(squeeze(plv_phase(f,:,:,:)),cfg.study.srate,cfg.study.lat_sim,PLI_win,PLI_win_overlap,chan_contrasts,surg_flag,num_resamps);
            pli_data(f,:,:)=PLI.PLI; dpli_data(f,:,:)=PLI.dPLI;
            pli_surg=PLI.PLI_surg; dpli_surg=PLI.dPLI_surg;
             pli_lat=PLI.lat; ss=find(pli_lat<=0); ss=find(pli_lat<=h.cfg.study.base_int(1)); if isempty(ss); bs(1)=1; else; bs(1)=ss(end); end
            ss=find(pli_lat<=h.cfg.study.base_int(2)); bs(2)=ss(end);
            base_samps_pli=bs(1):bs(2);
        else
            pli_data = []; pli_surg = [];
            dpli_data = []; dpli_surg = [];
        end
        
           % calculating PLV/PLI surrogate stats
        if ~isempty(plv_surg)
            plv_surg_based = bsxfun(@minus,plv_surg,nanmean(plv_surg(:,:,h.cfg.study.base_samps),3)); % baselining
           plv_surg_based_mean(f,:,:) = squeeze(nanmean(plv_surg_based,1));
           plv_surg_based_std(f,:,:) = squeeze(nanstd(plv_surg_based,[],1));
        end
        if ~isempty(pli_surg) 
            pli_surg_based = bsxfun(@minus,pli_surg,nanmean(pli_surg(:,:,base_samps_pli),3)); % baselining
           pli_surg_based_mean(f,:,:) = squeeze(nanmean(pli_surg_based,1));
           pli_surg_based_std(f,:,:) = squeeze(nanstd(pli_surg_based,[],1));
            dpli_surg_based = bsxfun(@minus,dpli_surg,nanmean(pli_surg(:,:,base_samps_pli),3)); % baselining
           dpli_surg_based_mean(f,:,:) = squeeze(nanmean(dpli_surg_based,1));
           dpli_surg_based_std(f,:,:) = squeeze(nanstd(dpli_surg_based,[],1));
        end
        
    end
    delete(hw);
    h.inv_soln(h.current_inv_soln).TFR_results.plv_surg_based_mean = plv_surg_based_mean; clear plv_surg_based_mean;
    h.inv_soln(h.current_inv_soln).TFR_results.plv_surg_based_std = plv_surg_based_std; clear plv_surg_based_std;
    h.inv_soln(h.current_inv_soln).TFR_results.pli_surg_based_mean = pli_surg_based_mean; clear pli_surg_based_mean
    h.inv_soln(h.current_inv_soln).TFR_results.pli_surg_based_std = pli_surg_based_std; clear pli_surg_based_std;
    h.inv_soln(h.current_inv_soln).TFR_results.dpli_surg_based_mean = dpli_surg_based_mean; clear dpli_surg_based_mean
    h.inv_soln(h.current_inv_soln).TFR_results.dpli_surg_based_std = dpli_surg_based_std; clear dpli_surg_based_std;

    
    %% exporting to inv_soln
    % Power
    h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt = avg_seed_wt; h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_evk = avg_seed_wt_evk; h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_ind = avg_seed_wt_ind;
    h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt = avg_comp_wt; h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_evk = avg_comp_wt_evk; h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_ind = avg_comp_wt_ind;
    
    
    % calculate PLV for porjected sensor noise to be used in statistics
    if h.radio_3D_include_noise_PLV.Value==1 && calc_flag==1; sm_calc_PLV_PLI_noise;  end 

    
    % PLV
    if any( strcmp(h.listbox_inv_FC_type.String(h.listbox_inv_FC_type.Value),'PLV') )
        plv_based=bsxfun(@minus,plv_data,nanmean(plv_data(:,:,cfg.study.base_samps),3));
        h.inv_soln(h.current_inv_soln).TFR_results.plv_data = plv_data; h.inv_soln(h.current_inv_soln).TFR_results.pli_data = pli_data;
        h.inv_soln(h.current_inv_soln).TFR_results.plv_based = plv_based;
    end
    
    % PLI
    if any( strcmp(h.listbox_inv_FC_type.String(h.listbox_inv_FC_type.Value),'PLI') )
        pli_based=bsxfun(@minus,pli_data,nanmean(pli_data(:,:,base_samps_pli),3));
        dpli_based=bsxfun(@minus,dpli_data,nanmean(dpli_data(:,:,base_samps_pli),3));
        h.inv_soln(h.current_inv_soln).TFR_results.pli_based = pli_based; h.inv_soln(h.current_inv_soln).TFR_results.dpli_based = dpli_based;
        h.inv_soln(h.current_inv_soln).TFR_results.pli_lat = pli_lat;
    end
    h.inv_soln(h.current_inv_soln).TFR_results.TFR_freqs = F2;
    h.inv_soln(h.current_inv_soln).TFR_results.coi_wt2 = coi_wt2;
    h.inv_soln(h.current_inv_soln).TFR_results.PLV_freqs = F2;
    
    
   %% setting current contrasts for plotting and selecting TFR/PLVs
   sm_update_plv_contrast_order();
    %%
    fprintf('Analysis Completed! Time to compute = %.f sec\n',toc);
end
%%

if exist('hm','var'); delete(hm); end
if exist('hm1','var'); delete(hm1); end
if exist('hm2','var'); delete(hm2); end



if h.monte_carlo_flag ~= 1
    h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
end


return
