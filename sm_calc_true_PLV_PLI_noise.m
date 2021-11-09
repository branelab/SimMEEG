function sm_calc_true_PLV_PLI_noise(varargin)
% This program used matlab's "cwt.m" --> future versions will incorporate "ft_freqanalysis.m" because it has a lot more functionality and evidence for use with M/EEG
%
%

global h
calc_flag=1;

tic;
if calc_flag==1
    %% initializing variables
    plv_surg_based_mean = []; plv_surg_based_std = [];
    pli_surg_based_mean = []; pli_surg_based_std = [];
    dpli_surg_based_mean = []; dpli_surg_based_std = [];
    
    
    %% setting visibiltiy of waitfor panel
    if h.monte_carlo_flag == 1
        h.waitfor_txt.String = sprintf('Time-Frequency Analyses\nProjected Noise of True Data\n\n Calculating ...'); drawnow;
    else
        h.waitfor_panel.Visible='on';
        h.waitfor_txt.String = sprintf('Time-Frequency Analyses\nProjected Noise of True Data\n\n Calculating ...'); drawnow;
    end
    
    hm1 = msgbox(sprintf('Running Time-Frequency Analyses on True Source Locations\n\nThis can take time.\nPlease be patient.\n\nCalculating ...'));
    hm1.Position(3:4)=[300 150]; htext = findobj(hm1, 'Type', 'Text'); htext.FontSize = 10; htext.HorizontalAlignment = 'left'; drawnow; % setting fontsize to being readable
    TB = str2num(h.edit_wavelet_TB.String); % wavelet parameter --> The larger the time-bandwidth parameter, the more spread out the wavelet is in time and narrower the wavelet is in frequency.
    
    
    %% checking if freqs/octave are within 4 to 48
    if str2num(h.edit_inv_wvlt_num_voices.String)<4
        h.edit_inv_wvlt_num_voices.String = 4;
        warndlg(sprintf('Number of frequencies/octave must be between 4 and 48\n Now set to 4'),'Fres/Ocatve too low');
    elseif str2num(h.edit_inv_wvlt_num_voices.String)>48
        h.edit_inv_wvlt_num_voices.String = 48;
        warndlg(sprintf('Number of frequencies/octave must be between 4 and 48\n Now set to 48'),'Fres/Ocatve too high');
    end
    
    
    
    %% %%%%%%%%%%%%%%%%%%% Time-Frequency Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TFR & PLV/PLI parameters - for standard wavelet analyses using matlab ---> alternative would be to use "ft_freqanalysis.m" but will leave this to later because this will take some time to test conversion of data structures and parameters
%     true_swf = h.sim_data.sig_final;
    %% INVERSE projecting sensor noise to the true source locations to perform TFR analyses
    sens_data = h.sim_data.sens_noise_final;
    lf = h.anatomy.leadfield;   % using leadfield that was used to forward project the source data to sensors
    vx_idx = h.cfg.source.vx_idx; source_ori = h.cfg.source.vx_ori; 
    lf_gain = str2num(h.edit_leadfield_gain.String);
    true_swf = project_inverse_SimSignals(sens_data,h.anatomy.leadfield,vx_idx,source_ori,lf_gain);
    
    
    %% get setting from "Analysis" subpanel on "Source Modeling" Tab
    flimits = str2num(h.edit_inv_plv_freq_bands.String);
    min_max_freq=flimits;
    
    %% calculate wavelets (total & induced) under signal final
    %% "true_swf" wavelets
    clear true_wt true_wt_ind true_wt_evk;
    true_ind = bsxfun(@minus,true_swf,nanmean(true_swf,3)); % induced by subtracting mean across trials (i.e., evoked response)
    fprintf('Calculating wavelets for true source waveforms...\n')
    for v=1:size(true_swf,2)    % num true sources
        %% Wavelets - Total Power
        for t=1:size(true_swf,3)
            [true_wt(:,:,v,t),F,coi_true_wt]=cwt(squeeze(true_swf(:,v,t)),'morse',h.cfg.study.srate,'TimeBandwidth',TB,'FrequencyLimits',flimits,'VoicesPerOctave',str2num(h.edit_inv_wvlt_num_voices.String)); % total power
            [true_wt_ind(:,:,v,t),F,coi_true_wt]=cwt(squeeze(true_ind(:,v,t)),'morse',h.cfg.study.srate,'TimeBandwidth',TB,'FrequencyLimits',flimits,'VoicesPerOctave',str2num(h.edit_inv_wvlt_num_voices.String));    %induced power
        end
        [true_wt_evk(:,:,v),F,coi_true_wt]=cwt(squeeze(nanmean(true_swf(:,v,:),3)),'morse',h.cfg.study.srate,'TimeBandwidth',TB,'FrequencyLimits',flimits,'VoicesPerOctave',str2num(h.edit_inv_wvlt_num_voices.String));    % evoked power
    end
    
    %% cwt has outputs freq indices in descending order
    F2=flipud(F); true_wt2=flipud(true_wt); true_wt2_ind=flipud(true_wt_ind); true_wt2_evk=flipud(true_wt_evk);
    true_wt3=abs(true_wt2); % converting to real
    true_wt3_ind=abs(true_wt2_ind); % converting to real
    true_wt3_evk=abs(true_wt2_evk); % converting to real
    
    % Decibel
    % all
    db_true_wt = 10*bsxfun(@minus,log10(true_wt3),log10(nanmean(true_wt3(:,h.cfg.study.base_samps,:,:),2)));
    true_wt_based = bsxfun(@minus,db_true_wt,nanmean(db_true_wt(:,h.cfg.study.base_samps,:,:),2));   % decibel baselined
    % induced
    db_true_wt = 10*bsxfun(@minus,log10(true_wt3_ind),log10(nanmean(true_wt3_ind(:,h.cfg.study.base_samps,:,:),2)));
    true_wt_ind_based = bsxfun(@minus,db_true_wt,nanmean(db_true_wt(:,h.cfg.study.base_samps,:,:),2));   % decibel baselined
    % evoked
    db_true_wt = 10*bsxfun(@minus,log10(true_wt3_evk),log10(nanmean(true_wt3_evk(:,h.cfg.study.base_samps,:,:),2)));
    true_wt_evk_based = bsxfun(@minus,db_true_wt,nanmean(db_true_wt(:,h.cfg.study.base_samps,:,:),2));   % decibel baselined
    
    avg_true_wt=nanmean(true_wt_based,4);
    avg_true_wt_ind=nanmean(true_wt_ind_based,4);
    avg_true_wt_evk=nanmean(true_wt_evk_based,4);
    
    %% PLV/PLI calculations based on wavelets
    fprintf('Calculating PLV ...\n')
    
    plv_phase=angle(true_wt2);
    chan_contrasts = nchoose2(1:3);
    
    F_plv=F2;
    coi_wt2=coi_true_wt; coi_wt2(coi_true_wt>max(F2))=nan; coi_wt2(coi_true_wt<min(F2))=nan;
    
    num_resamps = str2num(h.edit_inv_plv_surrogates.String);
    if num_resamps<=0; surg_flag=0; else; surg_flag=1; end
    
    clear plv_data pli_data dpli_data plv_surg pli_surg spli_surg;
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
            PLI_win=range(h.cfg.study.lat_sim)/50; PLI_win_overlap=PLI_win/2;
            [PLI]=sm_calc_PLI_ath(squeeze(plv_phase(f,:,:,:)),h.cfg.study.srate,h.cfg.study.lat_sim,PLI_win,PLI_win_overlap,chan_contrasts,surg_flag,num_resamps);
            pli_data(f,:,:)=PLI.PLI; dpli_data(f,:,:)=PLI.dPLI;
            pli_surg=PLI.PLI_surg; dpli_surg=PLI.dPLI_surg;
            pli_lat=PLI.lat;  ss=find(pli_lat<=0); ss=find(pli_lat<=h.cfg.study.base_int(1)); if isempty(ss); bs(1)=1; else; bs(1)=ss(end); end
            ss=find(pli_lat<=h.cfg.study.base_int(2)); bs(2)=ss(end);
            base_samps_pli=bs(1):bs(2);
        else
            pli_data = []; pli_surg=[];
            dpli_data = []; dpli_surg=[];
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
    h.cfg.source.TFR_results.Noise.plv_surg_based_mean = plv_surg_based_mean; clear plv_surg_based_mean;
    h.cfg.source.TFR_results.Noise.plv_surg_based_std = plv_surg_based_std; clear plv_surg_based_std;
    h.cfg.source.TFR_results.Noise.pli_surg_based_mean = pli_surg_based_mean; clear plv_surg_based_mean;
    h.cfg.source.TFR_results.Noise.pli_surg_based_std = pli_surg_based_std; clear pli_surg_based_std;
    h.cfg.source.TFR_results.Noise.dpli_surg_based_mean = dpli_surg_based_mean; clear dpli_surg_based_mean
    h.cfg.source.TFR_results.Noise.dpli_surg_based_std = dpli_surg_based_std; clear dpli_surg_based_std;
    
    % exporting to cfg.source
    h.cfg.source.TFR_results.Noise.avg_true_wt = avg_true_wt; h.cfg.source.TFR_results.Noise.avg_true_wt_evk = avg_true_wt_evk; h.cfg.source.TFR_results.Noise.avg_true_wt_ind = avg_true_wt_ind;
    if any( strcmp(h.listbox_inv_FC_type.String(h.listbox_inv_FC_type.Value),'PLV') )
        plv_based=bsxfun(@minus,plv_data,nanmean(plv_data(:,:,h.cfg.study.base_samps),3));
        h.cfg.source.TFR_results.Noise.plv_based = plv_based;
        h.cfg.source.TFR_results.Noise.plv_data = plv_data;
    end
    if any( strcmp(h.listbox_inv_FC_type.String(h.listbox_inv_FC_type.Value),'PLI') )
        pli_based=bsxfun(@minus,pli_data,nanmean(pli_data(:,:,base_samps_pli),3));
        dpli_based=bsxfun(@minus,dpli_data,nanmean(dpli_data(:,:,base_samps_pli),3));
        h.cfg.source.TFR_results.Noise.pli_data = pli_data;
        h.cfg.source.TFR_results.Noise.pli_based = pli_based; h.cfg.source.TFR_results.Noise.dpli_based = dpli_based;
        h.cfg.source.TFR_results.Noise.pli_lat = pli_lat;
    end
    h.cfg.source.TFR_results.Noise.TFR_freqs = F2;
    h.cfg.source.TFR_results.Noise.coi_wt2 = coi_wt2;
    h.cfg.source.TFR_results.Noise.PLV_freqs = F2;
    
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
