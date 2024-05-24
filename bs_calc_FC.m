function bs_calc_FC(varargin)
global h

if h.monte_carlo_flag == 1
    h.waitfor_txt.String = sprintf('Time-Frequency Analyses of Peak Data\n\n Calculating ...'); drawnow;
else
    h.waitfor_panel.Visible='on';
    h.waitfor_txt.String = sprintf('Time-Frequency Analyses of Peak Data\n\n Calculating ...'); drawnow;
end

hm1 = msgbox(sprintf('Running Time-Frequency Analyses of Peak Data\n\n Calculating ...')); 

if length(unique(h.current_3D_peak_idx(1:3)))<3    % means that there is replication of the found peak in the peaks found
    vx_pos = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx(1:3),:);
    hm = msgbox(sprintf('Peaks nearest to the True sources are the same\n   Nearest to Source 1 = %.f %.f %.f mm (index=%.f)\n   Nearest to Source 2 = %.f %.f %.f mm (index=%.f)\n   Nearest to Source 3 = %.f %.f %.f mm (index=%.f)\n\nTime-Frequency results will be incorrect for the repeated peaks\n',vx_pos(1,:), h.current_3D_peak_idx(1), vx_pos(2,:), h.current_3D_peak_idx(2),vx_pos(3,:), h.current_3D_peak_idx(3)));
end

TB = str2num(h.edit_wavelet_TB.String); % wavelet parameter --> The larger the time-bandwidth parameter, the more spread out the wavelet is in time and narrower the wavelet is in frequency.

% setting peaks to true source locations
% h.current_3D_peak_idx = h.sim_data.cfg.source.vx_idx;

cfg=h.sim_data.cfg;
%% Plotting --> Confirming
if length(h.current_3D_peak_idx)>3
       hm2=msgbox(sprintf('\nSelected closest peak locations to True Sources.\n\nConducting PLV/PLI analyses.\n')); 
    h.current_3D_peak_idx = h.current_3D_peak_idx(1:length(h.sim_data.cfg.source.vx_idx)); 
end

if ~isempty(h.current_3D_peak_idx) && length(h.current_3D_peak_idx)<=3
   
    switch h.inv_soln(h.current_inv_soln).Type
        case  {'SPA' 'LCMV (FT)' 'SAM (FT)' 'sLORETA (FT)' 'dics (FT)' 'pcc (FT)' 'SAM (FT)'}     % scalar inverse solutions
            h.current_peak_swf_trials = nan(size(h.sim_data.sens_final,1), length(h.current_3D_peak_idx),size(h.sim_data.sens_final,3));
            for t=1:size(h.sim_data.sens_final,3)
                h.current_peak_swf_trials(:,:,t) = [h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx)'*squeeze(h.sim_data.sens_final(:,:,t))']';
                h.current_peak_swf_trials(:,:,t) = bsxfun(@minus, h.current_peak_swf_trials(:,:,t), nanmean(h.current_peak_swf_trials(h.sim_data.cfg.study.base_samps,:,t)));
            end
        case {'SIA' 'MIA''sMCMV' 'bRAPBeam' 'TrapMUSIC'}    % multi-source beamformers
            h.current_peak_swf_trials = nan(size(h.sim_data.sens_final,1), length(h.current_3D_peak_idx),size(h.sim_data.sens_final,3));
            for t=1:size(h.sim_data.sens_final,3)
                h.current_peak_swf_trials(:,:,t) = [h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx)'*squeeze(h.sim_data.sens_final(:,:,t))']';
%                 % % using nulled wts calculated from Ninv Noise Covarianace to get better supression of noise because it doesn't include possible active source as in Rinv
%                 h.current_peak_swf_trials(:,:,t) = [h.inv_soln(h.current_inv_soln).soln.nulled_wts(:,h.current_3D_peak_idx)'*squeeze(h.sim_data.sens_final(:,:,t))']';
                h.current_peak_swf_trials(:,:,t) = bsxfun(@minus, h.current_peak_swf_trials(:,:,t), nanmean(h.current_peak_swf_trials(h.sim_data.cfg.study.base_samps,:,t)));
            end
        case {'MNE (FT)' 'eLORETA (FT)' 'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}     % vector inverse solutions
            % picking orientation with maximal response in active interval to generate a source waveform
             h.current_peak_swf = h.inv_soln(h.current_inv_soln).soln.avg.pow(h.inv_soln(h.current_inv_soln).leadfield.inside==1,:)';
            
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
%                   try
                      swf(:,ox,:) = (squeeze(wts(:,h.current_3D_peak_idx,ox))'*squeeze(nanmean(h.sim_data.sens_final,3))')';
%                   catch
%                       wts = reshape(wts,[size(wts,1)/size(h.sim_data.sens_final,2) size(h.sim_data.sens_final,2) size(wts,2)]);
%                       h.inv_soln(h.current_inv_soln).soln.wts = permute(wts,[1 3 2]);
%                       swf(:,ox,:) = [squeeze(h.inv_soln(h.current_inv_soln).soln.wts(h.current_3D_peak_idx,ox,:))*squeeze(nanmean(h.sim_data.sens_final,3))']';
%                   end
                
               swf(:,ox,:) = bsxfun(@minus, swf(:,ox,:), nanmean(swf(h.sim_data.cfg.study.base_samps,ox,:)));
              end
              
              [~,max_ori] = max(squeeze(rms(swf(h.sim_data.cfg.study.bl_bmf.act_samps,:,:),1)));
              
              h.current_peak_swf_trials = nan(size(h.sim_data.sens_final,1), length(h.current_3D_peak_idx),size(h.sim_data.sens_final,3));
              for v=1:length(h.current_3D_peak_idx)
                  for t=1:size(h.sim_data.sens_final,3)
                      h.current_peak_swf_trials(:,v,t) = [squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx(v),max_ori(v)))'*squeeze(h.sim_data.sens_final(:,:,t))']';
                  end
              end
                h.current_peak_swf_trials = bsxfun(@minus, h.current_peak_swf_trials, nanmean(h.current_peak_swf_trials(h.sim_data.cfg.study.base_samps,:,:)));
           
    end
    
   

    %% %%%%%%%%%%%%%%%%%%% Time-Frequency Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TFR & PLV/PLI parameters
    [num_chans,num_freqs,num_minmaxfr]=size(cfg.source.sig_freqs);
    lat=cfg.study.lat_sim;
    min_max_freq=cfg.study.plot_freq_int;
    %% calculate wavelets (total & induced) under signal final
    clear wt wt_ind wt_evk;
    sig_ind = bsxfun(@minus,h.current_peak_swf_trials,nanmean(h.current_peak_swf_trials,3)); % induced by subtracting mean across trials (i.e., evoked response)
    fprintf('Calculating wavelets ...\n')
    for v=1:num_chans
        %% Wavelets - Total Power
%         wt_param=[3 30]; %[3 60];
        for t=1:size(h.current_peak_swf_trials,3)
%             [wt(:,:,v,t),F,coi_wt]=cwt(squeeze(h.current_peak_swf_trials(:,v,t)),'morse',cfg.study.srate,'WaveletParameters',wt_param); % total power
%             [wt_ind(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_ind(:,v,t)),'morse',cfg.study.srate,'WaveletParameters',wt_param);    %induced power
            [wt(:,:,v,t),F,coi_wt]=cwt(squeeze(h.current_peak_swf_trials(:,v,t)),'morse',cfg.study.srate,'TimeBandwidth',TB); % total power
            [wt_ind(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_ind(:,v,t)),'morse',cfg.study.srate,'TimeBandwidth',TB);    %induced power
        end
%         [wt_evk(:,:,v),F,coi_wt]=cwt(squeeze(nanmean(h.current_peak_swf_trials(:,v,:),3)),'morse',cfg.study.srate,'WaveletParameters',wt_param);    % evoked power
        [wt_evk(:,:,v),F,coi_wt]=cwt(squeeze(nanmean(h.current_peak_swf_trials(:,v,:),3)),'morse',cfg.study.srate,'TimeBandwidth',TB);    % evoked power
    end
    F2=flipud(F); wt2=flipud(wt); wt2_ind=flipud(wt_ind); wt2_evk=flipud(wt_evk);
    ss=find(cfg.study.lat_sim<=cfg.study.base_int(1)); bs(1)=ss(end);
    ss=find(cfg.study.lat_sim<=cfg.study.base_int(2)); bs(2)=ss(end);
    base_samps=bs(1):bs(2);
    wt3=abs(wt2); % converting to real 
    wt3_ind=abs(wt2_ind); % converting to real
    wt3_evk=abs(wt2_evk); % converting to real
    %     wt_based=20*log10(bsxfun(@rdivide,wt3,nanmean(wt3(:,base_samps,:),2))); % dB
%     % dividing by baseline then multiply 100 to get percent then baseline
%     wt_based=bsxfun(@rdivide,wt3,nanmean(wt3(:,base_samps,:,:),2))*100;   % percentage
%     wt_ind_based=bsxfun(@rdivide,wt3_ind,nanmean(wt3_ind(:,base_samps,:,:),2))*100;   % percentage
%     wt_evk_based=bsxfun(@rdivide,wt3_evk,nanmean(nanmean(wt3(:,base_samps,:),2),3))*100;   % percentage
%     % baselining
%     wt_based=bsxfun(@minus,wt_based,nanmean(wt_based(:,base_samps,:,:),2));   % percentage
%     wt_ind_based=bsxfun(@minus,wt_ind_based,nanmean(wt_ind_based(:,base_samps,:,:),2));   % percentage
%     wt_evk_based=bsxfun(@minus,wt_evk_based,nanmean(wt_evk_based(:,base_samps,:),2));   % percentage
%    
%     
    
    % Decibel
    % all
    db_wt = 10*bsxfun(@minus,log10(wt3),log10(nanmean(wt3(:,base_samps,:,:),2)));
    wt_based = bsxfun(@minus,db_wt,nanmean(db_wt(:,base_samps,:,:),2));   % decibel baselined
    % induced
    db_wt = 10*bsxfun(@minus,log10(wt3_ind),log10(nanmean(wt3_ind(:,base_samps,:,:),2)));
    wt_ind_based = bsxfun(@minus,db_wt,nanmean(db_wt(:,base_samps,:,:),2));   % decibel baselined
    % evoked
    db_wt = 10*bsxfun(@minus,log10(wt3_evk),log10(nanmean(wt3_evk(:,base_samps,:,:),2)));
    wt_evk_based = bsxfun(@minus,db_wt,nanmean(db_wt(:,base_samps,:,:),2));   % decibel baselined
    
    avg_wt=nanmean(wt_based,4);
    avg_wt_ind=nanmean(wt_ind_based,4);
    avg_wt_evk=nanmean(wt_evk_based,4);
    %% PLV/PLI calculations based on wavelets
    fprintf('Calculating PLV & PLI ...\n')
    sf=find(F2<=min_max_freq(1)); if isempty(sf); sf=1;end
    ef=find(F2<=min_max_freq(2)); if isempty(ef); ef=length(F2);end
    f_samps=sf(end):ef(end);
    phase_data=angle(wt2(f_samps,:,:,:));
    
    F_plv=F2(f_samps);
    coi_wt2=coi_wt; coi_wt2(coi_wt>max(F2(f_samps)))=nan; coi_wt2(coi_wt<min(F2(f_samps)))=nan;
    clear plv_data pli_data;
    chan_contrasts=nchoosek(1:size(h.current_peak_swf_trials,2),2); surg_flag=0; num_resamps=1;
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
    ss=find(pli_lat<=cfg.study.base_int(1)); if isempty(ss); bs(1)=1; else; bs(1)=ss(end); end
    ss=find(pli_lat<=cfg.study.base_int(2)); bs(2)=ss(end);
    base_samps_pli=bs(1):bs(2);
    pli_based=bsxfun(@minus,pli_data,nanmean(pli_data(:,:,base_samps_pli),3));
    dpli_based=bsxfun(@minus,dpli_data,nanmean(dpli_data(:,:,base_samps_pli),3));
    
 if h.monte_carlo_flag == 0  
    %% %%%%%%%%%%%%%%%%%%% Plotting Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initializing plotting parameters
    % num_iter=1500;  % number of iterations to find PLVs
    pos=round(cfg.source.sig_start*cfg.study.srate)-(cfg.study.lat_sim(1)*cfg.study.srate); % 3sigs x Nfreqs
    num_clmns=num_chans; num_rows=num_freqs;
    min_max=[-100 100]; % time domain scale as percent of baseline    %[-abs(max(max(max(h.current_peak_swf_trials)))) abs(max(max(max(h.current_peak_swf_trials))))]*100;
   min_max2 = [-1 1]; %[-3 3]; %str2num(h.edit_plot_caxis.String); % wavelet color axis scale as percent of baseline
       min_max2_evk = [-10 10]; % [-30  30]; %str2num(h.edit_plot_caxis.String); % wavelet color axis scale as percent of baseline
       min_max2_ind = [-1 1]; % [-12 12]; %str2num(h.edit_plot_caxis.String); % wavelet color axis scale as percent of baseline
    min_max3=[-max(max(max(abs(h.current_peak_swf_trials)))) max(max(max(abs(h.current_peak_swf_trials))))]*110; % wavelet color axis scale as percent of baseline
    plv_caxis=[-.25 .25]; pli_caxis=[-.25 .25]; dpli_caxis=pli_caxis; %[-0.25 0.25];
    %     plv_caxis=[0 1]; pli_caxis=[0 1]; dpli_caxis=[-.5 .5]; %[-0.25 0.25];
    
    mrk_clr=h.src_clr; % lines(length(h.current_3D_peak_idx));
    plv_clr=[.7 0 .9; 1 0 1; 1 .6 0];
    xtik=[-.4:.2:1.2];
    f_size=10; % font size for axis & title
    f_size2=8; % font size for legend
    ln_style={'-' '-' '-'};
    t0=find(cfg.study.lat_sim<=0); t0=t0(end);
    
    %% figure(997): Signal final waves
    figure; set(gcf,'color','w'); clf;
    ax=subplot_axes(4,num_clmns,.06,.05,0,0,0);
    for v=1:num_chans
        %% Time-domain waves
        axes(ax(v)); cla;  hold on; axis on;
        p1=plot(cfg.study.lat_sim,squeeze(h.current_peak_swf_trials(:,v,:))*100,'color',[1 1 1]*.6);
        p2=plot(cfg.study.lat_sim,squeeze(nanmean(h.current_peak_swf_trials(:,v,:),3))*100,'color',mrk_clr(v,:),'linewidth',2);
        plot([0 0],[min_max3],'k--');
        axis([cfg.study.plot_time_int min_max3]);  set(gca,'XTick',xtik);
        title(sprintf('Source %.f',h.current_3D_peak_idx(v)),'Color',mrk_clr(v,:)); set(gca,'Fontsize',f_size); box on;
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
        title(sprintf('Total Power: Source %.f',h.current_3D_peak_idx(v)),'Color',mrk_clr(v,:));
        if v==1; ylabel('Frequency (Hz)'); end
        axes(ax(v+6)); cla;  hold on; axis on;
        surf(cfg.study.lat_sim,F2,squeeze(avg_wt_evk(:,:,v))); view(0,90); shading interp; colormap(jet); axis tight;
        plot3(cfg.study.lat_sim,coi_wt,ones(size(coi_wt)),'color',[1 1 1]*.7,'linewidth',2)
        plot3([0 0],[min_max_freq],[min_max2_evk(2) min_max2_evk(2)],'k--');
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(min_max2_evk); set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Evoked Power: Source %.f',h.current_3D_peak_idx(v)),'Color',mrk_clr(v,:));
        if v==1; ylabel('Frequency (Hz)'); end
        axes(ax(v+9)); cla;  hold on; axis on;
        surf(cfg.study.lat_sim,F2,squeeze(avg_wt_ind(:,:,v))); view(0,90); shading interp; colormap(jet); axis tight;
        plot3(cfg.study.lat_sim,coi_wt,ones(size(coi_wt)),'color',[1 1 1]*.7,'linewidth',2)
        plot3([0 0],[min_max_freq],[min_max2_ind(2) min_max2_ind(2)],'k--');
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(min_max2_ind); set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Induced Power: Source %.f',h.current_3D_peak_idx(v)),'Color',mrk_clr(v,:));
        if v==1; ylabel('Frequency (Hz)'); end
        xlabel('Time (sec');
    end
    ax1=axes('Position',[.84 ax(6).Position(2) .1 ax(12).Position(4)]); axis off; hc=colorbar('peer',ax1,'Location','EastOutside');  ax1.Position(3)=.1; ylabel(hc,'Power (dB re:baseline)'); caxis(min_max2); hc.Label.Position=[2 0 0];
    ax2=axes('Position',[.84 ax(9).Position(2) .1 ax(12).Position(4)]); axis off; hc=colorbar('peer',ax2,'Location','EastOutside');  ax2.Position(3)=.1; ylabel(hc,'Power (dB re:baseline)'); caxis(min_max2_evk); hc.Label.Position=[2 0 0];
    ax3=axes('Position',[.84 ax(12).Position(2) .1 ax(12).Position(4)]); axis off; hc=colorbar('peer',ax3,'Location','EastOutside');  ax3.Position(3)=.1; ylabel(hc,'Power (dB re:baseline)'); caxis(min_max2_ind); hc.Label.Position=[2 0 0];
    
    %% figure(998): PLV & PLI plots
    figure; clf; set(gcf,'color','w');
    ax=subplot_axes(3,num_clmns,.06,.05,0,0,0);
    for vx=1:length(chan_contrasts)
        axes(ax(vx)); cla;  hold on; axis on;
        surf(cfg.study.lat_sim,F_plv,squeeze(plv_based(:,vx,:))); view(0,90); shading interp; colormap(jet);
        %         surf(cfg.study.lat_sim,F_plv,squeeze(plv_data(:,vx,:))); view(0,90); shading interp; colormap(jet);
        plot3(cfg.study.lat_sim,coi_wt2,ones(size(coi_wt2))*plv_caxis(2),'color',[1 1 1]*.7,'linewidth',2);
        plot3([0 0],[min_max_freq],[1 1],'k--');
        title(sprintf('PLV Source %.f vs %.f',h.current_3D_peak_idx(chan_contrasts(vx,:))),'Color',plv_clr(vx,:));
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(plv_caxis); set(gca,'XTick',xtik,'Fontsize',f_size);
        if vx==1; ylabel('Freq (Hz)'); end
        
        
        axes(ax(vx+3)); cla;  hold on; axis on;
        surf(pli_lat,F_plv,squeeze(pli_based(:,vx,:))); view(0,90); shading interp; colormap(jet);
        %         surf(pli_lat,F_plv,squeeze(pli_data(:,vx,:))); view(0,90); shading interp; colormap(jet);
        plot3(cfg.study.lat_sim,coi_wt2,ones(size(coi_wt2))*plv_caxis(2),'color',[1 1 1]*.7,'linewidth',2);
        plot3([0 0],[min_max_freq],[1 1],'k--');
        title(sprintf('PLI Source %.f vs %.f',h.current_3D_peak_idx(chan_contrasts(vx,:))),'Color',plv_clr(vx,:));
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(pli_caxis); set(gca,'XTick',xtik,'Fontsize',f_size);
        if vx==1; ylabel('Freq (Hz)'); end
        
        axes(ax(vx+6)); cla;  hold on; axis on;
        surf(pli_lat,F_plv,squeeze(dpli_based(:,vx,:))); view(0,90); shading interp; colormap(jet);
        %          surf(pli_lat,F_plv,squeeze(dpli_data(:,vx,:))-0.5); view(0,90); shading interp; colormap(jet);
        plot3(cfg.study.lat_sim,coi_wt2,ones(size(coi_wt2))*plv_caxis(2),'color',[1 1 1]*.7,'linewidth',2);
        %         surf(cfg.study.lat_sim,Fcoh,squeeze(nanmean(wcoh,3))); view(0,90); shading interp; colormap(jet); axis tight;
        plot3([0 0],[min_max_freq],[1 1],'k--');
        title(sprintf('dPLI Source %.f vs %.f',h.current_3D_peak_idx(chan_contrasts(vx,:))),'Color',plv_clr(vx,:));
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(dpli_caxis);  set(gca,'XTick',xtik,'Fontsize',f_size);
        xlabel('Time (sec)');
        if vx==1; ylabel('Freq (Hz)'); end
    end
    ax1=axes('Position',[.84 ax(3).Position(2) .1 ax(3).Position(4)]); axis off; hc=colorbar('peer',ax1,'Location','EastOutside');  ax1.Position(3)=.1; ylabel(hc,'PLV'); caxis(plv_caxis); hc.Label.Position=[2 0 0];
    ax2=axes('Position',[.84 ax(6).Position(2) .1 ax(6).Position(4)]); axis off; hc=colorbar('peer',ax2,'Location','EastOutside');  ax2.Position(3)=.1; ylabel(hc,'PLI'); caxis(pli_caxis); hc.Label.Position=[2 0 0];
    ax3=axes('Position',[.84 ax(9).Position(2) .1 ax(9).Position(4)]); axis off; hc=colorbar('peer',ax3,'Location','EastOutside');  ax3.Position(3)=.1; ylabel(hc,'dPLI'); caxis(dpli_caxis); hc.Label.Position=[2 0 0];
    
        %% Figure(999): Modulation Index
    figure; clf; 
    [ax]=subplot_axes(3,4,.05,.05,0,0,0);
    a_idx = [1 2 5 6 9 10; 3 4 7 8 11 12];
    y_lim = [0 .1];
    for a=1:length(h.cfg.source.phase_amp_contrasts)
        % source indices of PAC contrasts
        v_fc = h.cfg.source.phase_amp_contrasts(a,1);  % source index for carrier source
        v_fm = h.cfg.source.phase_amp_contrasts(a,2);  % source index for modulator source
        % carrier and modulator indices for target frequency
        fc_idx = find(h.cfg.source.sig_phase_amp_freq_idx(a,:)>0);  % carrier freq index
        fm_idx = h.cfg.source.sig_phase_amp_freq_idx(a,fc_idx);     % modulator freq index
        % carrier and modulater frequency
        fc = nanmean( squeeze(h.cfg.source.sig_freqs(1,fc_idx,:)));
        fm = nanmean( squeeze(h.cfg.source.sig_freqs(1,fm_idx,:)));
       
        % finding wavelet index for fc and fm
        [~,fc_wt]=min(abs(F2-fc)); % freq index of wavelet that is closest to carrier freq 
        [~,fm_wt]=min(abs(F2-fm)); % freq index of wavelet that is closest to carrier freq 
        % sepearting sig and prepost intervals
        sig_int = [h.cfg.source.sig_start(v_fc,fc_idx) h.cfg.source.sig_start(v_fc,fc_idx)+h.cfg.source.sig_durs(v_fc,fc_idx)];
        ss = floor(h.cfg.study.srate *(sig_int-h.cfg.study.lat_sim(1)))+1; sig_samps = ss(1):ss(2);
        prepost_int1 = [h.cfg.study.lat_sim(1) h.cfg.source.sig_start(v_fc,fc_idx)];
        ss1 = floor(h.cfg.study.srate *(prepost_int1-h.cfg.study.lat_sim(1)))+1; 
        prepost_int2 = [sig_int(2) h.cfg.study.lat_sim(end)];
        ss2 = floor(h.cfg.study.srate *(prepost_int2-h.cfg.study.lat_sim(1)))+1; 
        prepost_samps = [ss1(1):ss1(2) ss2(1):ss2(2)] ;

        % signal PAC
        fc_amp = abs(squeeze(wt2(fc_wt,sig_samps,v_fc,:))); fc_amp = reshape(fc_amp,[numel(fc_amp) 1]); 
        fm_phase = angle(squeeze(wt2(fm_wt,sig_samps,v_fm,:))); fm_phase = reshape(fm_phase,[numel(fm_phase) 1]); 
        nbin=36; phase_bin=linspace(-pi,pi,nbin);
        [MI,distKL,amplP,amplQ,binEdges,binCenters]=modulationIndex(fm_phase,fc_amp,nbin);
        b2=bar(ax(a_idx(1,a)),(phase_bin/(2*pi))*360,squeeze(amplP)); b2.BarWidth=1; b2.FaceColor=h.src_clr(v_fc,:); b2.EdgeColor=[1 1 1]*0;
        ax(a_idx(1,a)).YLim = y_lim; 
        title(ax(a_idx(1,a)),sprintf('Signal: Source %.f modulated by Source %.f ', h.current_3D_peak_idx(h.PAC_source_contrasts(a,:)) ));
        text(ax(a_idx(1,a)),(phase_bin(1)/(2*pi))*360,y_lim(2)*.95,sprintf('Modulation Index = %.4f',MI))
        
        % prepost PAC
        fc_amp = abs(squeeze(wt2(fc_wt,prepost_samps,v_fc,:))); fc_amp = reshape(fc_amp,[numel(fc_amp) 1]);
        fm_phase = angle(squeeze(wt2(fm_wt,prepost_samps,v_fm,:))); fm_phase = reshape(fm_phase,[numel(fm_phase) 1]);
        nbin=36; phase_bin=linspace(-pi,pi,nbin);
        [MI,distKL,amplP,amplQ,binEdges,binCenters]=modulationIndex(fm_phase,fc_amp,nbin);
        b2=bar(ax(a_idx(2,a)),(phase_bin/(2*pi))*360,squeeze(amplP)); b2.BarWidth=1; b2.FaceColor=h.src_clr(v_fc,:); b2.EdgeColor=[1 1 1]*0;
        ax(a_idx(2,a)).YLim = y_lim;
        title(ax(a_idx(2,a)),sprintf('Prepost: Source %.f modulated by Source %.f ',  h.current_3D_peak_idx(h.PAC_source_contrasts(a,:)) ));
        text(ax(a_idx(2,a)),(phase_bin(1)/(2*pi))*360,y_lim(2)*.95,sprintf('Modulation Index = %.4f',MI))

    end
    ax(9).YLabel.String = 'Normalized Amplitude Value'; %ax(a_idx(1,5)).YLabel.String = 'Normalized Gamma-Band Amplitude Value'; ax(a_idx(1,7)).YLabel.String = 'Normalized Gamma-Band Amplitude Value';
    ax(9).XLabel.String = 'Modulator Phase (degrees)';
 else
     % exporting for saving during Monte Carlo
     h.current_avg_wt = avg_wt; h.current_avg_wt_evk = avg_wt_evk; h.current_avg_wt_ind = avg_wt_ind;
     h.current_plv_based = plv_based; h.current_pli_based = pli_based; h.current_dpli_based = dpli_based;
     h.current_plv_data = plv_data; h.current_pli_data = pli_data;
     h.current_TFR_freqs = F2;
     h.current_coi_wt2 = coi_wt2;
     h.current_PLV_freqs = F_plv;
     h.current_pli_lat = pli_lat;
     
 end
    
else
    msgbox(sprintf('\nNo Peak Data\n\nPlease change Image Threshold\n')); 
end
if exist('hm'); close(hm); end
if exist('hm1'); close(hm1); end
if exist('hm2'); close(hm2); end

if h.monte_carlo_flag ~= 1
    h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
end

