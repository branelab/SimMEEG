function plot_SimSignals(varargin)
global h

% hm = msgbox(sprintf('Time-Frequency Analyses of Source Data\n\n Calculating ...'));
h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Time-Frequency Analyses of Source Data\n\n Calculating ...'); drawnow;

TB = str2num(h.edit_wavelet_TB.String); % wavelet parameter --> The larger the time-bandwidth parameter, the more spread out the wavelet is in time and narrower the wavelet is in frequency.


%% Plotting --> Confirming
if ~isempty(h.sim_data)
    sig_final = h.sim_data.sig_final;
    sig_wav = h.sim_data.sig_wav;
    sig_win = h.sim_data.sig_win;
    
    prepost_wav = h.sim_data.prepost_wav;
    prepost_win = h.sim_data.prepost_win;
    
    if isempty(sig_wav); sig_wav = zeros(size(sig_final)); end
    if isempty(sig_win); sig_win = zeros(size(sig_final)); end
    if isempty(prepost_wav); prepost_wav = zeros(size(sig_final)); end
    if isempty(prepost_win); prepost_win = zeros(size(sig_final)); end
    
    % padd zeros for ARM sources
    dims1 = size(sig_wav); dims2 = size(sig_final);
    if any(dims1(1:3)~=dims2(1:3)); vx = dims1(2)+1:dims2(2); sig_wav(:,vx,:) = zeros(dims1(1),length(vx),dims1(3)); end
    dims1 = size(sig_win); dims2 = size(sig_final);
    if any(dims1(1:3)~=dims2(1:3)); vx = dims1(2)+1:dims2(2); sig_win(:,vx,:) = zeros(dims1(1),length(vx),dims1(3)); end
    
    dims1 = size(prepost_wav); dims2 = size(sig_final);
    if any(dims1(1:3)~=dims2(1:3)); vx = dims1(2)+1:dims2(2); prepost_wav(:,vx,:) = zeros(dims1(1),length(vx),dims1(3)); end
    dims1 = size(prepost_win); dims2 = size(sig_final);
    if any(dims1(1:3)~=dims2(1:3)); vx = dims1(2)+1:dims2(2); prepost_win(:,vx,:) = zeros(dims1(1),length(vx),dims1(3)); end
        
    cfg = h.sim_data.cfg;
    cfg.study.plot_time_int = str2num(h.edit_plot_time_int.String);
    cfg.study.plot_freq_int = str2num(h.edit_plot_freq_int.String);
    cfg.study.plot_caxis = str2num(h.edit_plot_caxis.String);
    
    %% %%%%%%%%%%%%%%%%%%% Time-Frequency Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TFR & PLV/PLI parameters
    [num_chans,num_freqs,num_minmaxfr]=size(cfg.source.sig_freqs);
    num_chans = size(sig_final,2);
    lat=cfg.study.lat_sim;
    min_max_freq=cfg.study.plot_freq_int;
    %% calculate wavelets (total & induced) under signal final
    clear wt wt_ind wt_evk;
    sig_ind=bsxfun(@minus,sig_final,nanmean(sig_final,3)); % induced by subtracting mean across trials (i.e., evoked response)
    fprintf('Calculating wavelets ...\n')
    for v=1:num_chans
        %% Wavelets - Total Power
        %         wt_param=[3 30]; %[3 60];
        for t=1:size(sig_final,3)
            %             [wt(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_final(:,v,t)),'morse',cfg.study.srate,'WaveletParameters',wt_param); % total power
            %             [wt_ind(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_ind(:,v,t)),'morse',cfg.study.srate,'WaveletParameters',wt_param);    %induced power
            [wt(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_final(:,v,t)),'morse',cfg.study.srate,'TimeBandwidth',TB); % total power
            [wt_ind(:,:,v,t),F,coi_wt]=cwt(squeeze(sig_ind(:,v,t)),'morse',cfg.study.srate,'TimeBandwidth',TB);    %induced power
        end
        %         [wt_evk(:,:,v),F,coi_wt]=cwt(squeeze(nanmean(sig_final(:,v,:),3)),'morse',cfg.study.srate,'WaveletParameters',wt_param);    % evoked power
        [wt_evk(:,:,v),F,coi_wt]=cwt(squeeze(nanmean(sig_final(:,v,:),3)),'morse',cfg.study.srate,'TimeBandwidth',TB);    % evoked power
    end
    F2=flipud(F); wt2=flipud(wt); wt2_ind=flipud(wt_ind); wt2_evk=flipud(wt_evk);
    ss=find(cfg.study.lat_sim<=cfg.study.base_int(1)); bs(1)=ss(end);
    ss=find(cfg.study.lat_sim<=cfg.study.base_int(2)); bs(2)=ss(end);
    base_samps=bs(1):bs(2);
    wt3=abs(wt2); % converting to real
    wt3_ind=abs(wt2_ind); % converting to real
    wt3_evk=abs(wt2_evk); % converting to real
    %     wt_based=20*log10(bsxfun(@rdivide,wt3,nanmean(wt3(:,base_samps,:),2))); % dB
    % dividing by baseline then multiply 100 to get percent then baseline
    %     wt_based=bsxfun(@rdivide,wt3,nanmean(wt3(:,base_samps,:,:),2))*100;   % percentage
    %     wt_ind_based=bsxfun(@rdivide,wt3_ind,nanmean(wt3_ind(:,base_samps,:,:),2))*100;   % percentage
    %     wt_evk_based=bsxfun(@rdivide,wt3_evk,nanmean(nanmean(wt3(:,base_samps,:),2),3))*100;   % percentage
    %     % baselining
    %     wt_based=bsxfun(@minus,wt_based,nanmean(wt_based(:,base_samps,:,:),2));   % percentage
    %     wt_ind_based=bsxfun(@minus,wt_ind_based,nanmean(wt_ind_based(:,base_samps,:,:),2));   % percentage
    %     wt_evk_based=bsxfun(@minus,wt_evk_based,nanmean(wt_evk_based(:,base_samps,:),2));   % percentage
    
    
    % Decibel
    % all
    db_wt = 10*bsxfun(@minus,log10(wt3),log10(nanmean(wt3(:,h.cfg.study.base_samps,:,:),2)));
    wt_based = bsxfun(@minus,db_wt,nanmean(db_wt(:,h.cfg.study.base_samps,:,:),2));   % decibel baselined
    % induced
    db_wt = 10*bsxfun(@minus,log10(wt3_ind),log10(nanmean(wt3_ind(:,h.cfg.study.base_samps,:,:),2)));
    wt_ind_based = bsxfun(@minus,db_wt,nanmean(db_wt(:,h.cfg.study.base_samps,:,:),2));   % decibel baselined
    % evoked
    db_wt = 10*bsxfun(@minus,log10(wt3_evk),log10(nanmean(wt3_evk(:,h.cfg.study.base_samps,:,:),2)));
    wt_evk_based = bsxfun(@minus,db_wt,nanmean(db_wt(:,h.cfg.study.base_samps,:,:),2));   % decibel baselined
    
    % averaging across trials
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
    chan_contrasts=nchoosek(1:size(sig_final,2),2); surg_flag=0; num_resamps=1;
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
    
    %% %%%%%%%%%%%%%%%%%%% Plotting Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initializing plotting parameters
    % num_iter=1500;  % number of iterations to find PLVs
    pos=round(cfg.source.sig_start*cfg.study.srate)-(cfg.study.lat_sim(1)*cfg.study.srate); % 3sigs x Nfreqs
    num_clmns=num_chans; num_rows=num_freqs;
    min_max = str2num(h.edit_plot_caxis.String); % time domain scale as percent of baseline    %[-abs(max(max(max(sig_final)))) abs(max(max(max(sig_final))))]*100;
    min_max2 = [-3 3]; %str2num(h.edit_plot_caxis.String); % wavelet color axis scale as percent of baseline
    min_max2_evk = [-30  30]; %str2num(h.edit_plot_caxis.String); % wavelet color axis scale as percent of baseline
    min_max2_ind = [-12 12]; %str2num(h.edit_plot_caxis.String); % wavelet color axis scale as percent of baseline
    min_max3=[-max(max(max(abs(sig_final)))) max(max(max(abs(sig_final))))]*110; % wavelet color axis scale as percent of baseline
    %     plv_caxis=[-.5 .5]; pli_caxis=[-.5 .5]; dpli_caxis=pli_caxis; %[-0.25 0.25];
    min_max2 = [-3 3]; %str2num(h.edit_plot_caxis.String); % wavelet color axis scale as percent of baseline
    min_max2_evk = [-15  15]; %str2num(h.edit_plot_caxis.String); % wavelet color axis scale as percent of baseline
    min_max2_ind = [-3 3]; %str2num(h.edit_plot_caxis.String); % wavelet color axis scale as percent of baseline
    plv_caxis=[-.25 .25]; pli_caxis=[-.25 .25]; dpli_caxis=pli_caxis; %[-0.25 0.25];
    
    mrk_clr=h.cfg.source.src_clr;
    plv_clr=[.7 0 .9; 1 0 1; 1 .6 0]; 
    xtik=[-.4:.2:1.2];
    f_size=10; % font size for axis & title
    f_size2=8; % font size for legend
    ln_style={'-' '-' '-'};
    t0=find(cfg.study.lat_sim<=0); t0=t0(end);
    
    if size(find(sig_wav~=0),1)~=0  % then synthetic signal was generated, else real source data were loaded
        %% figure(995): Signal & Prepost example waves separate plots
        figure(995); clf; set(gcf,'color','w');
        ax=subplot_axes(num_rows+1,num_clmns,.06,.06,0,0,0);
        v=1; tx=[1 2]; % trial# for example waves
        %     min_max4=[-max(abs(squeeze(nansum(sig_win(:,v,tx,:),4)+nansum(prepost_win(:,v,tx,:),4)))) max(abs(squeeze(nansum(sig_win(:,v,tx,:),4)+nansum(prepost_win(:,v,tx,:),4))))]*110;
        min_max4=[-max(max(abs(squeeze(nansum(sig_win(:,v,tx,:),4)+nansum(prepost_win(:,v,tx,:),4))))) max(max(abs(squeeze(nansum(sig_win(:,v,tx,:),4)+nansum(prepost_win(:,v,tx,:),4)))))]*110;
        if min_max4(1)==min_max4(2); min_max4=[-110 110]; end
        a=0;
        
        
        for f=1:num_freqs
            a=a+1;
            axes(ax(a)); cla; hold on; axis on;
            plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
            p1=plot(cfg.study.lat_sim,squeeze(sig_wav(:,v,tx,f))*100,'color',mrk_clr(v,:),'linewidth',1); %(2).LineWidth=2;
            p2=plot(cfg.study.lat_sim,squeeze(sig_win(:,v,tx,f))*100,'-','color',mrk_clr(v,:),'linewidth',2);
            axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
            title(sprintf('Source %.f Signal(%.1f to %.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:)); box on;
            a=a+1;
            axes(ax(a)); cla; hold on; axis on; box on;
            plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
            p1=plot(cfg.study.lat_sim,squeeze(prepost_wav(:,v,tx,f))*100,'color',mrk_clr(v,:),'linewidth',1); %p1(2).LineWidth=2;
            p2=plot(cfg.study.lat_sim,squeeze(prepost_win(:,v,tx,f))*100,'-','color',mrk_clr(v,:),'linewidth',2);
            axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
            title(sprintf('Source %.f Prepost (%.1f to %.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:)); box on;
            
            a=a+1;
            axes(ax(a)); cla; hold on; axis on; box on;
            plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
            p1=plot(cfg.study.lat_sim,squeeze(sig_wav(:,v,tx,f)+prepost_wav(:,v,tx,f))*100,'color',mrk_clr(v,:),'linewidth',1); %p1(2).LineWidth=2;
            p2=plot(cfg.study.lat_sim,squeeze(sig_win(:,v,tx,f)+prepost_win(:,v,tx,f))*100,'-','color',mrk_clr(v,:),'linewidth',2);
            axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
            title(sprintf('Source %.f (%.1f to %.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:)); box on;
        end
        %% final signal
        a=a+1;
        axes(ax(a)); cla; hold on; axis on;
        plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
        p1=plot(cfg.study.lat_sim,squeeze(nansum(sig_wav(:,v,tx,:),4))*100,'color',mrk_clr(v,:),'linewidth',1); %p1(2).LineWidth=2;
        p2=plot(cfg.study.lat_sim,squeeze(nansum(sig_win(:,v,tx,:),4))*100,'-','color',mrk_clr(v,:),'linewidth',2);
        axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Source %.f Signal Sum',v),'Color',mrk_clr(v,:)); box on;
        a=a+1;
        axes(ax(a)); cla; hold on; axis on; box on;
        plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
        p1=plot(cfg.study.lat_sim,squeeze(nansum(prepost_wav(:,v,tx,:),4))*100,'color',mrk_clr(v,:),'linewidth',1); %p1(2).LineWidth=2;
        p2=plot(cfg.study.lat_sim,squeeze(nansum(prepost_win(:,v,tx,:),4))*100,'-','color',mrk_clr(v,:),'linewidth',2);
        axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Source %.f PrePost Sum',v),'Color',mrk_clr(v,:)); box on;
        
        a=a+1;
        axes(ax(a)); cla; hold on; axis on; box on;
        plot([0 0],min_max4,'k--','linewidth',1); plot(cfg.study.plot_time_int,[0 0],'k-','linewidth',1);
        p1=plot(cfg.study.lat_sim,squeeze(nansum(sig_wav(:,v,tx,:),4)+nansum(prepost_wav(:,v,tx,:),4))*100,'color',mrk_clr(v,:),'linewidth',1); %p1(2).LineWidth=2;
        p2=plot(cfg.study.lat_sim,squeeze(nansum(sig_win(:,v,tx,:),4)+nansum(prepost_win(:,v,tx,:),4))*100,'-','color',mrk_clr(v,:),'linewidth',2);
        axis([cfg.study.plot_time_int min_max4]);  set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Source %.f Sum',v),'Color',mrk_clr(v,:)); box on;
        
        %% figure(996): Signal & Prepost waves overlaid
        figure(996); clf; set(gcf,'color','w');
        ax=subplot_axes(num_rows,num_clmns,.06,.05,0,0,0);
        a=0;
        for f=1:num_freqs
            for v=1:num_chans
                a=a+1;
                axes(ax(a)); cla; hold on; axis on;
                p1=plot(cfg.study.lat_sim,squeeze(prepost_wav(:,v,:,f))*100,'color',[1 1 1]*.6);
                p2=plot(cfg.study.lat_sim,squeeze(sig_wav(:,v,:,f))*100,'color',mrk_clr(v,:));
                p3=plot(cfg.study.lat_sim,squeeze(nanmean(prepost_wav(:,v,:,f),3))*100,'color',[1 1 1]*.4,'linewidth',2);
                p4=plot(cfg.study.lat_sim,squeeze(nanmean(sig_wav(:,v,:,f),3))*100,'color',mrk_clr(v,:)*.75,'linewidth',2);
                plot([0 0],[min_max],'k--');
                axis([cfg.study.plot_time_int min_max]);  set(gca,'XTick',xtik,'Fontsize',f_size);
                try title(sprintf('Source %.f (%.1f-%.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:)); box on;
                catch
                    title(sprintf('Source %.f',v),'Color',mrk_clr(v,:)); box on;
                end
                if v==1; ylabel('Amplitude (%)'); end
                if f==num_freqs; xlabel('Time (sec'); end
                
                if f==1
                    legend([p2(1), p4(1), p1(1),  p3(1) ], {'Signal','Avg Signal','Prepost','Avg Prepost'},'Location','NorthWest','Fontsize',f_size2)
                end
            end
        end
    end
    %% figure(997): Signal final waves
    figure(997); set(gcf,'color','w'); clf;
    ax=subplot_axes(4,num_clmns,.06,.05,0,0,0);
    for v=1:num_chans
        %% Time-domain waves
        axes(ax(v)); cla;  hold on; axis on;
        p1=plot(cfg.study.lat_sim,squeeze(sig_final(:,v,:))*100,'color',[1 1 1]*.6);
        p2=plot(cfg.study.lat_sim,squeeze(nanmean(sig_final(:,v,:),3))*100,'color',mrk_clr(v,:),'linewidth',2);
        plot([0 0],[min_max3],'k--');
        axis([cfg.study.plot_time_int min_max3]);  set(gca,'XTick',xtik);
        title(sprintf('Source %.f',v),'Color',mrk_clr(v,:)); set(gca,'Fontsize',f_size); box on;
        legend([p1(1),p2],{'Trials','Average'},'Location','NorthWest','FontSize',f_size2)
        if v==1; ylabel('Amplitude (%)'); end
    end
    %% Power wavelets
    for v=1:num_chans
        axes(ax(v+num_chans)); cla;  hold on; axis on;
        surf(cfg.study.lat_sim,F2,squeeze(avg_wt(:,:,v))); view(0,90); shading interp; colormap(jet); axis tight;
        plot3(cfg.study.lat_sim,coi_wt,ones(size(coi_wt)),'color',[1 1 1]*.7,'linewidth',2)
        plot3([0 0],[min_max_freq],[min_max2(2) min_max2(2)],'k--');
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(min_max2); set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Total Power: Source %.f',v),'Color',mrk_clr(v,:));
        if v==1; ylabel('Frequency (Hz)'); end
        axes(ax(v+(2*num_chans))); cla;  hold on; axis on;
        surf(cfg.study.lat_sim,F2,squeeze(avg_wt_evk(:,:,v))); view(0,90); shading interp; colormap(jet); axis tight;
        plot3(cfg.study.lat_sim,coi_wt,ones(size(coi_wt)),'color',[1 1 1]*.7,'linewidth',2)
        plot3([0 0],[min_max_freq],[min_max2_evk(2) min_max2_evk(2)],'k--');
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(min_max2_evk); set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Evoked Power: Source %.f',v),'Color',mrk_clr(v,:));
        if v==1; ylabel('Frequency (Hz)'); end
        axes(ax(v+(3*num_chans))); cla;  hold on; axis on;
        surf(cfg.study.lat_sim,F2,squeeze(avg_wt_ind(:,:,v))); view(0,90); shading interp; colormap(jet); axis tight;
        plot3(cfg.study.lat_sim,coi_wt,ones(size(coi_wt)),'color',[1 1 1]*.7,'linewidth',2)
        plot3([0 0],[min_max_freq],[min_max2_ind(2) min_max2_ind(2)],'k--');
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(min_max2_ind); set(gca,'XTick',xtik,'Fontsize',f_size);
        title(sprintf('Induced Power: Source %.f',v),'Color',mrk_clr(v,:));
        if v==1; ylabel('Frequency (Hz)'); end
        xlabel('Time (sec');
    end
    ax1=axes('Position',[.84 ax(6).Position(2) .1 ax(12).Position(4)]); axis off; hc=colorbar('peer',ax1,'Location','EastOutside');  ax1.Position(3)=.1; ylabel(hc,'Power (dB re:baseline)'); caxis(min_max2); hc.Label.Position=[2 0 0];
    ax2=axes('Position',[.84 ax(9).Position(2) .1 ax(12).Position(4)]); axis off; hc=colorbar('peer',ax2,'Location','EastOutside');  ax2.Position(3)=.1; ylabel(hc,'Power (dB re:baseline)'); caxis(min_max2_evk); hc.Label.Position=[2 0 0];
    ax3=axes('Position',[.84 ax(12).Position(2) .1 ax(12).Position(4)]); axis off; hc=colorbar('peer',ax3,'Location','EastOutside');  ax3.Position(3)=.1; ylabel(hc,'Power (dB re:baseline)'); caxis(min_max2_ind); hc.Label.Position=[2 0 0];
    
    %% figure(998): PLV & PLI plots
     plv_clr = repmat([0 0 0],length(chan_contrasts),1);
    xn = ceil(length(chan_contrasts)^.5);
    figure(9981); clf; set(gcf,'color','w');    % PLV
    ax1=subplot_axes(xn,xn,.06,.05,0,0,0);
    figure(9982); clf; set(gcf,'color','w');    % PLI
    ax2=subplot_axes(xn,xn,.06,.05,0,0,0);
    figure(9983); clf; set(gcf,'color','w');    % dPLI
     ax3=subplot_axes(xn,xn,.06,.05,0,0,0);
   for vx=1:length(chan_contrasts)
        axes(ax1(vx)); cla;  hold on; axis on;
        surf(cfg.study.lat_sim,F_plv,squeeze(plv_based(:,vx,:))); view(0,90); shading interp; colormap(jet);
        %         surf(cfg.study.lat_sim,F_plv,squeeze(plv_data(:,vx,:))); view(0,90); shading interp; colormap(jet);
        plot3(cfg.study.lat_sim,coi_wt2,ones(size(coi_wt2)),'color',[1 1 1]*.7,'linewidth',2);
        plot3([0 0],[min_max_freq],[1 1],'k--');
        title(sprintf('PLV Source %.f vs %.f',chan_contrasts(vx,:)),'Color',plv_clr(vx,:));
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(plv_caxis); set(gca,'XTick',xtik,'Fontsize',f_size);
        if vx==1; ylabel('Freq (Hz)'); end
        
        axes(ax2(vx)); cla;  hold on; axis on;
        surf(pli_lat,F_plv,squeeze(pli_based(:,vx,:))); view(0,90); shading interp; colormap(jet);
        %         surf(pli_lat,F_plv,squeeze(pli_data(:,vx,:))); view(0,90); shading interp; colormap(jet);
        plot3(cfg.study.lat_sim,coi_wt2,ones(size(coi_wt2)),'color',[1 1 1]*.7,'linewidth',2);
        plot3([0 0],[min_max_freq],[1 1],'k--');
        title(sprintf('PLI Source %.f vs %.f',chan_contrasts(vx,:)),'Color',plv_clr(vx,:));
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(pli_caxis); set(gca,'XTick',xtik,'Fontsize',f_size);
        if vx==1; ylabel('Freq (Hz)'); end
        
        axes(ax3(vx)); cla;  hold on; axis on;
        surf(pli_lat,F_plv,squeeze(dpli_based(:,vx,:))); view(0,90); shading interp; colormap(jet);
        %          surf(pli_lat,F_plv,squeeze(dpli_data(:,vx,:))-0.5); view(0,90); shading interp; colormap(jet);
        plot3(cfg.study.lat_sim,coi_wt2,ones(size(coi_wt2)),'color',[1 1 1]*.7,'linewidth',2);
        %         surf(cfg.study.lat_sim,Fcoh,squeeze(nanmean(wcoh,3))); view(0,90); shading interp; colormap(jet); axis tight;
        plot3([0 0],[min_max_freq],[1 1],'k--');
        title(sprintf('dPLI Source %.f vs %.f',chan_contrasts(vx,:)),'Color',plv_clr(vx,:));
        axis([cfg.study.plot_time_int cfg.study.plot_freq_int]); caxis(dpli_caxis);  set(gca,'XTick',xtik,'Fontsize',f_size);
        xlabel('Time (sec)');
        if vx==1; ylabel('Freq (Hz)'); end
    end
    ax1=axes('Position',[.84 ax(3).Position(2) .1 ax(3).Position(4)]); axis off; hc=colorbar('peer',ax1,'Location','EastOutside');  ax1.Position(3)=.1; ylabel(hc,'PLV'); caxis(plv_caxis); hc.Label.Position=[2 0 0];
    ax2=axes('Position',[.84 ax(6).Position(2) .1 ax(6).Position(4)]); axis off; hc=colorbar('peer',ax2,'Location','EastOutside');  ax2.Position(3)=.1; ylabel(hc,'PLI'); caxis(pli_caxis); hc.Label.Position=[2 0 0];
    ax3=axes('Position',[.84 ax(9).Position(2) .1 ax(9).Position(4)]); axis off; hc=colorbar('peer',ax3,'Location','EastOutside');  ax3.Position(3)=.1; ylabel(hc,'dPLI'); caxis(dpli_caxis); hc.Label.Position=[2 0 0];
    
    %     %% figure(1000) & figure(1001): Polar plot of phases
    %     bin_wdth=20;
    %     % get histcounts to set maximum polarhistogram values.
    %     for f=1:num_freqs
    %         for v=1:num_chans
    %             ph=histcounts(squeeze(cfg.prepost_phase(v,:,f)),bin_wdth);
    %             p1_max(v,f,1)=max(ph);
    %             ph=histcounts(squeeze(cfg.sig_phase(v,:,f)),bin_wdth);
    %             p1_max(v,f,2)=max(ph);
    %
    %             pd1=cfg.prepost_phase(chan_contrasts(v,1),:,f)-cfg.prepost_phase(chan_contrasts(v,2),:,f);
    %             pd2=cfg.sig_phase(chan_contrasts(v,1),:,f)-cfg.sig_phase(chan_contrasts(v,2),:,f);
    %             ph=histcounts(pd1,bin_wdth);
    %             pd_max(v,f,1)=max(ph);
    %             ph=histcounts(pd2,bin_wdth);
    %             pd_max(v,f,2)=max(ph);
    %         end
    %     end
    %     bin_axis=[0 360 0 max(max(max(p1_max)))+1]; bin_axis2=[0 360 0 max(max(max(pd_max)))+1];
    %     ab=0; ac=0;
    %     figure(1000); clf; set(gcf,'color','w');
    %     ax1=subplot_axes(num_rows,6,0.05,0.05,0,0.05,0); ax1_pos={ax1.Position}; clf;
    %     figure(1001); clf; set(gcf,'color','w');
    %     ax2=subplot_axes(num_rows,6,0.05,0.05,0,0,0);  ax2_pos={ax2.Position}; clf;
    %
    %     for f=1:num_freqs
    %         % Signal phases
    %         figure(1000);
    %         for v=1:num_chans
    %             ab=ab+1;
    %             axes('Position',ax1_pos{ab}); cla;
    %             %         polarscatter(squeeze(cfg.prepost_phase(v,:,f)),ones(size(squeeze(cfg.prepost_phase(v,:,f)))),'markerfacecolor',mrk_clr(v,:),'markeredgecolor',mrk_clr(v,:)*.5);
    %             polarhistogram(squeeze(cfg.prepost_phase(v,:,f)),bin_wdth,'facecolor',mrk_clr(v,:),'edgecolor',mrk_clr(v,:)*.5)
    %             title(sprintf('PrePost Phases\nSource %.f\n(%.1f-%.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:),'Fontsize',f_size);
    %             axis(bin_axis);
    %             axes('Position',ax1_pos{ab+3}); cla;
    %             %         polarscatter(squeeze(cfg.sig_phase(v,:,f)),ones(size(squeeze(cfg.sig_phase(v,:,f)))),'markerfacecolor',mrk_clr(v,:),'markeredgecolor',mrk_clr(v,:)*.5);
    %             polarhistogram(squeeze(cfg.sig_phase(v,:,f)),bin_wdth,'facecolor',mrk_clr(v,:),'edgecolor',mrk_clr(v,:)*.5)
    %             title(sprintf('Signal Phases\nSource %.f\n(%.1f-%.1f Hz)',v,squeeze(cfg.source.sig_freqs(v,f,:))),'Color',mrk_clr(v,:),'Fontsize',f_size);
    %             axis(bin_axis);
    %         end
    %         ab=ab+3;
    %         % Signal phase difference for PLV
    %         figure(1001);
    %         for vx=1:num_chans
    %             ac=ac+1;
    %             axes('Position',ax2_pos{ac}); cla;
    %             pd=cfg.prepost_phase(chan_contrasts(vx,1),:,f)-cfg.prepost_phase(chan_contrasts(vx,2),:,f);
    %             polarhistogram(pd,bin_wdth,'facecolor',plv_clr(vx,:),'edgecolor',plv_clr(vx,:)*.5)
    %             %         polarscatter(pd,ones(size(pd)),'markerfacecolor',plv_clr(vx,:),'markeredgecolor',plv_clr(vx,:)*.5);
    %             title(sprintf('PrePost Phase Diff\nSource %.f-%.f \n(%.1f-%.1f Hz)',chan_contrasts(vx,:),squeeze(cfg.source.sig_freqs(v,f,:))),'Color',plv_clr(vx,:),'Fontsize',f_size);
    %             plv12=abs((sum(exp(1i*(squeeze(pd))),2))')./size(cfg.sig_phase,2);
    %             my_sine = round(sin(pd)*(10^6))/(10^6); sign_test = sign(my_sine); pli12 = squeeze(abs(mean(squeeze(nanmean(sign_test, 2)), 3)));
    %             Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; dpli12 = 2*(squeeze(mean(mean(Y, 2), 4))-.5);
    %             text((250/360)*2*pi,max(max(max(pd_max)))*2.1,sprintf('PLV = %.2f\nPLI = %.2f\ndPLI = %.2f',plv12,pli12,dpli12),'Color',plv_clr(vx,:));
    %             axis(bin_axis2);
    %
    %             axes('Position',ax2_pos{ac+3}); cla;
    %             pd=cfg.sig_phase(chan_contrasts(vx,1),:,f)-cfg.sig_phase(chan_contrasts(vx,2),:,f);
    %             polarhistogram(pd,bin_wdth,'facecolor',plv_clr(vx,:),'edgecolor',plv_clr(vx,:)*.5)
    %             %       polarscatter(pd,ones(size(pd)),'markerfacecolor',plv_clr(vx,:),'markeredgecolor',plv_clr(vx,:)*.5);
    %             title(sprintf('Signal Phase Diff\n Source %.f-%.f \n(%.1f-%.1f Hz)',chan_contrasts(vx,:),squeeze(cfg.source.sig_freqs(v,f,:))),'Color',plv_clr(vx,:),'Fontsize',f_size);
    %             plv12=abs((sum(exp(1i*(squeeze(pd))),2))')./size(cfg.sig_phase,2);
    %             my_sine = round(sin(pd)*(10^6))/(10^6); sign_test = sign(my_sine); pli12 = squeeze(abs(mean(squeeze(nanmean(sign_test, 2)), 3)));
    %             Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5; dpli12 = 2*(squeeze(mean(mean(Y, 2), 4))-.5);
    %             text((250/360)*2*pi,max(max(max(pd_max)))*2.1,sprintf('PLV = %.2f\nPLI = %.2f\ndPLI = %.2f',plv12,pli12,dpli12),'Color',plv_clr(vx,:));
    %             axis(bin_axis2);
    %         end
    %         ac=ac+3;
    %
    %     end
    
    %% Figure(999): Modulation Index
    figure(999); clf;
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
        title(ax(a_idx(1,a)),sprintf('Signal: Source %.f modulated by Source %.f ', h.PAC_source_contrasts(a,:)));
        text(ax(a_idx(1,a)),(phase_bin(1)/(2*pi))*360,y_lim(2)*.95,sprintf('Modulation Index = %.4f',MI))
        
        % prepost PAC
        fc_amp = abs(squeeze(wt2(fc_wt,prepost_samps,v_fc,:))); fc_amp = reshape(fc_amp,[numel(fc_amp) 1]);
        fm_phase = angle(squeeze(wt2(fm_wt,prepost_samps,v_fm,:))); fm_phase = reshape(fm_phase,[numel(fm_phase) 1]);
        nbin=36; phase_bin=linspace(-pi,pi,nbin);
        [MI,distKL,amplP,amplQ,binEdges,binCenters]=modulationIndex(fm_phase,fc_amp,nbin);
        b2=bar(ax(a_idx(2,a)),(phase_bin/(2*pi))*360,squeeze(amplP)); b2.BarWidth=1; b2.FaceColor=h.src_clr(v_fc,:); b2.EdgeColor=[1 1 1]*0;
        ax(a_idx(2,a)).YLim = y_lim;
        title(ax(a_idx(2,a)),sprintf('Prepost: Source %.f modulated by Source %.f ', h.PAC_source_contrasts(a,:)));
        text(ax(a_idx(2,a)),(phase_bin(1)/(2*pi))*360,y_lim(2)*.95,sprintf('Modulation Index = %.4f',MI))
        
    end
    ax(9).YLabel.String = 'Normalized Amplitude Value'; %ax(a_idx(1,5)).YLabel.String = 'Normalized Gamma-Band Amplitude Value'; ax(a_idx(1,7)).YLabel.String = 'Normalized Gamma-Band Amplitude Value';
    ax(9).XLabel.String = 'Modulator Phase (degrees)';
    
else
    msgbox(sprintf('\nNo Simulated Data exists\n\nPlease Run Simulation\n'));
end
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
% close(hm)
