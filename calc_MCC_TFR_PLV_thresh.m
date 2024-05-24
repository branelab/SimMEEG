function [metrics] = calc_MCC_TFR_PLV_thresh(FC_true, FC_seed, FC_comp, plv_comp_indices, plv_comp_idx, plv_true_locs, plv_seed_locs, plv_comp_locs, vol, freq_int, time_int, ROI_freq_int, ROI_time_int, thresh_vals, freqs, lat, plot_flag, fc_scale, fc_map_scale, ROI_only_flag)
% This calculates MCC for all the sample points and all frequencies
% within the FC data scanning through threshold plv/pli values "thresh_vals"
%
% INPUT:
%   FC_true = true FC [freqs x comparisons x samps]
%   FC_seed = observed FC for the seed FCs (only include Hits) with same order as FC true [freqs x comparisons x samps]
%   FC_comp = observed FC for other comparisons (Correct Rejections) [freqs x comparisons x samps]
%   plv_comp_indices = indices of plv_comp_locs that correspond to the plv_comp_idx; 
%   plv_comp_idx = index of FC comparisons within the plv_comp_locs
%   thresh_vals = (t) FC thresholds for calculating MCC. This will be for negative and positive values for searching for hits, misses, CRs, FAs.
%   freq_int = freq of interest interval for summing across to get a global metrics across this interval.
%   time_int = time of interest interval for summing across to get a global metrics across this interval.
%   lat = latency of samps
%   plot_flag = (1) plots some results (0) no plots
%   ROI_only_flag = (0) calculate stats across all freq_int and time_int and calculate for ROI (1) calculate only for ROI
%
% OUTPUT:
%    metrics.perf(t).MCC                = Mathew's Correlation Coefficient
%                   .true_positives     = True positives (Hits) for each thresh_vals
%                   .false_positives    = False positives (False Alarms) for each thresh_vals
%                   .true_negatives     = True positives (Correct Rejections) for each thresh_vals
%                   .false_negatives    = True positives (Misses) for each thresh_vals
%

% thresh_vals = 0:.1:1;
% FC_data = inv_soln(a).TFR_results.plv_based;
% FC_true = cfg.source.TFR_results.plv_based;
num_seeds = size(plv_seed_locs,1);
plv_true_idx = nchoosek(1:size(plv_true_locs,1),2);

if num_seeds>1
plv_seed_idx = nchoosek(1:num_seeds,2);
else
plv_seed_idx = [];    
end

% find and maintain order of seeds within FC_comp
plv_comp_seed_idx = [];
for v=1:size(plv_seed_idx,1)
    if ~isempty(find(ismember(plv_comp_idx,plv_seed_idx(v,:),'rows')))
    plv_comp_seed_idx(v) = find(ismember(plv_comp_idx,plv_seed_idx(v,:),'rows'));
    end
end

%% freq/time interval of interest for calculating metrics
ss = find(freqs>=freq_int(1)); fx(1) = ss(1); ss = find(freqs<=freq_int(2)); fx(2) = ss(end); fsamps = fx(1):fx(2); 
ss = find(lat>=time_int(1)); sx(1) = ss(1); ss = find(lat<=time_int(2)); sx(2) = ss(end); tsamps = sx(1):sx(2); 
%% ROI freq/time interval of interest for summing across to calculate metrics
ss = find(freqs>=ROI_freq_int(1)); rfx(1) = ss(1); ss = find(freqs<=ROI_freq_int(2)); rfx(2) = ss(end); ROI_fsamps = rfx(1):rfx(2); 
ss = find(lat>=ROI_time_int(1)); rsx(1) = ss(1); ss = find(lat<=ROI_time_int(2)); rsx(2) = ss(end); ROI_tsamps = rsx(1):rsx(2); 


if plot_flag==1
    figure(9999); clf; set(gcf,'color','w');
    ax = subplot_axes(2,5,.035,.125,0,0,0);
end

% thresh_vals = fliplr(thresh_vals);
for t=1:length(thresh_vals)
    %% Thresholding FC data
    %% Seed Hits --> when true>thresh and seed>thresh or true<-thresh and seed<-thresh
    xhit = (FC_true>thresh_vals(t) & FC_seed>thresh_vals(t)) | (FC_true<-thresh_vals(t) & FC_seed<-thresh_vals(t) );
    FC_hit = double(xhit>0);  % Hits for positive & negative FCs
    hit = FC_hit; %hit(FC_hit==0)=nan;
%     true_total = sum(sum(sum(FC_true>thresh_vals(t)))); % total number of FC_true that exceeds threshold
    
    %    v=2; figure(2); clf; surf(squeeze(hit(:,v,:))); view(0,90); shading flat; colormap(jet); axis tight; colorbar;

    %% Seed Misses --> when true>thresh and seed<thresh or true<-thresh and seed>-thresh
    x = (FC_true>=thresh_vals(t) & FC_seed<=thresh_vals(t) & FC_seed>=-thresh_vals(t)) | (FC_true<=-thresh_vals(t) & FC_seed>=-thresh_vals(t) & FC_seed<=thresh_vals(t));
    FC_miss = double(x>0);  % Hits for positive & negative FCs
    miss = FC_miss; %miss(FC_miss==0)=nan;
    %     figure(3); clf;  surf(squeeze(miss(:,v,:))); view(0,90); shading flat; colormap(jet); axis tight; colorbar;

    %% Seed False Alarms--> when true<thresh and seed>thresh or true>-thresh and seed<-thresh
    xfa = (FC_true<=thresh_vals(t) & FC_true>=-thresh_vals(t) & FC_seed>thresh_vals(t)) | (FC_true<=thresh_vals(t) & FC_true>=-thresh_vals(t) & FC_seed<-thresh_vals(t));
    FC_fa = double(xfa>0);  % Hits for positive & negative FCs
    fa = FC_fa; %fa(FC_fa==0)=nan;
    %     figure(4); clf;  surf(squeeze(fa(:,v,:))); view(0,90); shading flat; colormap(jet); axis tight; colorbar;
    
    
    %% Seed Correct Rejections --> when true<thresh and seed<thresh or true>-thresh and seed>-thresh
%     x = (FC_true<=thresh_vals(t) & FC_seed<=thresh_vals(t) & FC_true>=-thresh_vals(t) & FC_seed>=-thresh_vals(t) ) | ( xhit==0 & xfa==0 ); % if no FC_true exists and no FAs then they samples are CRs
%     FC_cr = double(x>0);  % Hits for positive & negative FCs
%     cr = FC_cr; %cr(FC_cr==0)=nan;
    FC_cr = ones(size(FC_seed))-(hit + miss + fa);  % all possible minus hit+miss+fa
    cr = FC_cr; %cr(FC_cr==0)=nan;
    %     figure(5); clf;  surf(squeeze(cr(:,v,:))); view(0,90); shading flat; colormap(jet); axis tight; colorbar;



    %% Comparison False Alarms--> when true<thresh and seed>thresh or true>-thresh and seed<-thresh
    x = FC_comp>=thresh_vals(t) | FC_comp<=-thresh_vals(t);
    comp_fa = double(x>0);  % Hits for positive & negative FCs
    %    figure(6); clf; surf(squeeze(comp_fa(:,v,:))); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis([0 1]);
   
    %% Comparison Correct Rejections --> when true<thresh and seed<thresh or true>-thresh and seed>-thresh
%     x = FC_comp<=thresh_vals(t) & FC_comp>=-thresh_vals(t);
%     comp_cr = double(x>0);  % Hits for positive & negative FCs
    comp_cr = ones(size(FC_comp)) - comp_fa; 
    %    figure(7); clf; surf(squeeze(comp_cr(:,v,:))); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis([0 1]);

    %% summing up data
    hit_sum = squeeze(nansum(hit(fsamps,:,tsamps),2));
    miss_sum = squeeze(nansum(miss(fsamps,:,tsamps),2));
    cr_sum = squeeze(nansum(cr(fsamps,:,tsamps),2) + nansum(comp_cr(fsamps,:,tsamps),2));
%     fa_sum = squeeze(nansum(fa(fsamps,:,tsamps),2) + nansum(comp_fa(fsamps,:,tsamps),2));
    fa_sum = squeeze(nansum(fa(fsamps,:,tsamps),2) + nansum(comp_fa(fsamps,:,tsamps),2));

    % Plot Testing
    %        figure(2); clf; surf(lat(tsamps),freqs(fsamps),squeeze(hit_sum)); view(0,90); shading flat; colormap(jet); axis tight; colorbar;
    %        figure(2); clf; surf(lat(tsamps),freqs(fsamps),squeeze(miss_sum)); view(0,90); shading flat; colormap(jet); axis tight; colorbar;
    %        figure(3); clf; surf(lat(tsamps),freqs(fsamps),squeeze(cr_sum)); view(0,90); shading flat; colormap(jet); axis tight; colorbar;
    %        figure(4); clf; surf(lat(tsamps),freqs(fsamps),squeeze(fa_sum)); view(0,90); shading flat; colormap(jet); axis tight; colorbar;


    %% MCC Performance for FC
    if ROI_only_flag == 0
        [metrics.perf(t)]=calc_classifier_performance(hit_sum,miss_sum,cr_sum,fa_sum);
        metrics.struct_xx_sum = '[freq x samps x threshold]';
        metrics.hit(:,:,t) = hit_sum;
        metrics.miss(:,:,t) = miss_sum;
        metrics.cr(:,:,t) = cr_sum;
        metrics.fa(:,:,t) = fa_sum;
    end
    
    %% freq/time interval of interest for summing to calculate metrics
    hit_sum2 = squeeze( sum( sum( sum( hit(ROI_fsamps,:,ROI_tsamps) ) ) ) );
    miss_sum2 = squeeze( sum( sum( sum( miss(ROI_fsamps,:,ROI_tsamps) ) ) ) );
    cr_sum2 = squeeze( sum( sum( sum( cr(ROI_fsamps,:,ROI_tsamps) ) ) ) ) + squeeze( sum( sum( sum( comp_cr(ROI_fsamps,:,ROI_tsamps) ) ) ) );
    fa_sum2 = squeeze( sum( sum( sum( fa(ROI_fsamps,:,ROI_tsamps) ) ) ) ) + squeeze( sum( sum( sum( comp_fa(ROI_fsamps,:,ROI_tsamps) ) ) ) );

    [metrics.ROI_perf(t)]=calc_classifier_performance(hit_sum2,miss_sum2,cr_sum2,fa_sum2);
    metrics.struct_xx_sum = '[freq x samps x threshold]';
    metrics.ROI_hit(:,:,t) = hit_sum2;
    metrics.ROI_miss(:,:,t) = miss_sum2;
    metrics.ROI_cr(:,:,t) = cr_sum2;
    metrics.ROI_fa(:,:,t) = fa_sum2;


    %% plotting combined data
    if plot_flag==1
        %         comb = (hit*2) + (miss*-1) + (cr*0) + (fa*-2);

        plv_axis = fc_scale;
        mcc_axis = [-1 1]; 
        min_max = [-.6 .6];

        %% plot true FC TFR averaged across
        axes(ax(3)); cla; hold on; axis on; box on;
        xdata = FC_true; xdata(abs(xdata)<thresh_vals(t)) = nan;
        xdata = double(squeeze(nanmean(xdata,2)));
        surf(lat,freqs,xdata); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis(plv_axis);
        plot3([0 0],[freqs([1 end])],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k--','LineWidth',2);
        plot3([ROI_time_int(1) ROI_time_int(1); ROI_time_int(2) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % horizontals of ROI rectangle
        plot3([ROI_time_int(1) ROI_time_int(2); ROI_time_int(1) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int]',[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % verticals of ROI rectangle
        text(ROI_time_int(1),ROI_freq_int(2)+(.05*range(ax(3).YLim)),'ROI','FontWeight','bold');
        axis([lat([1 end]) freqs([1 end])']);
        ax(2).XTick = ax(1).XTick;
        title ('True FC');
        ylabel('Frequency (Hz)');
        
        %% plot True FC waves
        clear p1 p2 p3 p4 p5 p6 p7
        axes(ax(1)); cla; hold on; axis on; box on;
        for v=1:num_seeds
                   p1(v,:) = plot(lat,squeeze(FC_true(ROI_fsamps,v,:)),'g','linewidth',2);
        end
        axis([lat([1 end]) min_max])
        plot(lat([1 end]),[thresh_vals(t) thresh_vals(t)],'k-','LineWidth',2); plot(lat([1 end]),[-thresh_vals(t) -thresh_vals(t)],'k-','LineWidth',2);
        plot([0 0],[-ax(1).YLim(end)*2 ax(1).YLim(end)*2],'k--','LineWidth',2);
%         legend([p3(1,1);p4(1,1);p5(1);p6(1,1);p1(1,1);p2(1,1);p7(1,1);],{'Hit' 'CR' 'Miss' 'FA' 'true' 'seed' 'other' },'Location','southwest','NumColumns',2,'Orientation','horizontal','location','best');
        legend([p1(1,1)],{'True'},'Location','southwest','NumColumns',1,'Orientation','horizontal','location','southeast');
        xlabel('Time (sec)'); ylabel('FC');
        title(sprintf('ROI True FC waves (thresh = %.2f)',thresh_vals(t)));
        ax(1).Position(3:4) = ax(2).Position(3:4); 
        %% plot Observed FC waves
        axes(ax(2)); cla; hold on; axis on; box on;
        clear p1 p2 p3 p4 p5 p6 p7
        %         for v=1:size(FC_comp,2); p7 = plot(lat,squeeze(FC_comp(:,v,:)),'color',[1 0 0]); end
        xdata = reshape(FC_comp(ROI_fsamps,:,:),[length(ROI_fsamps)*size(comp_fa,2) size(FC_comp,3)]);
        p7 = plot(lat,xdata,'color',[1 0 0]); 
        ydata = reshape(comp_fa(ROI_fsamps,:,ROI_tsamps),[length(ROI_fsamps)*size(comp_fa,2) length(ROI_tsamps)]);
        p6 = plot(lat(ROI_tsamps),ydata*min_max(2)*.8,'rs','MarkerFaceColor','r');
        for v=1:num_seeds
            xhit = squeeze(hit(ROI_fsamps,v,ROI_tsamps)); xhit(xhit==0) = nan; 
            xcr = squeeze(cr(ROI_fsamps,v,ROI_tsamps)); xcr(xcr==0) = nan; 
            xmiss = squeeze(miss(ROI_fsamps,v,ROI_tsamps)); xmiss(xmiss==0) = nan; 
            xfa = squeeze(fa(ROI_fsamps,v,ROI_tsamps)); xfa(xfa==0) = nan; 
            p3 = plot(lat(ROI_tsamps),xhit*min_max(2)*.95,'bs','MarkerFaceColor','b');
            p4 = plot(lat(ROI_tsamps),xcr*min_max(2)*.9,'ks','MarkerFaceColor','k');
            p5 = plot(lat(ROI_tsamps),xmiss*min_max(2)*.85,'s','color','g','MarkerFaceColor','g');
            p6 = plot(lat(ROI_tsamps),xfa*min_max(2)*.8,'rs','MarkerFaceColor','r');

            p2(v,:) = plot(lat,squeeze(FC_seed(ROI_fsamps,v,:)),'b','linewidth',2);
        end
        axis([lat([1 end]) min_max])
        plot(lat([1 end]),[thresh_vals(t) thresh_vals(t)],'k-','LineWidth',2); plot(lat([1 end]),[-thresh_vals(t) -thresh_vals(t)],'k-','LineWidth',2);
        plot([0 0],[-ax(1).YLim(end)*2 ax(1).YLim(end)*2],'k--','LineWidth',2);
%         legend([p3(1,1);p4(1,1);p5(1);p6(1,1);p1(1,1);p2(1,1);p7(1,1);],{'Hit' 'CR' 'Miss' 'FA' 'true' 'seed' 'other' },'Location','southwest','NumColumns',2,'Orientation','horizontal','location','best');
        legend([p3(1,1);p4(1,1);p5(1);p6(1,1); p2(1,1); p7(1,1)],{'Hit' 'CR' 'Miss' 'FA' 'Seed' 'Other' },'Location','southwest','NumColumns',2,'Orientation','horizontal','location','southeast');
        xlabel('Time (sec)'); ylabel('FC');
        title(sprintf('ROI Observed FC waves (thresh = %.2f)',thresh_vals(t)));
        ax(1).Position(3:4) = ax(2).Position(3:4); 
        %% plot hit FC TFR averaged across
        axes(ax(4)); cla; hold on; axis on; box on;
        xdata = FC_seed; 
        xdata(hit==0)=nan;
%         xdata(abs(xdata)<thresh_vals(t)) = nan;
        xdata = double(squeeze(nanmean(xdata,2)));
        surf(lat,freqs,xdata); view(0,90); shading flat; colormap(jet); axis tight; caxis(plv_axis);
        cx = colorbar; cx.Label.String = 'FC value';  
        plot3([0 0],[freqs([1 end])],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k--','LineWidth',2);
       plot3([ROI_time_int(1) ROI_time_int(1); ROI_time_int(2) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % horizontals of ROI rectangle
        plot3([ROI_time_int(1) ROI_time_int(2); ROI_time_int(1) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int]',[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % verticals of ROI rectangle
        text(ROI_time_int(1),ROI_freq_int(2)+(.05*range(ax(3).YLim)),'ROI','FontWeight','bold');
        axis([lat([1 end]) freqs([1 end])']);
        ax(3).XTick = ax(1).XTick;
        title ('True-Positive FC');
        %% plot false-alarm FC TFR averaged across 
        axes(ax(5)); cla; hold on; axis on; box on;
        xdata = FC_comp; 
        xdata(comp_fa==0) = nan; %xdata(abs(xdata)<thresh_vals(t)) = nan;
        xdata = double(squeeze(nanmean(xdata,2)));
        surf(lat,freqs,xdata); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis(plv_axis);
        plot3([0 0],[freqs([1 end])],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k--','LineWidth',2);
       plot3([ROI_time_int(1) ROI_time_int(1); ROI_time_int(2) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % horizontals of ROI rectangle
        plot3([ROI_time_int(1) ROI_time_int(2); ROI_time_int(1) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int]',[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % verticals of ROI rectangle
        text(ROI_time_int(1),ROI_freq_int(2)+(.05*range(ax(3).YLim)),'ROI','FontWeight','bold');
        axis([lat([1 end]) freqs([1 end])']);
        ax(3).XTick = ax(1).XTick;
        title ('False-Positive FC');
        %% plot MCC TFR      
        xdata = double(metrics.perf(t).MCC); %xdata(isnan(xdata)) = 0;
        axes(ax(8)); cla; hold on; axis on; box on;
        surf(lat(tsamps),freqs(fsamps),xdata); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis(mcc_axis);
        plot3([0 0],[freqs([1 end])],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k--','LineWidth',2);
       plot3([ROI_time_int(1) ROI_time_int(1); ROI_time_int(2) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % horizontals of ROI rectangle
        plot3([ROI_time_int(1) ROI_time_int(2); ROI_time_int(1) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int]',[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % verticals of ROI rectangle
        text(ROI_time_int(1),ROI_freq_int(2)+(.05*range(ax(3).YLim)),'ROI','FontWeight','bold');
        axis([lat([1 end]) freqs([1 end])']);
        ax(5).XTick = ax(1).XTick;
        title ('Mathew''s Correlation Coefficient (MCC)');
        ylabel('Frequency (Hz)');
        %% plot TPR TFR
        axes(ax(9)); cla; hold on; axis on; box on;
        %              surf(lat,freqs,squeeze(comb(:,v,:))); view(0,90); shading flat; colormap(jet); axis tight; caxis([-2 2]);
        xdata = squeeze(metrics.perf(t).TPR); xdata(squeeze(sum(hit(fsamps,:,tsamps),2)==0))=nan;
        surf(lat(tsamps),freqs(fsamps),xdata); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis([0 1]);
        plot3([0 0],[freqs([1 end])],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k--','LineWidth',2);
       plot3([ROI_time_int(1) ROI_time_int(1); ROI_time_int(2) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % horizontals of ROI rectangle
        plot3([ROI_time_int(1) ROI_time_int(2); ROI_time_int(1) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int]',[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % verticals of ROI rectangle
        text(ROI_time_int(1),ROI_freq_int(2)+(.05*range(ax(3).YLim)),'ROI','FontWeight','bold');
        axis([lat([1 end]) freqs([1 end])']);
         ax(3).XTick = ax(1).XTick;
       title ('True Positive Rate');
        %% plot FPR TFR
        axes(ax(10)); cla; hold on; axis on; box on;
        %              surf(lat(tsamps),freqs(fsamps),squeeze(comb(:,v,:))); view(0,90); shading flat; colormap(jet); axis tight; caxis([-2 2]);
        xdata = squeeze(metrics.perf(t).FPR); 
        x = logical(squeeze(nansum(comp_fa(fsamps,:,tsamps),2)));
        xdata(x==0) = nan; % masking
        surf(lat(tsamps),freqs(fsamps),xdata); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis([0 1]);
        plot3([0 0],[freqs([1 end])],[ax(4).CLim(end)*2 ax(4).CLim(end)*2],'k--','LineWidth',2);
       plot3([ROI_time_int(1) ROI_time_int(1); ROI_time_int(2) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % horizontals of ROI rectangle
        plot3([ROI_time_int(1) ROI_time_int(2); ROI_time_int(1) ROI_time_int(2)],[ROI_freq_int; ROI_freq_int]',[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k-','LineWidth',2);    % verticals of ROI rectangle
        text(ROI_time_int(1),ROI_freq_int(2)+(.05*range(ax(3).YLim)),'ROI','FontWeight','bold');
        axis([lat([1 end]) freqs([1 end])']);
        ax(4).XTick = ax(1).XTick;
        title ('False Positive Rate');


        %% plot FC graph averaged across time freq that had present MCC values
        axes(ax(6)); cla; hold on; axis off; box off;
        opt.vol_nums = 1:length(vol); vw_angle = [-90 90];
        fc_scale = plv_axis;
             %% true FC map
        fc_contrasts = plv_true_idx;
%         mcc_idx = double(metrics.perf(t).MCC); mcc_idx(abs(mcc_idx)>0) = 1; mcc_idx(isnan(mcc_idx))=0;
%         mcc_idx = permute(repmat(mcc_idx,[1 1 size(fc_contrasts,1)]),[1 3 2]);
%         fc_data = FC_true; %fc_data(mcc_idx==0)=nan; % masking based on existing MCC data
        fc_data = FC_true; fc_data(abs(fc_data)<thresh_vals(t)) = nan;
        fc_data = squeeze( nanmean( nanmean( fc_data(ROI_fsamps,:,ROI_tsamps),1) ,3));
        bl_plot_mesh(vol,opt); view(vw_angle);
        [p1, s1] = bl_plot_FC_graph(fc_data, vol, plv_true_locs, fc_contrasts, 0, fc_map_scale);
        p1=handle(p1); for p=1:length(p1); p1(p).Color = 'g'; p1(p).Color(4) = 1;end
        s1.MarkerEdgeColor = [0 0 0]; s1.MarkerFaceColor = [0 0 0];
 
       %% plot False-Alarm FC graph averaged across time freq that had present MCC values
        axes(ax(7)); cla; hold on; axis off; box off;
        opt.vol_nums = 1:length(vol); vw_angle = [-90 90];
        fc_scale = plv_axis;
        bl_plot_mesh(vol,opt); view(vw_angle);
        %% false-alarm FC map
        fc_contrasts = plv_comp_indices;
%         fc_data = FC_comp; 
%         fc_data(comp_fa==0) = nan; %xdata(abs(xdata)<thresh_vals(t)) = nan;
%         fc_data = squeeze( nanmean( nanmean( fc_data(ROI_fsamps,:,ROI_tsamps),1) ,3));
        fc_data = FC_comp; fc_data(abs(fc_data)<thresh_vals(t)) = nan;
        fc_data = squeeze( nanmean( nanmean( fc_data(ROI_fsamps,:,ROI_tsamps),1) ,3));
        [p1, s1] = bl_plot_FC_graph(fc_data, vol, plv_comp_locs, fc_contrasts, 0, fc_map_scale);
        p1=handle(p1); for p=1:length(p1); try; p1(p).Color = 'r'; p1(p).LineStyle = '-'; p1(p).Color(4) = .5; end; end
         title(sprintf('True ROI FC map'));
         axis tight; 
      %% hit FC map
        fc_contrasts = plv_seed_idx;
%         mcc_idx = double(metrics.perf(t).MCC); mcc_idx(abs(mcc_idx)>0) = 1; mcc_idx(isnan(mcc_idx))=0;
%         mcc_idx = permute(repmat(mcc_idx,[1 1 size(fc_contrasts,1)]),[1 3 2]);
%         fc_data = FC_seed; fc_data(mcc_idx==0)=nan; % masking based on existing MCC data
%         fc_data = squeeze( nanmean( nanmean( fc_data(ROI_fsamps,:,ROI_tsamps),1) ,3));
        fc_data = FC_seed; fc_data(abs(fc_data)<thresh_vals(t)) = nan;
        fc_data = squeeze( nanmean( nanmean( fc_data(ROI_fsamps,:,ROI_tsamps),1) ,3));
        [p1, s1] = bl_plot_FC_graph(fc_data, vol, plv_seed_locs, fc_contrasts, 0, fc_map_scale);
        p1=handle(p1); for p=1:length(p1); p1(p).Color = 'b'; p1(p).LineStyle = '-'; p1(p).Color(4) = 1;end
        title(sprintf('Observed ROI FC map\n(MCC = %.2f)',metrics.ROI_perf(t).MCC));
        axis tight; 

    end

end
% for t=1:length(thresh_vals); figure(103); clf; surf(metrics(a).perf(t).MCC); view(0,90); shading flat; colormap(jet); colorbar; axis tight; axis([1 size(metrics(a).perf(t).MCC,2) 1 size(metrics(a).perf(t).MCC,1)]); caxis([-1 1]); title(thresh_vals(t)); pause; end

%% Mean-Squared and Absolute Errors for FC waves within ROI
metrics.ROI_FC_wave_error_mse = squeeze( nanmean( nanmean( (FC_true(fsamps,:,tsamps)-FC_seed(fsamps,:,tsamps)).^2 ,1) ,3) );  % mean squared error
metrics.ROI_FC_wave_error_abs = squeeze( nanmean( nanmean( abs(FC_true(fsamps,:,tsamps)-FC_seed(fsamps,:,tsamps)) ,1) ,3) );  % mean absolute error
metrics.thresh_vals = thresh_vals; 




