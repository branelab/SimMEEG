function [metrics] = calc_MCC_TFR_PLV_surg_noise(FC_true, FC_true_noise_surg, FC_true_noise_surg_std, FC_seed, FC_comp, FC_comp_noise_surg, FC_comp_noise_surg_std, std_step_size, std_stop, plv_comp_indices, plv_comp_idx, plv_true_locs, plv_seed_locs, plv_comp_locs, vol, freq_int, time_int, ROI_freq_int, ROI_time_int, freqs, lat, plot_flag, fc_scale, fc_map_scale, ROI_only_flag)
% This calculates MCC for all the sample points and all frequencies
% within the FC data scanning through threshold based on standard deviations from values "thresh_vals"
%
% INPUT:
%   FC_true = true FC [freqs x comparisons x samps]
%   FC_seed = observed FC for the seed FCs (Hits) with same order as FC true [freqs x comparisons x samps]
%   FC_comp = observed FC for other comparisons (Correct Rejections) [freqs x comparisons x samps]
%   FC_comp_noise_surg = mean of surrogate (noise) for other comparisons for thresholding the data to calculate metrics [freqs x comparisons x samps]
%   FC_comp_noise_surg_std = standard deviation of surrogate (noise) for other comparisons for thresholding the data to calculate metrics [freqs x comparisons x samps]
%   std_step_size = step size of standard deviations starting from 0 until only correct rejections are found:  thresholds = FC_comp_noise_surg +/- (std_val * FC_comp_noise_surg_std).
%   plv_comp_idx = indices of plv_comp_locs that correspond to the plv_comp_idx;
%   plv_comp_indices = indices of FC comparisons within the plv_comp_locs
%   thresh_vals = (t) FC thresholds for calculating MCC. This will be for negative and positive values for searching for hits, misses, CRs, FAs.
%   freq_int = freq of interest interval for summing across to get a global metrics across this interval.
%   time_int = time of interest interval for summing across to get a global metrics across this interval.
%   lat = latency of samps
%   plot_flag = (1) plots some results (0) no plots
%   ROI_only_flag = (0) calculate stats across all freq_int and time_int and calculate for ROI (1) calculate only for ROI
%
% OUTPUT:
%    metrics.accruacy(t).hit_idx        = accuracy for hits, miss, cr_seed, cr_comp, fa_seed, and fa_comp for each [freqs x comparisons x samps].
%
%    metrics.perf(std_iter).MCC                = Mathew's Correlation Coefficient
%                   .true_positives     = True positives (Hits) for each thresh_vals
%                   .false_positives    = False positives (False Alarms) for each thresh_vals
%                   .true_negatives     = True positives (Correct Rejections) for each thresh_vals
%                   .false_negatives    = True positives (Misses) for each thresh_vals
%
%
%   metrics.hit                         = hits summed across seeds [freqs x samps x std_iter threshold]
%          .miss                        = miss summed across seeds [freqs x samps x std_iter threshold]
%          .cr                          = cr summed across seeds [freqs x samps x std_iter threshold]
%          .fa                          = fa summed across seeds [freqs x samps x std_iter threshold]


% thresh_vals = 0:.1:1;
% FC_data = inv_soln(a).TFR_results.plv_based;
% FC_true = cfg.source.TFR_results.plv_based;
num_seeds = size(FC_true,2);
plv_true_idx = nchoosek(1:size(FC_true,2),2);
plv_seed_idx = nchoosek(1:num_seeds,2);
% find and maintain order of seeds within FC_comp
for v=1:size(plv_seed_idx,1)
    plv_comp_seed_idx(v) = find(ismember(plv_comp_idx,plv_seed_idx(v,:),'rows'));
end

%% freq/time interval of interest for calculating metrics
ss = find(freqs>=freq_int(1)); fx(1) = ss(1); ss = find(freqs<=freq_int(2)); fx(2) = ss(end); fsamps = fx(1):fx(2); 
ss = find(lat>=time_int(1)); sx(1) = ss(1); ss = find(lat<=time_int(2)); sx(2) = ss(end); tsamps = sx(1):sx(2); 
%% ROI freq/time interval of interest for summing across to calculate metrics
ss = find(freqs(fsamps)>=ROI_freq_int(1)); rfx(1) = ss(1); ss = find(freqs(fsamps)<=ROI_freq_int(2)); rfx(2) = ss(end); ROI_fsamps = rfx(1):rfx(2); 
ss = find(lat(tsamps)>=ROI_time_int(1)); rsx(1) = ss(1); ss = find(lat(tsamps)<=ROI_time_int(2)); rsx(2) = ss(end); ROI_tsamps = rsx(1):rsx(2); 

metrics.lat = lat(tsamps);
wbar = waitbar(0); 

cr_seed_idx = zeros(size(FC_true(fsamps,:,tsamps))); % correct rejections found = 1;
std_val = 0; std_iter = 0;
while nansum(nansum(nansum(cr_seed_idx)))< numel(cr_seed_idx) || std_val>=std_stop
    std_iter = std_iter + 1;

    metrics.accuracy(std_iter).std_val = std_val;
    metrics.std_val(std_iter) = std_val;

    FC_comp_noise_surg_thresh = FC_comp_noise_surg + (std_val*FC_comp_noise_surg_std);
    FC_true_noise_surg_thresh = FC_true_noise_surg + (std_val*FC_true_noise_surg_std);

    %% Thresholding FC data

    %% hits
    metrics.accuracy(std_iter).hit_idx =   double (squeeze( (FC_true(fsamps,:,tsamps) > FC_true_noise_surg_thresh(fsamps,:,tsamps)) & (FC_seed(fsamps,:,tsamps) > FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx,tsamps)) ) | ...
        squeeze( (FC_true(fsamps,:,tsamps) < -FC_true_noise_surg_thresh(fsamps,:,tsamps)) & (FC_seed(fsamps,:,tsamps) < -FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx,tsamps)) ) );
    %% Miss
    metrics.accuracy(std_iter).miss_idx =   double (squeeze( (FC_true(fsamps,:,tsamps) > FC_true_noise_surg_thresh(fsamps,:,tsamps)) & (FC_seed(fsamps,:,tsamps) < FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx,tsamps)) & (FC_seed(fsamps,:,tsamps) > -FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx,tsamps)) ) | ...
        squeeze( (FC_true(fsamps,:,tsamps) < -FC_true_noise_surg_thresh(fsamps,:,tsamps)) & (FC_seed(fsamps,:,tsamps) < FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx,tsamps)) & (FC_seed(fsamps,:,tsamps) > -FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx,tsamps)) ) );
    %% CR seed
    metrics.accuracy(std_iter).cr_seed_idx =   double (squeeze( (FC_true(fsamps,:,tsamps) < FC_true_noise_surg_thresh(fsamps,:,tsamps)) & (FC_true(fsamps,:,tsamps) > -FC_true_noise_surg_thresh(fsamps,:,tsamps))  & ...
        (FC_seed(fsamps,:,tsamps) < FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx,tsamps)) & (FC_seed(fsamps,:,tsamps) > -FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx,tsamps)) ) ) ;
    cr_seed_idx = metrics.accuracy(std_iter).cr_seed_idx; % needed for stoping criteria that all time-freq samples have been correctly rejected
    %% FA seed
    metrics.accuracy(std_iter).fa_seed_idx =   double (squeeze( (FC_seed(fsamps,:,tsamps) > FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx,tsamps)) & (FC_true(fsamps,:,tsamps) < FC_true_noise_surg_thresh(fsamps,:,tsamps)) & (FC_true(fsamps,:,tsamps) > -FC_true_noise_surg_thresh(fsamps,:,tsamps)) ) | ...
        squeeze( (FC_seed(fsamps,:,tsamps) < -FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx,tsamps)) & (FC_true(fsamps,:,tsamps) < FC_true_noise_surg_thresh(fsamps,:,tsamps)) & (FC_true(fsamps,:,tsamps) > -FC_true_noise_surg_thresh(fsamps,:,tsamps)) ) );

    %% FA comp
    metrics.accuracy(std_iter).fa_comp_idx =   double (squeeze( (FC_comp(fsamps,:,tsamps) > FC_comp_noise_surg_thresh(fsamps,:,tsamps)) ) | ...
        squeeze( (FC_comp(fsamps,:,tsamps) < -FC_comp_noise_surg_thresh(fsamps,:,tsamps)) ) );
    %% CR comp
    metrics.accuracy(std_iter).cr_comp_idx = (metrics.accuracy(std_iter).fa_comp_idx-1)~=0; % not a false alarm then must be a correct rejection

    %% metrics for each time/freq sample
    if length(unique(fx))>1
    hit_sum = squeeze(nansum(metrics.accuracy(std_iter).hit_idx,2));
    miss_sum = squeeze(nansum(metrics.accuracy(std_iter).miss_idx,2));
    cr_sum = squeeze(nansum(metrics.accuracy(std_iter).cr_seed_idx,2) + nansum(metrics.accuracy(std_iter).cr_comp_idx,2));
    fa_sum = squeeze(nansum(metrics.accuracy(std_iter).fa_seed_idx,2) + nansum(metrics.accuracy(std_iter).fa_comp_idx,2));
    else
    hit_sum = squeeze(nansum(metrics.accuracy(std_iter).hit_idx,1));
    miss_sum = squeeze(nansum(metrics.accuracy(std_iter).miss_idx,1));
    cr_sum = squeeze(nansum(metrics.accuracy(std_iter).cr_seed_idx,1) + nansum(metrics.accuracy(std_iter).cr_comp_idx,1));
    fa_sum = squeeze(nansum(metrics.accuracy(std_iter).fa_seed_idx,1) + nansum(metrics.accuracy(std_iter).fa_comp_idx,1));
    end
    
    [metrics.perf(std_iter)]=calc_classifier_performance(hit_sum,miss_sum,cr_sum,fa_sum);
    metrics.struct_xx_sum = '[freq x samps x std_iter threshold]';
    metrics.hit(:,:,std_iter) = hit_sum;
    metrics.miss(:,:,std_iter) = miss_sum;
    metrics.cr(:,:,std_iter) = cr_sum;
    metrics.fa(:,:,std_iter) = fa_sum;


    %% ROI time-freq metrics summed across time_int and freq_int 
     if length(unique(fx))>1
   hit_sum2 = squeeze(nansum(nansum(nansum(metrics.accuracy(std_iter).hit_idx(ROI_fsamps,:,ROI_tsamps),2))));
    miss_sum2 = squeeze(nansum(nansum(nansum(metrics.accuracy(std_iter).miss_idx(ROI_fsamps,:,ROI_tsamps),2))));
    cr_sum2 = squeeze(nansum(nansum(nansum(metrics.accuracy(std_iter).cr_seed_idx(ROI_fsamps,:,ROI_tsamps),2))) + nansum(nansum(nansum(metrics.accuracy(std_iter).cr_comp_idx(ROI_fsamps,:,ROI_tsamps),2))));
    fa_sum2 = squeeze(nansum(nansum(nansum(metrics.accuracy(std_iter).fa_seed_idx(ROI_fsamps,:,ROI_tsamps),2))) + nansum(nansum(nansum(metrics.accuracy(std_iter).fa_comp_idx(ROI_fsamps,:,ROI_tsamps),2))));
     else
   hit_sum2 = squeeze(nansum(nansum(metrics.accuracy(std_iter).hit_idx(:,ROI_tsamps),1)));
    miss_sum2 = squeeze(nansum(nansum(metrics.accuracy(std_iter).miss_idx(:,ROI_tsamps),1)));
    cr_sum2 = squeeze(nansum(nansum(metrics.accuracy(std_iter).cr_seed_idx(:,ROI_tsamps),1)) + nansum(nansum(metrics.accuracy(std_iter).cr_comp_idx(:,ROI_tsamps),1)));
    fa_sum2 = squeeze(nansum(nansum(metrics.accuracy(std_iter).fa_seed_idx(:,ROI_tsamps),1)) + nansum(nansum(metrics.accuracy(std_iter).fa_comp_idx(:,ROI_tsamps),1)));
     end
     
    [metrics.ROI_perf(std_iter)]=calc_classifier_performance(hit_sum2,miss_sum2,cr_sum2,fa_sum2);
    metrics.ROI_struct_xx_sum = '[std_iter threshold]';
    metrics.ROI_hit(std_iter) = hit_sum2;
    metrics.ROI_miss(std_iter) = miss_sum2;
    metrics.ROI_cr(std_iter) = cr_sum2;
    metrics.ROI_fa(std_iter) = fa_sum2;
    %% Plotting
    


    %% plotting combined data
    if plot_flag==1
    figure(101); clf; set(gcf,'color','w');
    ax = subplot_axes(2,4,.035,.075,0,0,0);

        plv_axis = fc_scale;
        mcc_axis = [-1 1];
        min_max = [-.6 .6];

        %% plot true FC TFR averaged across
        axes(ax(2)); cla; hold on; axis on; box on;

        xdata = FC_true(fsamps,:,tsamps); xdata(metrics.accuracy(std_iter).hit_idx==0) = nan;
        xdata = double(squeeze(nanmean(xdata,2)));
        surf(lat(tsamps),freqs(fsamps),xdata); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis(plv_axis);
        plot3([0 0],[freqs(fsamps([1 end]))],[ax(2).CLim(end)*2 ax(2).CLim(end)*2],'k--','LineWidth',2);
        axis([lat(tsamps([1 end])) freqs(fsamps([1 end]))']);
        ax(2).XTick = ax(1).XTick;
        title ('True FC');
        ylabel('Frequency (Hz)');

        %% plot FC waves
        clear p1 p2 p3 p4 
        axes(ax(1)); cla; hold on; axis on; box on;
        for v=1:size(FC_comp,2)
            p1(v,:) = plot(lat(tsamps),[squeeze(FC_comp_noise_surg_thresh(fsamps,v,tsamps)); -squeeze(FC_comp_noise_surg_thresh(fsamps,v,tsamps))],'color',[0 0 0]);
            p2(v,:) = plot(lat(tsamps),squeeze(FC_comp(fsamps,v,tsamps)),'r');
        end
        for v=1:num_seeds
%             p1(v,:) = plot(lat(tsamps),[squeeze(FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx(v),tsamps)); -squeeze(FC_comp_noise_surg_thresh(fsamps,plv_comp_seed_idx(v),tsamps))],'color',[0 0 0]);
%             p2(v,:) = plot(lat(tsamps),squeeze(FC_comp(fsamps,plv_comp_seed_idx(v),tsamps)),'color',[1 0 0]);
            p3(v,:) = plot(lat(tsamps),squeeze(FC_true(fsamps,v,tsamps)),'g');
            p4(v,:) = plot(lat(tsamps),squeeze(FC_seed(fsamps,v,tsamps)),'b');
        end
%       p3 = plot(lat(tsamps), squeeze(nansum(nansum(metrics.accuracy(std_iter).hit_idx))/(size(metrics.accuracy(std_iter).hit_idx,1)*size(metrics.accuracy(std_iter).hit_idx,2)))+(.8*fc_scale(1)),'go'); % False positive rate
%     p4 = plot(lat(tsamps), squeeze(nansum(nansum(metrics.accuracy(std_iter).miss_idx))/(size(metrics.accuracy(std_iter).miss_idx,1)*size(metrics.accuracy(std_iter).miss_idx,2)))+(.84*fc_scale(1)),'bo'); % False positive rate
%     p5 = plot(lat(tsamps), squeeze(nansum(nansum(metrics.accuracy(std_iter).cr_seed_idx))/(size(metrics.accuracy(std_iter).cr_seed_idx,1)*size(metrics.accuracy(std_iter).cr_seed_idx,2)))+(.82*fc_scale(1)),'gs'); % False positive rate
%     p6 = plot(lat(tsamps), squeeze(nansum(nansum(metrics.accuracy(std_iter).fa_comp_idx))/(size(metrics.accuracy(std_iter).fa_comp_idx,1)*size(metrics.accuracy(std_iter).fa_comp_idx,2)))+(.82*fc_scale(1)),'r*'); % False positive rate
%     p7 = plot(lat(tsamps), squeeze(nansum(nansum(metrics.accuracy(std_iter).fa_seed_idx))/(size(metrics.accuracy(std_iter).fa_seed_idx,1)*size(metrics.accuracy(std_iter).fa_seed_idx,2)))+(.82*fc_scale(1)),'rs'); % False positive rate

        axis([lat(tsamps([1 end])) min_max])
%         plot(lat([1 end]),[thresh_vals(t) thresh_vals(t)],'k--','LineWidth',2); plot(lat([1 end]),[-thresh_vals(t) -thresh_vals(t)],'k--','LineWidth',2);
        plot([0 0],[-ax(1).YLim(end)*2 ax(1).YLim(end)*2],'k--','LineWidth',2);
%         legend([p3(1,1);p4(1,1);p5(1);p6(1,1);p1(1,1);p2(1,1);p7(1,1);],{'Hit' 'CR' 'Miss' 'FA' 'true' 'seed' 'other' },'Location','southwest','NumColumns',2,'Orientation','horizontal');
          legend([p3(1,1);p4(1,1);p2(1);p1(1,1)],{'true' 'seed' 'other' 'thresholds' },'Location','southwest','NumColumns',2,'Orientation','horizontal');
        xlabel('Time (sec)'); ylabel('FC');
        title(sprintf('FC waves (Std Dev = %.2f)',std_val));
        ax(1).Position(3:4) = ax(2).Position(3:4);


        %% plot hit FC TFR averaged across
        axes(ax(3)); cla; hold on; axis on; box on;
        xdata = FC_seed(fsamps,:,tsamps);
        xdata(metrics.accuracy(std_iter).hit_idx==0)=nan;
        %         xdata(abs(xdata)<thresh_vals(t)) = nan;
        xdata = double(squeeze(nanmean(xdata,2)));
        surf(lat(tsamps),freqs(fsamps),xdata); view(0,90); shading flat; colormap(jet); axis tight; caxis(plv_axis);
        cx = colorbar; cx.Label.String = 'FC value';
        plot3([0 0],[freqs([1 end])],[ax(2).CLim(end)*2 ax(2).CLim(end)*2],'k--','LineWidth',2);
        axis([lat(tsamps([1 end])) freqs(fsamps([1 end]))']);
        ax(3).XTick = ax(1).XTick;
        title ('True-Positive FC');

        %% plot false-alarm FC TFR averaged across
        axes(ax(4)); cla; hold on; axis on; box on;
        xdata = FC_comp(fsamps,:,tsamps);
        xdata(metrics.accuracy(std_iter).fa_comp_idx==0) = nan; %xdata(abs(xdata)<thresh_vals(t)) = nan;
        xdata = double(squeeze(nanmean(xdata,2)));
        surf(lat(tsamps),freqs(fsamps),xdata); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis(plv_axis);
        plot3([0 0],[freqs([1 end])],[ax(2).CLim(end)*2 ax(2).CLim(end)*2],'k--','LineWidth',2);
        axis([lat(tsamps([1 end])) freqs(fsamps([1 end]))']);
        ax(3).XTick = ax(1).XTick;
        title ('False-Positive FC');

        %% plot MCC TFR
        xdata = double(metrics.perf(std_iter).MCC); %xdata(isnan(xdata)) = 0;
        axes(ax(6)); cla; hold on; axis on; box on;
        surf(lat(tsamps),freqs(fsamps),xdata); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis(mcc_axis);
        plot3([0 0],[freqs([1 end])],[ax(2).CLim(end)*2 ax(2).CLim(end)*2],'k--','LineWidth',2);
                axis([lat(tsamps([1 end])) freqs(fsamps([1 end]))']);
        ax(5).XTick = ax(1).XTick;
        title ('Mathew''s Correlation Coefficient (MCC)');
        ylabel('Frequency (Hz)');


        %% plot TPR TFR
        axes(ax(7)); cla; hold on; axis on; box on;
        %              surf(lat(tsamps),freqs(fsamps),squeeze(comb(:,v,:))); view(0,90); shading flat; colormap(jet); axis tight; caxis([-2 2]);
        xdata = squeeze(metrics.perf(std_iter).TPR); xdata(squeeze(nansum(metrics.accuracy(std_iter).hit_idx,2))==0)=nan;
        surf(lat(tsamps),freqs(fsamps),xdata); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis([0 1]);
        plot3([0 0],[freqs([1 end])],[ax(3).CLim(end)*2 ax(3).CLim(end)*2],'k--','LineWidth',2);
                axis([lat(tsamps([1 end])) freqs(fsamps([1 end]))']);
        ax(3).XTick = ax(1).XTick;
        title ('True Positive Rate');

        %% plot FPR TFR
        axes(ax(8)); cla; hold on; axis on; box on;
        %              surf(lat(tsamps),freqs(fsamps),squeeze(comb(:,v,:))); view(0,90); shading flat; colormap(jet); axis tight; caxis([-2 2]);
        xdata = squeeze(metrics.perf(std_iter).FPR); xdata(squeeze(nansum(metrics.accuracy(std_iter).fa_comp_idx,2))==0)=nan;
        surf(lat(tsamps),freqs(fsamps),xdata); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis([0 1]);
        plot3([0 0],[freqs([1 end])],[ax(4).CLim(end)*2 ax(4).CLim(end)*2],'k--','LineWidth',2);
                axis([lat(tsamps([1 end])) freqs(fsamps([1 end]))']);
        ax(4).XTick = ax(1).XTick;
        title ('False Positive Rate');


        %% plot FC graph averaged across time freq that had present MCC values
        axes(ax(5)); cla; hold on; axis off; box off;
        opt.vol_nums = 1:length(vol); vw_angle = [-90 90];
        fc_scale = plv_axis;
        %% false-alarm FC map
        fc_contrasts = plv_comp_indices;
        fc_data = FC_comp(fsamps,:,tsamps);
        fc_data(metrics.accuracy(std_iter).fa_comp_idx==0) = nan; %xdata(abs(xdata)<thresh_vals(t)) = nan;
        fc_data = squeeze( nanmean( nanmean( fc_data,1) ,3));
        [p1, s1] = bl_plot_FC_graph(fc_data, vol, plv_comp_locs, fc_contrasts, 0, fc_map_scale);
        p1=handle(p1); for p=1:length(p1); try; p1(p).Color = 'r'; p1(p).LineStyle = '-'; p1(p).Color(4) = .25; end; end
        title(sprintf('FC map'));
        %% true FC map
        fc_contrasts = plv_true_idx;
        mcc_idx = double(metrics.perf(std_iter).MCC); mcc_idx(abs(mcc_idx)>0) = 1; mcc_idx(isnan(mcc_idx))=0;
        mcc_idx = permute(repmat(mcc_idx,[1 1 size(fc_contrasts,1)]),[1 3 2]);
        fc_data = FC_true(fsamps,:,tsamps); fc_data(mcc_idx==0)=nan; % masking based on existing MCC data
        fc_data = squeeze( nanmean( nanmean( fc_data,1) ,3));
        bl_plot_mesh(vol,opt); view(vw_angle);
        [p1, s1] = bl_plot_FC_graph(fc_data, vol, plv_true_locs, fc_contrasts, 0, fc_map_scale);
        p1=handle(p1); for p=1:length(p1); p1(p).Color = 'g'; p1(p).Color(4) = 1;end
        s1.MarkerEdgeColor = [0 1 0]; s1.MarkerFaceColor = [0 1 0];
        title(sprintf('FC map averaged over existing MCC samples'));
        %% hit FC map
        fc_contrasts = plv_seed_idx;
        mcc_idx = double(metrics.perf(std_iter).MCC); mcc_idx(abs(mcc_idx)>0) = 1; mcc_idx(isnan(mcc_idx))=0;
        mcc_idx = permute(repmat(mcc_idx,[1 1 size(fc_contrasts,1)]),[1 3 2]);
        fc_data = FC_seed(fsamps,:,tsamps); fc_data(mcc_idx==0)=nan; % masking based on existing MCC data
        fc_data = squeeze( nanmean( nanmean( fc_data,1) ,3));
        [p1, s1] = bl_plot_FC_graph(fc_data, vol, plv_seed_locs, fc_contrasts, 0, fc_map_scale);
        p1=handle(p1); for p=1:length(p1); p1(p).Color = 'b'; p1(p).LineStyle = ':'; p1(p).Color(4) = 1;end
        title(sprintf('FC map'));



    end
    %% stepping up std deviation value
    std_val = std_val + std_step_size;
    if std_val > std_stop; break; end
%     fprintf('StdDev = %.2f, Percent CRs = %.3f\n',std_val, 100*nansum(nansum(nansum(cr_seed_idx)))/numel(cr_seed_idx))
    waitbar(nansum(nansum(nansum(cr_seed_idx)))/numel(cr_seed_idx), wbar, sprintf('%.2f %% Seeded CRs found at StdDev = %.2f',100*nansum(nansum(nansum(cr_seed_idx)))/numel(cr_seed_idx), std_val));


end

if ROI_only_flag == 1   % removing results to reduce memory storage
    metrics = rmfield(metrics,'accuracy');
    metrics = rmfield(metrics,'perf');
    metrics = rmfield(metrics,'struct_xx_sum');
    metrics = rmfield(metrics,'hit');
    metrics = rmfield(metrics,'miss');
    metrics = rmfield(metrics,'cr');
    metrics = rmfield(metrics,'fa');
end

close(wbar)

% for t=1:length(thresh_vals); figure(103); clf; surf(metrics(a).perf(std_iter).MCC); view(0,90); shading flat; colormap(jet); colorbar; axis tight; axis([1 size(metrics(a).perf(std_iter).MCC,2) 1 size(metrics(a).perf(std_iter).MCC,1)]); caxis([-1 1]); title(thresh_vals(t)); pause; end

%% Mean-Squared and Absolute Errors for FC waves within ROI
metrics.ROI_FC_wave_error_mse = squeeze( nanmean( nanmean( (FC_true(fsamps,:,tsamps)-FC_seed(fsamps,:,tsamps)).^2 ,1) ,3) );  % mean squared error
metrics.ROI_FC_wave_error_abs = squeeze( nanmean( nanmean( abs(FC_true(fsamps,:,tsamps)-FC_seed(fsamps,:,tsamps)) ,1) ,3) );  % mean absolute error

