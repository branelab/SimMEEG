function [metrics] = calc_MCC_TFR_thresh(FC_true, FC_seed, FC_comp, freq_int, time_int, ROI_freq_int, ROI_time_int, thresh_vals, freqs, lat, plot_flag, ROI_only_flag)
% This calculates MCC for all the sample points and all frequencies
% within the FC data scanning through threshold plv/pli values "thresh_vals"
%
% INPUT:
%   FC_true = true FC [freqs x comparisons x samps]
%   FC_seed = observed FC for the seed FCs (Hits) with same order as FC true [freqs x comparisons x samps]
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
num_seeds = size(FC_true,2);


%% freq/time interval of interest for calculating metrics
ss = find(freqs>=freq_int(1)); fx(1) = ss(1); ss = find(freqs<=freq_int(2)); fx(2) = ss(end); fsamps = fx(1):fx(2); 
ss = find(lat>=time_int(1)); sx(1) = ss(1); ss = find(lat<=time_int(2)); sx(2) = ss(end); tsamps = sx(1):sx(2); 
%% ROI freq/time interval of interest for summing across to calculate metrics
ss = find(freqs>=ROI_freq_int(1)); rfx(1) = ss(1); ss = find(freqs<=ROI_freq_int(2)); rfx(2) = ss(end); ROI_fsamps = rfx(1):rfx(2); 
ss = find(lat>=ROI_time_int(1)); rsx(1) = ss(1); ss = find(lat<=ROI_time_int(2)); rsx(2) = ss(end); ROI_tsamps = rsx(1):rsx(2); 


if plot_flag==1
    clf; set(gcf,'color','w');
    ax = subplot_axes(2,5,.035,.125,0,0,0);
end

for t=1:length(thresh_vals)
    %% Thresholding FC data
    %% Seed Hits --> when true>thresh and seed>thresh or true<-thresh and seed<-thresh
    xhit = (FC_true>thresh_vals(t) & FC_seed>=thresh_vals(t)) | (FC_true<-thresh_vals(t) & FC_seed<=-thresh_vals(t) );
    FC_hit = double(xhit>0);  % Hits for positive & negative FCs
    hit = FC_hit; %hit(FC_hit==0)=nan;
    %    figure(1); clf; surf(squeeze(hit(:,3,:))); view(0,90); shading flat; colormap(jet); axis tight; colorbar;

    %% Seed Misses --> when true>thresh and seed<thresh or true<-thresh and seed>-thresh
    x = (FC_true>thresh_vals(t) & FC_seed<=thresh_vals(t) & FC_seed>=-thresh_vals(t)) | (FC_true<-thresh_vals(t) & FC_seed>=-thresh_vals(t) & FC_seed<=thresh_vals(t));
    FC_miss = double(x>0);  % Hits for positive & negative FCs
    miss = FC_miss; %miss(FC_miss==0)=nan;
    %     figure(2); clf;  surf(squeeze(miss(:,3,:))); view(0,90); shading flat; colormap(jet); axis tight; colorbar;

    %% Seed False Alarms--> when true<thresh and seed>thresh or true>-thresh and seed<-thresh
    xfa = (FC_true<=thresh_vals(t) & FC_true>=-thresh_vals(t) & FC_seed>thresh_vals(t)) | (FC_true<=thresh_vals(t) & FC_true>=-thresh_vals(t) & FC_seed<-thresh_vals(t));
    FC_fa = double(xfa>0);  % Hits for positive & negative FCs
    fa = FC_fa; %fa(FC_fa==0)=nan;
    %     figure(4); clf;  surf(squeeze(fa(:,3,:))); view(0,90); shading flat; colormap(jet); axis tight; colorbar;
    %% Seed Correct Rejections --> when true<thresh and seed<thresh or true>-thresh and seed>-thresh
    x = (FC_true<=thresh_vals(t) & FC_seed<=thresh_vals(t) & FC_true>=-thresh_vals(t) & FC_seed>=-thresh_vals(t) ) | ( xhit==0 & xfa==0 ); % if no FC_true exists and no FAs then they samples are CRs
    FC_cr = double(x>0);  % Hits for positive & negative FCs
    cr = FC_cr; %cr(FC_cr==0)=nan;
    %     figure(3); clf;  surf(squeeze(cr(:,3,:))); view(0,90); shading flat; colormap(jet); axis tight; colorbar;


    %% Comparison Correct Rejections --> when true<thresh and seed<thresh or true>-thresh and seed>-thresh
    x = FC_comp<=thresh_vals(t) & FC_comp>=-thresh_vals(t);
    comp_cr = double(x>0);  % Hits for positive & negative FCs
    %    figure(5); clf; surf(squeeze(comp_cr(:,3,:))); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis([0 1]);

    %% Comparison False Alarms--> when true<thresh and seed>thresh or true>-thresh and seed<-thresh
    x = FC_comp>=thresh_vals(t) | FC_comp<=-thresh_vals(t);
    comp_fa = double(x>0);  % Hits for positive & negative FCs
    %    figure(6); clf; surf(squeeze(comp_fa(:,3,:))); view(0,90); shading flat; colormap(jet); colorbar; axis tight; caxis([0 1]);

    %% summing up data
    hit_sum = squeeze(nansum(hit(fsamps,:,tsamps),2));
    miss_sum = squeeze(nansum(miss(fsamps,:,tsamps),2));
    cr_sum = squeeze(nansum(cr(fsamps,:,tsamps),2) + nansum(comp_cr(fsamps,:,tsamps),2));
    fa_sum = squeeze(nansum(fa(fsamps,:,tsamps),2) + nansum(comp_fa(fsamps,:,tsamps),2));

    % Plot Testing
    %        figure(1); clf; surf(lat(tsamps),freqs(fsamps),squeeze(hit_sum)); view(0,90); shading flat; colormap(jet); axis tight; colorbar;
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

end

%% Mean-Squared and Absolute Errors for FC waves within ROI
metrics.ROI_FC_wave_error_mse = squeeze( nanmean( nanmean( (FC_true(fsamps,:,tsamps)-FC_seed(fsamps,:,tsamps)).^2 ,1) ,3) );  % mean squared error
metrics.ROI_FC_wave_error_abs = squeeze( nanmean( nanmean( abs(FC_true(fsamps,:,tsamps)-FC_seed(fsamps,:,tsamps)) ,1) ,3) );  % mean absolute error
metrics.thresh_vals = thresh_vals; 




