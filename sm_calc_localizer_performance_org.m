function sm_calc_localizer_performance(varargin)
%function [metrics]=Sim_study_metrics_v2(eeg_sim,ssp,true_locs,true_idx,inv_soln,HeadModelMat,inv_type,act_samps,inv_data,ssp_data);
% This program will calculate the SimStudy metrics to evaluate the LCMV & MCMV performances.
%

global h
%   dist_thresh = [X Y Z] distance (mm) threshold to search for Hits relative to true source locations

dist_thresh = str2num(h.edit_inv_dist_search_thresh.String);    % [X Y Z] distance (mm) threshold to search for Hits relative to true source locations
if isempty(dist_thresh)
        dist_thresh = [20 20 20];
        h.edit_inv_dist_search_thresh.String = sprintf('%.f %.f %.f',dist_thresh);
else
    if nansum(dist_thresh)==0
        dist_thresh = [20 20 20];
        h.edit_inv_dist_search_thresh.String = sprintf('%.f %.f %.f',dist_thresh);
    end
    if length(dist_thresh)~=3
        dist_thresh = repmat(dist_thresh(1),3,1);    % [X Y Z] distance (mm) threshold to search for Hits relative to true source locations
            h.edit_inv_dist_search_thresh.String = sprintf('%.f %.f %.f',dist_thresh);
    end
end

act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
metrics.true_idx = h.sim_data.cfg.source.vx_idx;
metrics.dist_thresh = (nansum(dist_thresh.^2)).^.5;  % e.g., 15 mm distance threshold is the radius from the ssp location as defined above and not the --> sqrt((.015^2 + .015^2 + .015^2));    % 15 mm distance threshold in any one direction will yield a distance error  = 26 mm

%% creating an empty array if no inv_soln sources exist
    true_locs = h.sim_data.cfg.source.vx_locs;  % true source locations

if ~isempty(h.inv_soln(h.current_inv_soln).peak_voxels)
    
    
    %% localizer performance  = hits / (FA+Miss)
    %  Hit = source within 15 mm distance from ssp dipole and has the lowest metrics.amp_mse --> meaning it best matches the ssp waveform.
    
    
    % removing duplicates that can exist if true peaks not found
    [x,p_idx] = unique(h.inv_soln(h.current_inv_soln).peak_voxels(:,5),'stable');
    peak_voxels = h.inv_soln(h.current_inv_soln).peak_voxels(p_idx,:);
    peak_idx = h.inv_soln(h.current_inv_soln).peak_voxels(p_idx,5);
    peak_locs = peak_voxels(:,1:3);  % found peak locations
    [v_idx]=find_nearest_voxel(peak_locs,true_locs);    % find nearest true_source for each peak source
    metrics.nearest_true_idx = metrics.true_idx(v_idx);
    %% first criteria find peaks within dist_thresh surrounding true_locs
    metrics.loc_error = nan(size(peak_locs,1),size(true_locs,1));
    for v=1:size(true_locs,1)
        xd = bsxfun(@minus,peak_locs,true_locs(v,:));
        xd = (nansum(xd.^2,2)).^.5;
        metrics.loc_error(:,v) = xd;
    end
    peak_hit = metrics.loc_error<metrics.dist_thresh;
    
    
        %% second criteria --> if more than 1 found peak within dist_thresh search radius then choose the closest
    hit_idx = nan(size(true_locs,1),1);
%     for v=1:size(hit_idx,1) % looping through real dipole sources
%         peak_hit_idx = find(peak_hit(:,v)==1, 1);
%         if ~isempty(peak_hit_idx)
%             if peak_hit_idx<=3  % one of the true peaks that was originally found to be near to it
%                 hit_idx(v) = peak_idx(peak_hit_idx);
%             end
%         end
%     end
    for v=1:size(peak_hit,2)
        if ~isempty( find(peak_hit(:,v)==1))
        x_idx = find(peak_hit(:,v)==1);
        [x,y] = min(metrics.loc_error(x_idx,v));
        hit_idx(v) = peak_idx(x_idx(y));
        end
    end
    miss_idx = metrics.true_idx(isnan(hit_idx));    % Missed true locations
     
    %% Find number of false-alarms when last true source found
     [fa_idx,f_idx] = setxor(peak_idx,hit_idx,'stable');  % False-Alarms
     [~,h_idx] = setxor(peak_idx,fa_idx,'stable'); 

     if ~isempty(h_idx) % if no true sources found then select False Alarms for all shown peak voxels
         fb_idx =  peak_voxels(f_idx,4)>=min(peak_voxels(h_idx,4)) ;   % voxel value of inv_soln
         FA_idx = peak_voxels(f_idx(fb_idx),5);
     else
         FA_idx = peak_voxels(:,5);
     end
else
    peak_voxels =[]; peak_idx =[];
    peak_locs = [];
     metrics.nearest_true_idx = nan(size(metrics.true_idx));
     metrics.loc_error = nan(size(metrics.true_idx));
     hit_idx = nan(size(metrics.true_idx));
    miss_idx = metrics.true_idx(isnan(hit_idx));    % Missed true locations
    FA_idx = [];
end
    

    
    %% not using for right now because SIA and MIA do a search first for MCMV_idx and do not want to overide this functionality
    %         if length(peak_hit_idx)>1   % using second cirteria to find hit
    %             d1 = real(h.norm_true_swf(act_samps,v));    % true source waveform
    %             d2 = real(h.current_norm_peak_swf(act_samps,peak_hit_idx)); % hit peak(s) waveform(s)
    %             wave_err = nanmean((bsxfun(@minus,d2,d1)).^2); % simple mean-squared error
    %         else
    %             Hit_idx(v) = peak_hit_idx;
    %
    %         end
    
   metrics.Hits = hit_idx;
   metrics.FA = FA_idx;
   metrics.Miss = miss_idx;
   metrics.num_CR = size(h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,1) - ( length(hit_idx) + length(FA_idx) + length(miss_idx) );

% else
%             metrics.true_idx=[];
%          metrics.dist_thresh=[];
%     metrics.nearest_true_idx=[];
%            metrics.loc_error=[]; 
%                 metrics.Hits=[]; 
%                   metrics.FA=[]; 
%                 metrics.Miss=[];
%                   metrics.num_CR=[];
% end

h.inv_soln(h.current_inv_soln).classifier_metrics = metrics;

% inv_soln's classifier performance
[perf]=sm_calc_classifier_performance(length(metrics.Hits(~isnan(metrics.Hits))),length(metrics.Miss),metrics.num_CR,length(metrics.FA));
if isnan(perf.MCC); perf.MCC=0; end
h.inv_soln(h.current_inv_soln).classifier_performance = perf;

%% printing to screen
h.text_inv_soln_performance.String = sprintf('Hits=%.f; Misses=%.f; FA=%.f; CR=%.f; MCC=%.2f',length(metrics.Hits(~isnan(metrics.Hits))),length(metrics.Miss),length(metrics.FA),metrics.num_CR,perf.MCC);

%% removing basket surround for missed true sources
for v=1:length(h.map3D_true_locs); h.map3D_true_locs(1,v).Visible = 'off'; end
hit_idx = find(~isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.Hits));
for v=1:length(hit_idx); h.map3D_true_locs(1,hit_idx(v)).Visible = 'on'; end




