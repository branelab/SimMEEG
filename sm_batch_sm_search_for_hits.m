function h = sm_batch_sm_search_for_hits(varargin)
 h= varargin{1};

%% found peak voxels
switch varargin{end}
    case 'initial search'
        found_idx = 1:size(h.inv_soln(h.current_inv_soln).peak_voxels,1);
        h.current_3D_peak_idx = h.inv_soln(h.current_inv_soln).peak_voxels(found_idx,5);
        h.current_3D_peak_voxels = h.inv_soln(h.current_inv_soln).peak_voxels(found_idx,:);
        
    case 'slider thresh'
        h.current_3D_thresh = h.slider_3D_image_thresh.Value;
        
        % Thresholding base on slider threshold
        h.current_inv_soln_show_peak_idx = find(h.inv_soln(h.current_inv_soln).peak_voxels(:,4)>h.current_3D_thresh);
        h.current_inv_soln_hide_peak_idx = find(h.inv_soln(h.current_inv_soln).peak_voxels(:,4)<=h.current_3D_thresh);  % peaks to be hidden when below threshold
        found_idx = h.current_inv_soln_show_peak_idx;
        
        h.current_3D_peak_idx = h.inv_soln(h.current_inv_soln).peak_voxels(found_idx,5);
        h.current_3D_peak_voxels = h.inv_soln(h.current_inv_soln).peak_voxels(found_idx,:);
end


%% Get True locs
if isfield(h,'sim_data')
    if isfield(h.sim_data,'cfg')
        true_locs = h.sim_data.cfg.source.vx_locs;
        true_idx = h.sim_data.cfg.source.vx_idx;
    else
        true_locs = h.cfg.source.vx_locs;
        true_idx = h.cfg.source.vx_idx;
    end
else
    true_locs = h.cfg.source.vx_locs;
    true_idx = h.cfg.source.vx_idx;
end
seed_idx = 1:size(true_locs,1);

%% Initializing variables
h.current_peak_hit_lf_idx = nan(1,length(true_idx));
h.current_peak_miss_lf_idx = nan(1,length(true_idx));
h.current_peak_fa_lf_idx = [];

if isempty(found_idx)  % No peaks found
    v_idx = nan(1,length(true_idx));
    peak_idx = nan(1,length(true_idx));
    p_idx = nan(1,length(true_idx));
    diff_idx = [];
else
    %% %%%%% Search Type: find nearest voxel to true source voxel or voxel with best swf match to true swf
    %% finding peak_hit but if nearest voxel is not within user-defined search radius then it gets "nan" = miss
    %             [v_idx]=find_nearest_voxel(vx_locs,h.current_3D_peak_voxels(:,1:3));
    search_rad = (sum(3*(str2num(h.edit_inv_dist_search_thresh.String)).^2)).^.5; % scalar distances = radius
    [norm_swf2, norm_true_swf] = sm_batch_sm_calc_normalized_swf(h,h.current_3D_peak_idx); % need to z-normalize waves because InvSoln have different scales
     norm_swf = nan(size(norm_swf2,1),max(found_idx));
     norm_swf(:,found_idx) = norm_swf2;
   clear diff_locs2 best_mse best_cor flip_idx
    for v=1:size(true_locs,1)
        % loc diff
        diff_locs = true_locs(v,:) - h.current_3D_peak_voxels(:,1:3); % vector loc difference
        diff_locs2(:,v) = sqrt(sum(diff_locs.^2,2));     % scalar loc difference
        % swf diff
        diff_mse = (nanmean((norm_swf(h.inv_soln(h.current_inv_soln).params.act_samps,:)-norm_true_swf(h.inv_soln(h.current_inv_soln).params.act_samps,v)).^2)).^.5;  %
        flip_mse = (nanmean((-norm_swf(h.inv_soln(h.current_inv_soln).params.act_samps,:)-norm_true_swf(h.inv_soln(h.current_inv_soln).params.act_samps,v)).^2)).^.5;  % need to flip because 180deg orientation flips can randomly occur with InvSoln
        [best_mse(:,v),flip_idx(:,v)] = min([diff_mse; flip_mse]); % flip_idx is for flipping swf ori 180deg
        best_mse(diff_locs2(:,v)>search_rad,v) = nan; % only peak_voxels within user-defined hit search radius
        
        best_cor(:,v) = abs(corr(norm_swf(h.inv_soln(h.current_inv_soln).params.act_samps,:),norm_true_swf(h.inv_soln(h.current_inv_soln).params.act_samps,v)));
        best_cor(diff_locs2(:,v)>search_rad,v) = nan; % only peak_voxels within user-defined hit search radius
        
        
        diff_locs2(diff_locs2(:,v)>search_rad,v) = nan; % only peak_voxels within user-defined hit search radius
    end
    flip_idx(flip_idx==2)=-1;
    
    %             best_mse = nan(4,3); best_mse([1 4],2) = nan; best_mse([1 4],1) = [2 7]; best_mse([1 4],3) = [3 8];
    peak_idx = h.current_3D_peak_idx;
    peak_hit = nan(size(seed_idx)); p_idx = peak_hit;
    
    switch h.menu_inv_hit_search_type.String{h.menu_inv_hit_search_type.Value}
        case 'Nearest'
            search_vals = diff_locs2;   % searches for minimum loc error (nearest) within the search radius
            [rw_idx,cl_idx] = sm_search_min_error_matrix(search_vals); % permutes matrices with duplicates removed to find minimum mse
            peak_hit(cl_idx) = peak_idx(rw_idx);
            p_idx(cl_idx) = rw_idx;
        case 'Wave Error'
            search_vals = best_mse;    % searches for best wave error (MSE) within the search radius
            [rw_idx,cl_idx] = sm_search_min_error_matrix(search_vals); % permutes matrices with duplicates removed to find minimum mse
            peak_hit(cl_idx) = peak_idx(rw_idx);
            p_idx(cl_idx) = rw_idx;
        case 'Wave Correlation'
            search_vals = 1-best_cor;    % searches for best wave correlation (finds minimum (1-corr) within the search radius
            [rw_idx,cl_idx] = sm_search_min_error_matrix(search_vals); % permutes matrices with duplicates removed to find minimum mse
            peak_hit(cl_idx) = peak_idx(rw_idx);
            p_idx(cl_idx) = rw_idx;
    end
    
    %%
    diff_idx = setdiff(1:length(peak_idx),p_idx);
    % re-ordering non-hit peaks so that they are largest to smallest for image value
    [~,diff_ord] = sort(h.current_3D_peak_voxels(diff_idx,4),'descend');
    try
        v_idx = [p_idx diff_idx(diff_ord)];
    catch
        v_idx = [p_idx' diff_idx(diff_ord)];
    end
    
    
    
end

h.current_3D_peak_idx = nan(size(v_idx));
peak_voxels = h.current_3D_peak_voxels;
h.current_3D_peak_voxels = nan(size(peak_voxels));

for v = 1:length(v_idx)
    vv_idx = v_idx(v);
    if ~isnan(vv_idx)
        h.current_3D_peak_idx(v) = peak_idx(vv_idx);
        h.current_3D_peak_voxels(v,:) = peak_voxels(vv_idx,:);  % re-ordered
    else
        
    end
end
h.current_3D_peak_voxels(h.current_3D_peak_voxels(:,5)==0,:) = nan; %


h.current_peak_hit_lf_idx(~isnan(p_idx)) = peak_idx(p_idx(~isnan(p_idx))); %leadfield voxel # for inv_soln headmodel with nans for misses
h.current_peak_miss_lf_idx(isnan(p_idx)) = true_idx(isnan(p_idx));  %leadfield voxel # for true source headmodel # with nans for misses
h.current_peak_fa_lf_idx = peak_idx(diff_idx); %leadfield voxel # for inv_soln headmodel



