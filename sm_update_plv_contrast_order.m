function sm_update_plv_contrast_order(varargin)
global h;

%% setting current contrasts for plotting and selecting TFR/PLVs
h.current_3D_plv_contrasts = h.inv_soln(h.current_inv_soln).plv_contrasts;
plv_seed_contrasts = nchoosek(h.inv_soln(h.current_inv_soln).plv_seed_idx,2);
x = h.current_3D_plv_contrasts; x(isnan(x)) = 0;
y = plv_seed_contrasts; y(isnan(y)) = 0;
h.current_3D_plv_contrasts_seed_idx =[];
for v=1:size(plv_seed_contrasts,1)
    if sum(y(v,:)==0)==0    % both indices are present
        h.current_3D_plv_contrasts_seed_idx(v) =  find(ismember(x,y(v,:),'rows'));
    else
%         h.current_3D_plv_contrasts_seed_idx(v) = v; % default is to take first 3 - but this is becuaase there are no voxel comparisons found in plv_contrasts
         xy =  find(ismember(x,y(v,:),'rows'));
         h.current_3D_plv_contrasts_seed_idx(v) = xy(1);
    end
end
h.inv_soln(h.current_inv_soln).plv_seed_idx_in_plv_contrast_idx = h.current_3D_plv_contrasts_seed_idx;

%         h.current_3D_plv_contrasts_seed_idx = find(h.current_3D_plv_contrasts_seed_idx==1);
% re-order so that seed_plv_idx are top 3
comp_idx = setdiff(1:size(h.current_3D_plv_contrasts,1),h.current_3D_plv_contrasts_seed_idx);
try
    h.current_3D_plv_contrasts_listbox_order = [h.current_3D_plv_contrasts_seed_idx; comp_idx'];
catch
    h.current_3D_plv_contrasts_listbox_order = [h.current_3D_plv_contrasts_seed_idx  comp_idx]';
end
update_listbox_plv_contrasts();



