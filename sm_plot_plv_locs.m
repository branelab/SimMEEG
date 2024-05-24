function sm_plot_plv_locs(varargin)
global h
    seed_idx = h.inv_soln(h.current_inv_soln).plv_seed_idx; seed_idx = seed_idx(~isnan(seed_idx));
    comp_idx = h.inv_soln(h.current_inv_soln).plv_comp_idx; comp_idx = comp_idx(~isnan(comp_idx));
    vx_pos = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos;

if isfield(h.inv_soln(h.current_inv_soln),'plv_seed_idx') && isfield(h,'seed_plot') && isfield(h,'comp_plot')
    % if plots exist then delete them
    delete(h.seed_plot); delete(h.comp_plot);
end

%% plot seed locations
if h.radio_3D_plot_plv_locs.Value==1
%              h.lf_plot = plot3(h.axes_3D_images,vx_pos(:,1),vx_pos(:,2),vx_pos(:,3),'k.','MarkerFaceColor','k');
    h.comp_plot = plot3(h.axes_3D_images,vx_pos(comp_idx,1),vx_pos(comp_idx,2),vx_pos(comp_idx,3),'ko','MarkerSize',4,'MarkerFaceColor','none');
    h.seed_plot = plot3(h.axes_3D_images,vx_pos(seed_idx,1),vx_pos(seed_idx,2),vx_pos(seed_idx,3),'ms','MarkerSize',8,'MarkerFaceColor','none');
end


