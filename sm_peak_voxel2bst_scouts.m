function sm_peak_voxel2bst_scouts
% scouts
xscout = struct('Vertices', [],'Seed', [],'Color', [],'Label', [],'Function', [],'Region', [],'Handles', []);

%% Seed Voxels - True locations
vx_idx = h.sim_data.cfg.source.vx_idx;
vv=0;
for v=1:length(vx_idx) 
    vv=vv+1;
    xscout(vv).Vertices = vx_idx(v); 
    xscout(vv).Seed = vx_idx(v); 
    xscout(vv).Label = sprintf('Seed %.f',v);
    xscout(vv).Function = 'Mean';
    xscout(vv).Region = ''; 
    xscout(vv).Handles = [];
    
end

%% Comparison Voxels - PLV contrasts
v_idx = h.inv_soln(h.current_inv_soln).plv_comp_idx;


[vx_idx] = find_nearest_voxel(h.anatomy.leadfield_eeg_vol(3).voxel_pos(v_idx,:), h.anatomy.leadfield_eeg_cortex(3).voxel_pos); 



for v=1:length(vx_idx) 
    vv=vv+1;
    xscout(vv).Vertices = vx_idx(v); 
    xscout(vv).Seed = vx_idx(v); 
    xscout(vv).Label = sprintf('Comp %.f',v);
    xscout(vv).Function = 'Mean';
    xscout(vv).Region = ''; 
    xscout(vv).Handles = [];
    
end



