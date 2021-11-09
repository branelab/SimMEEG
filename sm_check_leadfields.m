function leadfield_out = sm_check_leadfields(leadfield)
global h

leadfield_out = leadfield;

if size(leadfield_out.pos,1)<=size(leadfield_out.inside,1)  % pos is likely not a square matrix
    grid_locs = leadfield.pos;
    % setting square grid locations - need for surf.m
    t=1;   min_d=min(grid_locs(:,t)); max_d=max(grid_locs(:,t)); res_d=abs(grid_locs(1:end-1,t)-grid_locs(2:end,t))'; res_d=mode(res_d(res_d~=0)); % voxel resolution
    x_coords=min_d:res_d:max_d;
    t=2;   min_d=min(grid_locs(:,t)); max_d=max(grid_locs(:,t)); res_d=abs(grid_locs(1:end-1,t)-grid_locs(2:end,t))'; res_d=mode(res_d(res_d~=0)); % voxel resolution
    y_coords=min_d:res_d:max_d;
    t=3;   min_d=min(grid_locs(:,t)); max_d=max(grid_locs(:,t)); res_d=abs(grid_locs(1:end-1,t)-grid_locs(2:end,t))'; res_d=mode(res_d(res_d~=0)); % voxel resolution
    z_coords=min_d:res_d:max_d;
    dims=[length(x_coords) length(y_coords) length(z_coords)]; % cubed matrix
    
    % creating x,yz grid
    [x,y,z]=meshgrid(x_coords,y_coords,z_coords);
    x1=reshape(x,numel(x),1); 
    y1=reshape(y,numel(y),1); 
    z1=reshape(z,numel(z),1); 
    xyz = [x1 y1 z1];
    
    [vx_idx2]=find_nearest_voxel(leadfield.pos,xyz);  % finding closest dip on mesh net for cortically constrained modeling
    inside_idx = false(size(xyz,1),1);
    inside_idx(vx_idx2) = true;
    
    lf = cell(size(xyz,1),1);
    leadfield_out.inside = inside_idx;
    lf(inside_idx==1) = leadfield.leadfield'; 
    leadfield_out.leadfield = lf';
    leadfield_out.pos = xyz;
    
end


% figure; clf; hold on; scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'k.','MarkerFaceAlpha',.01); scatter3(xyz(inside_idx,1),xyz(inside_idx,2),xyz(inside_idx,3),'ro','MarkerFaceAlpha',.01);



