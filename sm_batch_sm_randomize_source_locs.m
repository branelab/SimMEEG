function h = sm_batch_sm_randomize_source_locs(h)
% This function randomizes source locations (and orientation if "Random")


num_sources = length(h.monte_params.cfg.source.vx_idx);
dist_thresh = str2num(h.edit_randomize_locs_dist_thresh.String);
rn = randperm(size(h.anatomy.leadfield.H,3));
seed_idx = rn(1);   % first voxel as seed then find other rnadom voxels that are > dist_thresh
vx_pos = h.anatomy.leadfield.voxel_pos;
fprintf('Locations randomized with spatial distances of %.f mm among sources\n',dist_thresh); 


k=0;
hm = waitbar(0,sprintf('Number of Distant Sources Found (Total = %.f)',num_sources));
%% Randomizing Locations
while k<1000
    k=k+1;
    rn = randperm(size(h.anatomy.leadfield.H,3));
    [vdist, ~]= pdist2(vx_pos(rn(1:num_sources),:),vx_pos(rn(1:num_sources),:),'euclidean','Smallest',2); % finds smallest Euclidean distane among sources
    waitbar(sum(vdist(2,:)>dist_thresh)/num_sources, hm)
    num_found = sum(vdist(2,:)>dist_thresh);
%     fprintf('Distant Sources Found = %.f\n',num_found);
    %% break loop
    if  num_found == num_sources
        vx_idx = rn(1:num_sources);
        break;
    end
    if k>=1000
        vx_idx = rn(1:num_sources);
        warndlg(sprintf('Only %.f/%.f sources were found to be %.f mm distant from each other\n\n Try reducing "Number Sources" or "Min Distance"',num_found,num_sources,dist_thresh));
    end
end
h.monte_params.cfg.source.vx_idx = vx_idx;
h.monte_params.cfg.source.vx_locs = h.anatomy.leadfield.voxel_pos(h.monte_params.cfg.source.vx_idx,:);

%% Randomizing Orientations
switch h.menu_ori_normal.String{h.menu_ori_normal.Value}
    case 'Random'
        fprintf('Orientiations set as random.\n');
        h.monte_params.cfg.source.vx_ori =[];
        for v=1:length(h.monte_params.cfg.source.vx_idx)
            az_el = deg2rad(randi([0 360],1,2)); [x,y,z] = sph2cart(az_el(1),az_el(2),1);
            h.monte_params.cfg.source.vx_ori(v,:) = [x y z];  % source orientations (X, Y, Z)
        end
    case 'Cortical Surface' % setting to normal to cortical surface will be done when calling "h.fcn_edit_source_CallBack"
        fprintf('Orientiations fit normal to cortical surface.\n');
        h.anatomy.leadfield.ori = normals(h.anatomy.leadfield.voxel_pos, h.anatomy.mesh_volumes(4).tri);
        for v=1:length(h.monte_params.cfg.source.vx_idx)
            h.monte_params.cfg.source.vx_ori(v,:) = h.anatomy.leadfield.ori(h.monte_params.cfg.source.vx_idx(v),:);
            [az,el] = cart2sph(h.monte_params.cfg.source.vx_ori(v,1),h.monte_params.cfg.source.vx_ori(v,2),h.monte_params.cfg.source.vx_ori(v,3));
            h.edit_source_ori(v).String = sprintf('%.f %.f',rad2deg(az),rad2deg(el));
        end
end

close(hm);
for v=1:3
    h.edit_source_locs(v).String = sprintf('%.f %.f %.f',h.monte_params.cfg.source.vx_locs(v,:));
    [az,el]=cart2sph(h.monte_params.cfg.source.vx_ori(v,1),h.monte_params.cfg.source.vx_ori(v,2),h.monte_params.cfg.source.vx_ori(v,3));
    h.edit_source_ori(v).String = sprintf('%.f %.f',rad2deg(az),rad2deg(el));
end

h.monte_params.cfg.study.source_locs_mm = h.anatomy.leadfield.voxel_pos(h.monte_params.cfg.source.vx_idx,:);
h.monte_params.cfg.study.source_locs = ft_warp_apply(inv(h.anatomy.mri.transform),h.monte_params.cfg.study.source_locs_mm);

% h.fcn_plot_3D_mri(''); 
% h.fcn_edit_source_CallBack(''); 
