% running TRAP MUSIC on SimMEEG data

[R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(h.sim_data.sens_final,h.cfg.study.bl_bmf.act_samps, h.cfg.study.bl_bmf.ctrl_samps);
vx_pos = h.anatomy.leadfield.voxel_pos; 

%% creating grid "vx_grid" as per Alex's description "nVoxels = nX * nY * nZ, and the ordering is such that the last index (i.e. Z direction) changes most fast, and the first (i.e. X) - most slow"
vx_res = 5;
xg = min(vx_pos(:,1))-vx_res:vx_res:max(vx_pos(:,1))+vx_res;
yg = min(vx_pos(:,2))-vx_res:vx_res:max(vx_pos(:,2))+vx_res;
zg = min(vx_pos(:,3))-vx_res:vx_res:max(vx_pos(:,3))+vx_res;

% C++ box order
v=0; vx_grid=[];
for d1=1:length(xg)
    for d2=1:length(yg)
        for d3=1:length(zg)
            v=v+1;
           vx_grid(v,:) = [xg(d1) yg(d2) zg(d3)];
        end
    end
end

% matlab box order
v=0; org_grid=[];
for d1=1:length(zg)
    for d2=1:length(yg)
        for d3=1:length(xg)
            v=v+1;
           org_grid(v,:) = [xg(d3) yg(d2) zg(d1)];
        end
    end
end

dims = [length(xg) length(yg) length(zg)]; 
% iV = idx3Dto1D(vx_grid, dims, 1)
v_idx = find_nearest_voxel(vx_pos,vx_grid);
grid_idx = zeros(size(vx_grid,1),1);
grid_idx(v_idx)=1; lstFlag = grid_idx;

%% Leadfield formatting
H = permute(h.anatomy.leadfield.H,[3 2 1]);
arrH = zeros(size(grid_idx,1),size(H,2),size(H,3));
arrH(v_idx,:,:) = H; 
% figure(999); clf; hold on; scatter3(vx_grid(grid_idx==0,1),vx_grid(grid_idx==0,2),vx_grid(grid_idx==0,3),'k.'); 
% scatter3(vx_grid(grid_idx==1,1),vx_grid(grid_idx==1,2),vx_grid(grid_idx==1,3),'ro');


maxSrc = 5; gap = 2; 
% subspace MCMV
sInSMCMV.beamType = 'MER';
sInSMCMV.R = R;
sInSMCMV.arrN = N;
sInSMCMV.arrH = arrH;
sInSMCMV.lstFlag = lstFlag;
sInSMCMV.dims = dims;
sInSMCMV.nSrc = maxSrc;
sInSMCMV.Cavg = Rbar;        % This is C-bar - only needed for MER
sInSMCMV.bMCMVS = true;    % Should be ALWAYS set to true
sInSMCMV.pVal = 1;        % p-value (kind of) for the peak to be considered real. Set to 1 to get all peaks
sInSMCMV.gap = gap;
sInSMCMV.bVerbose = true;
sInSMCMV.bPlotLambda = false;
sInSMCMV.bRAPBeam = false;   % Choose RAP (true) or SMCMV (false). You can try both
sInSMCMV.bDoTRAP = false;             % Don't ask - just set to false :)
sOutSMCMV = doSMCMV(sInSMCMV);

% TRAP MUSIC
sInTRAP.R = R;            % Full covariance
sInTRAP.arrN = N;    % Noise covariance
sInTRAP.arrH = arrH;    % Lead fields
sInTRAP.lstFlag = lstFlag; % Flags
sInTRAP.dims = dims;    % ROI dimensions in voxels
sInTRAP.nSrc = maxSrc;    % Max number of sources to extract
sInTRAP.gap = gap;        % Min allowed number nodes between the peaks (default 2)
sInTRAP.bPlotLambda = false;    % Set to false
sInTRAP.bVerbose = true;    % Enables more printouts

sOutTRAP = trapMUSIC(sInTRAP);

% reversing coordinatesback to anatomy
% rev_idx = find_nearest_voxel(vx_grid(grid_idx==1,:),vx_pos);

%% plotting
% img1 = sOutSMCMV.fImg(1,grid_idx==1)';
% img = nan(length(rev_idx),1);
% img(rev_idx)=img1;
% pk_flag=0; % just plotting the maps
% sFaceAlpha = .2;
% ori = zeros(length(find(grid_idx==1)),3);
% thresh_val = 0; 
% dist_limit = 15; 
% vx_locs = vx_pos;
% cmap = jet(255);
% min_max = [min(img) max(img)];
% vol = h.anatomy.mesh_volumes(3); 
% elec=[]; 
% vw_angle=[0 90];
% vol_types=1;
% grid_locs = org_grid; 
% slice_orient = [1 1 1];
% inside_idx = find(grid_idx==1); 
% [peak_voxels,vx_idx,s1,p1,s2,p2]=bl_plot_lcmv_peak_img_FT_new(img,ori,thresh_val,dist_limit,vx_locs,cmap,min_max,vol,elec,vw_angle,pk_flag,vol_types,grid_locs,inside_idx,slice_orient,sFaceAlpha);


%% just finding max peaks in images 
mcmv_voxel=[]; trap_voxel=[];
for v=1:maxSrc
    [~,pidx] = max(sOutSMCMV.fImg(v,grid_idx==1)');
    mcmv_voxel(v,:) = vx_grid(inside_idx(pidx),:);
    
    [~,pidx] = max(sOutTRAP.fImg(v,grid_idx==1)');
    trap_voxel(v,:) = vx_grid(inside_idx(pidx),:);
end

true_idx = h.cfg.source.vx_idx;
MIA_idx = h.inv_soln(2).soln.MCMV_idx;

figure(998); clf; set(gcf,'color','w'); hold on; 
opt.vol_nums=1; vol = h.anatomy.mesh_volumes(3); 
[p1]=bl_plot_mesh(vol,opt); view(-90,90);
scatter3(vx_pos(true_idx,1),vx_pos(true_idx,2),vx_pos(true_idx,3),'ks','SizeData',220,'linewidth',2); 
scatter3(vx_pos(MIA_idx,1),vx_pos(MIA_idx,2),vx_pos(MIA_idx,3),'ro','filled','SizeData',50); 
s1=scatter3(mcmv_voxel(:,1),mcmv_voxel(:,2),mcmv_voxel(:,3),'bo','filled','SizeData',50); 
scatter3(trap_voxel(:,1),trap_voxel(:,2),trap_voxel(:,3),'o','MarkerFaceColor',[0 .6 0],'MarkerEdgeColor','none','SizeData',50); 
legend({'' 'True' 'MIA' 'sMCMV' 'TRAP'},'Orientation','horizontal','Location','southoutside'); 


