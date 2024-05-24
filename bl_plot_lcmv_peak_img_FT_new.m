function [peak_voxels,vx_idx,s1,p1,s2,p2,s3]=bl_plot_lcmv_peak_img_FT_new(haxis, img,ori,thresh_val,dist_limit,vx_locs,cmap,min_max,vol,elec,vw_angle,pk_flag,vol_types,grid_locs,inside_idx,slice_orient,sFaceAlpha)
%function [peak_voxels,vx_idx,s1]=bl_plot_lcmv_peak_img_FT(img,ori,thresh_val,dist_limit,vx_locs,cmap,min_max,headmodel,sens,pk_flag);
%
% INPUT
%   img = amplitude values for each voxel [M x1]
%	thresh_val = amplitude threshold to mask img before searching for peak.
%	thresh_limit = minimum distance (in mm) between peaks
%   vx_locs = grid locations for voxels (in mm) [M x 3]
%   cmap = colormap
%   min_max = map scale for cdata
if nargin<10; pk_flag=0; end
% fprintf('\nPeak Searching Parameters:\n   Threshold value = %.2f',thresh_val);
% fprintf('\n   Minimum search distance = %.f mm\n',dist_limit);
s1=[];p1=[]; s2=[]; s3 = []; p2 =[];
haxis.NextPlot = 'replace';

voxel_vals=[vx_locs, img];
[peak_voxel,p_idx]=BRANELab_find_peak_voxel_thresh(voxel_vals,thresh_val,dist_limit);
if ~isempty(peak_voxel)
    peak_voxels=flipud(sortrows([peak_voxel p_idx],4));
    vx_idx=peak_voxels(:,5);
    if length(thresh_val)==1;   % single-state
        sig_idx=find(voxel_vals(:,4)>thresh_val);
        nan_idx=find(voxel_vals(:,4)<thresh_val);
    elseif length(thresh_val)==2;   % dual-state
        sig_idx1=find(voxel_vals(:,4)>thresh_val(2));
        sig_idx2=find(voxel_vals(:,4)<thresh_val(1));
        nan_idx=find(voxel_vals(:,4)<thresh_val(2) & voxel_vals(:,4)>thresh_val(1));
        sig_idx=[sig_idx1;sig_idx2];
    end
    voxel_vals(nan_idx,4)=nan;
    voxel_vals(voxel_vals(:,4)==0,4)=nan;
    ns=size(peak_voxels,1);
    ln_clr=lines(size(peak_voxels,1));
    %ln_clr=zeros(size(peak_voxels,1),3);   % for plotting Black & White images.
    %voxel_vals(:,1:3)=voxel_vals(:,1:3)/1000;
else
    peak_voxels=[]; vx_idx=[];
end
for vw=1%:3;
    haxis.CLim = min_max;                                            
    % s1=scatter3(haxis,voxel_vals(nan_idx,1),voxel_vals(nan_idx,2),voxel_vals(nan_idx,3),'.','markeredgecolor',[.7 .7 .7]);
    %s1=scatter3(haxis,voxel_vals(:,1),voxel_vals(:,2),voxel_vals(:,3),'.','Cdata',voxel_vals(:,4));
    if ~isempty(peak_voxel)
        if pk_flag==0  % thresholded map
            % plotting slices at found voxel locations
            xyz_pts=peak_voxel(:,1:3);
            vol_bmf=[]; %vol.bnd(3);
            if sum(slice_orient)>0; [s2]=bl_plot_slice_bmf(haxis, grid_locs,inside_idx,xyz_pts,img,slice_orient,sFaceAlpha,vol_bmf,min_max);             end
            
        elseif pk_flag==1  % peaks + thresholded map
            %sig_idx=setdiff(sig_idx,peak_voxels(:,5));
            %             s2=scatter3(haxis,voxel_vals(sig_idx,1),voxel_vals(sig_idx,2),voxel_vals(sig_idx,3),'s','filled','Cdata',voxel_vals(sig_idx,4));
            
            % plotting slices at found voxel locations
            xyz_pts=peak_voxel(:,1:3);
            vol_bmf=[]; %vol.bnd(3);
            if sum(slice_orient)>0; [s2]=bl_plot_slice_bmf(haxis,grid_locs,inside_idx,xyz_pts,img,slice_orient,sFaceAlpha,vol_bmf,min_max);             end
            haxis.NextPlot = 'add';
            for v=1:ns
                % changes marker size and orientation length to match amplitude of voxel.
                amp_gain=abs((voxel_vals(peak_voxels(v,5),4)/min_max(2)));
                s1(v)=scatter3(haxis,voxel_vals(peak_voxels(v,5),1),voxel_vals(peak_voxels(v,5),2),voxel_vals(peak_voxels(v,5),3),'o','MarkerEdgeColor',ln_clr(v,:),'MarkerFaceColor','none','sizedata',amp_gain*50,'linewidth',2);
                %keyboard;
                % plotting orientations
                vx_res=max(vx_locs(2,:)-vx_locs(1,:));
                vx_pos=vx_locs(peak_voxels(v,5),:);
                ori_pos=vx_pos+(2*amp_gain*vx_res*ori(peak_voxels(v,5),:));
                p1(v)=plot3(haxis,[vx_pos(1) ori_pos(1)],[vx_pos(2) ori_pos(2)],[vx_pos(3) ori_pos(3)],'color',ln_clr(v,:),'linewidth',2);
            end
            
        elseif pk_flag==2 % peaks only
            for v=1:ns
                % changes marker size and orientation length to match amplitude of voxel.
                amp_gain=abs((voxel_vals(peak_voxels(v,5),4)/min_max(2)));
                s1(v)=scatter3(haxis,voxel_vals(peak_voxels(v,5),1),voxel_vals(peak_voxels(v,5),2),voxel_vals(peak_voxels(v,5),3),'o','MarkerEdgeColor',ln_clr(v,:),'MarkerFaceColor','none','sizedata',amp_gain*50,'linewidth',2);
                %keyboard;
                % plotting orientations
                vx_res=max(vx_locs(2,:)-vx_locs(1,:));
                vx_pos=vx_locs(peak_voxels(v,5),:);
                ori_pos=vx_pos+(2*amp_gain*vx_res*ori(peak_voxels(v,5),:));
                p1(v)=plot3(haxis,[vx_pos(1) ori_pos(1)],[vx_pos(2) ori_pos(2)],[vx_pos(3) ori_pos(3)],'color',ln_clr(v,:),'linewidth',2);
            end
        elseif pk_flag==3   % peaks only with color as magnitude
            s1=scatter3(haxis,voxel_vals(sig_idx,1),voxel_vals(sig_idx,2),voxel_vals(sig_idx,3),'o','filled','Cdata',voxel_vals(sig_idx,4),'MarkerEdgeColor','k','sizedata',100);
            
        elseif pk_flag==4  % thresholded map but midline slices only
            for v=1:ns
                % changes marker size and orientation length to match amplitude of voxel.
                amp_gain=abs((voxel_vals(peak_voxels(v,5),4)/min_max(2)));
                s1(v)=scatter3(haxis,voxel_vals(peak_voxels(v,5),1),voxel_vals(peak_voxels(v,5),2),voxel_vals(peak_voxels(v,5),3),'o','MarkerEdgeColor',ln_clr(v,:),'MarkerFaceColor','none','sizedata',amp_gain*50,'linewidth',2);
                hold on;
                % plotting orientations
                vx_res=max(vx_locs(2,:)-vx_locs(1,:));
                vx_pos=vx_locs(peak_voxels(v,5),:);
                ori_pos=vx_pos+(2*amp_gain*vx_res*ori(peak_voxels(v,5),:));
                p1(v)=plot3(haxis,[vx_pos(1) ori_pos(1)],[vx_pos(2) ori_pos(2)],[vx_pos(3) ori_pos(3)],'color',ln_clr(v,:),'linewidth',2);
            end
            
            % plot slices at peak locations           
            xyz_pts=peak_voxel(:,1:3);
            vol_bmf=[]; %vol.bnd(3);
            project_flag = 0; % priject image to midline slices
            if sum(slice_orient)>0; [s2]=bl_plot_slice_bmf(haxis,grid_locs,inside_idx,xyz_pts,img,slice_orient,sFaceAlpha,vol_bmf,min_max,project_flag);             end

            % projecting data to midline slices
            mid_x = min(grid_locs(:,1)) + 0.5*(range(grid_locs(:,1))); 
            mid_y = min(grid_locs(:,2)) + 0.5*(range(grid_locs(:,2))); 
            mid_z = min(grid_locs(:,3)) + 0.5*(range(grid_locs(:,3))); 
            xyz_pts = [mid_x mid_y mid_z];

            vol_bmf=[]; %vol.bnd(3);
            project_flag = 1; % priject image to midline slices
            if sum(slice_orient)>0; [s3]=bl_plot_slice_bmf(haxis,grid_locs,inside_idx,xyz_pts,img,slice_orient,sFaceAlpha,vol_bmf,min_max,project_flag);             end
        end
        
    elseif isempty(peak_voxel) && pk_flag==0   % thresholded map only
        s1=scatter3(haxis,voxel_vals(sig_idx,1),voxel_vals(sig_idx,2),voxel_vals(sig_idx,3),'s','filled','Cdata',voxel_vals(sig_idx,4),'sizedata',50);
        
        %for v=1:ns; scatter3(haxis,voxel_vals(peak_voxels(v,5),1),voxel_vals(peak_voxels(v,5),2),voxel_vals(peak_voxels(v,5),3),'o','MarkerEdgeColor',ln_clr(v,:),'MarkerFaceColor',ln_clr(v,:),'sizedata',10,'linewidth',2); end
        if vw==2; title(sprintf('Voxel # %s',sprintf('%.f ',[vx_idx])),'color','k'); end
    elseif isempty(peak_voxel) && pk_flag==4   % map only midline
        
        % projecting data to midline slices
        mid_x = min(grid_locs(:,1)) + 0.5*(range(grid_locs(:,1)));
        mid_y = min(grid_locs(:,2)) + 0.5*(range(grid_locs(:,2)));
        mid_z = min(grid_locs(:,3)) + 0.5*(range(grid_locs(:,3)));
        xyz_pts = [mid_x mid_y mid_z];
        
        vol_bmf=[]; %vol.bnd(3);
        project_flag = 1; % priject image to midline slices
        if sum(slice_orient)>0; [s3]=bl_plot_slice_bmf(haxis,grid_locs,inside_idx,xyz_pts,img,slice_orient,sFaceAlpha,vol_bmf,min_max,project_flag);             end
        
    end
    
    
    if nargin>7
        if vol_types==1; %sum(vol_types==1)>0
            % Scalp
            %                  vol=vol.bnd(1); opt.vol_nums=1; vol.FaceColor=[.8 .6 .6]; vol.FaceAlpha=0.3; vol.EdgeColor='none'; bl_plot_mesh(vol,opt); bl_plot_mesh(vol)
            p2=bl_plot_mesh(vol,[],haxis);
            %             scalp_clr=[1 1 1]*.8; ft_plot_mesh(vol.bnd(1), 'facecolor',scalp_clr, 'facealpha', sFaceAlpha*.25, 'edgecolor', 'none', 'edgealpha', 'none'); % scalp
        end
        if vol_types==3; %sum(vol_types==3)>0
            % Brain Hull
            %              vol=vol.bnd(3); opt.vol_nums=1; vol.FaceColor=[.6 .6 .6]; vol.FaceAlpha=0.3; vol.EdgeColor='none'; bl_plot_mesh(vol,opt); bl_plot_mesh(vol)
            p2=bl_plot_mesh(vol,[],haxis);
            %             brain_clr=[1 1 1]*.8; ft_plot_mesh(vol.bnd(3),'facecolor',brain_clr,'facealpha', sFaceAlpha*.25,'edgecolor', 'none', 'edgealpha', 0.25); % brain
        end
        if ~isempty(elec);
            % Electrodes
            %         e_clr=[0 .2 .2];ft_plot_sens(elec,'label', 'label','fontcolor',e_clr,'elecshape','point','elecsize',5,'facecolor',e_clr);
        end
    end
    %     view(-90,90); axis off; axis tight;
    %     view(vw_angle(vw,:)); colormap(cmap); %colorbar;
end

view(haxis,vw_angle); haxis.Colormap = cmap;

