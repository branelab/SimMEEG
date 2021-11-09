function [h_slices]=bl_plot_slice_bmf(grid_locs,inside_idx,xyz_pts,img,slice_orient,sFaceAlpha,vol_bmf,min_max)
%
%
%
% slice_orient = [0 1 1] = [sagittal coronal axial] ; flags for turning off(0) or on(1) the slices at xyz_pts.
if nargin<6; sFaceAlpha=1; vol_bmf=[]; min_max=[0 1];
elseif nargin<7; vol_bmf=[]; min_max=[0 1];
elseif nargin<8; min_max=[0 1];
end

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

% need to be double precision
grid_locs = double(grid_locs); img=double(img); x=double(x); y=double(y); z=double(z); 

if length(grid_locs)~=numel(x)  % need to find img locations witihin meshgrid
    gd =grid_locs; 
    [img2]= griddata(gd(:,1),gd(:,2),gd(:,3),img,x,y,z,'natural');     % putting img into cubed 3D grid for slice.m
%     img2=permute(img2,[2 1 3]); % puts it into ft_plot_mesh X-->Y dimnesions
    [vx_idx2]=find_nearest_voxel(xyz_pts,grid_locs);  % finding closest dip on mesh net for cortically constrained modeling
    
else
    % creating img that is shaped in same coordinates at x,y,z --> needed for "slice.m".
    img_grid=nan(size(grid_locs,1),1);
    img_grid(inside_idx)=img;   % putting img into the grid img
    img2=reshape(img_grid,[dims(1) dims(2) dims(3)]);
    img2=permute(img2,[2 1 3]); % puts it into ft_plot_mesh X-->Y dimnesions
    [vx_idx2]=find_nearest_voxel(xyz_pts,grid_locs);  % finding closest dip on mesh net for cortically constrained modeling
    
end

if sum(slice_orient)~=0
    if slice_orient(1)==1; xslice=unique(grid_locs(vx_idx2,1)); vw_angle=[-90 0]; else xslice=[]; end  % sagittal
    if slice_orient(2)==1; yslice=unique(grid_locs(vx_idx2,2)); vw_angle=[0 0]; else yslice=[]; end  % coronal
    if slice_orient(3)==1; zslice=unique(grid_locs(vx_idx2,3)); vw_angle=[0 90]; else zslice=[]; end  % axial
    
    s1=slice(x,y,z,img2,xslice,yslice,zslice); %view(0,90);
    h_slices=s1;
    if sFaceAlpha==0
        % setting slice tansparency to interp so that smaller values have higher transparency.
        for t=1:size(s1); s1(t).FaceColor='interp'; s1(t).FaceAlpha='interp'; end% sFaceAlpha*( max_data(t)/max(max_data)); end
        caxis(min_max);
%         alpha_scale=(min_max(1):1e-1:min_max(2))/min_max(2);
%         img_scale=(min_max(1):1e-2:min_max(2))/min_max(2);
%         alpha_scale=zeros(length(0:1e-2:min_max(2)),1);
%         alpha_scale(end-length(img_scale)+1:end)=img_scale;
        alpha_scale=linspace(0,1,1000);
        min_img=min_max(1)/min_max(2);        
        alpha_scale(alpha_scale<min_img)=0;
        alpha('color'); alphamap(alpha_scale);  
        shading interp;
        caxis(min_max)
    else
        for t=1:size(s1); max_data(t)=nanmax(nanmax(abs(s1(t).CData))); end
        for t=1:size(s1); s1(t).FaceAlpha=sFaceAlpha; end %*( max_data(t)/max(max_data)); end
        shading interp;
        caxis(min_max)
    end
    
else
%     vw_angle=[0 90];
end

hold on;
if ~isempty(vol_bmf);    bl_plot_mesh(vol_bmf); end
% view(vw_angle);
