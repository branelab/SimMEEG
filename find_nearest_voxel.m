function [v_idx]=find_nearest_voxel(xyz_locs,xyz_voxels)
%function [v_idx]=find_nearest_voxel(xyz_locs,xyz_voxels);
%
%   xyz_locs = [# locations x X Y Z];
%   xyz_voxels [# voxels x X Y Z];
%   
%   finds the nearest voxel location based on least-squares fit between locations.

v_idx=[];

%% vectorized for speed
for v=1:size(xyz_locs,1)
    dp = repmat(xyz_locs(v,:),size(xyz_voxels,1),1);
    x = xyz_voxels-dp;
    [~,v_idx(v)] = min(sum((x.^2),2));
end



% for v=1:size(xyz_locs,1)
%     [xx,v_idx(v)]=min(sum((xyz_voxels-repmat(xyz_locs(v,:),[size(xyz_voxels,1) 1])).^2,2));
% end

