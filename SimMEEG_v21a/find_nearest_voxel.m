function [v_idx]=find_nearest_voxel(xyz_locs,xyz_voxels)
%function [v_idx]=find_nearest_voxel(xyz_locs,xyz_voxels);
%
%   xyz_locs = [# locations x X Y Z];
%   xyz_voxels [# voxels x X Y Z];
%
%   finds the nearest voxel location based on least-squares fit between locations.
v_idx=[];
parfor v=1:size(xyz_locs,1)
    [~,xy]=sort(sum((xyz_voxels-repmat(xyz_locs(v,:),[size(xyz_voxels,1) 1])).^2,2),'ascend');
    v_idx(v)=xy(1);
end




