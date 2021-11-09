function [peak_voxel,p_idx]=BRANELab_find_peak_voxel_thresh(voxel_vals,thresh_val,thresh_limit);
%function [peak_voxel,p_idx]=BRANELab_find_peak_voxel_thresh(voxel_vals,thresh_val,thresh_limit);
%
% This program will step through the voxels and find the peak absolute  amplitude
%	within a region separated by the thresh_limit value
%
% INPUT:
%   voxel_vals = xyz coordinates of each voxel (M x N = voxels x [x y z amp ...] ) must be at least an M x 4 matrix
%   thresh_val = value of voxel amplitude to threshold data before peak searching
%   thresh_limit = minimum distance between peak voxels (distance thus depends on voxel step size). This creates a bounding box with thresh limit distance.
%
% OUTPUT:
%   [peak_voxel] = peak voxels with all the columns in original voxel_vals
%
% written by A. Herdman June 30,2005
%
%   modified by A. Herdman June 6, 2015 - updated description to make more sense. 
%

voxel_abs=voxel_vals;
voxel_abs(:,4)=abs(voxel_vals(:,4));

if length(thresh_val)==1;
    [lv_idx]=find(voxel_abs(:,4)>thresh_val);
elseif length(thresh_val)==2;
    [lv_idx1]=find(voxel_abs(:,4)>thresh_val(2));
    [lv_idx2]=find(voxel_abs(:,4)<thresh_val(1));
    lv_idx=[lv_idx1;lv_idx2];
end

p_idx=[];
if isempty(lv_idx)~=1
    
    voxel_lv=voxel_abs(lv_idx,:);
    voxel_org=voxel_vals(lv_idx,:);
    
    p=0;
    
    for n=1:size(voxel_lv,1)
        vox_roi=voxel_lv(n,:);
        vox_roi_org=voxel_org(n,:);
        
        [roi_idx,z,zz]=find((voxel_lv(:,1)< vox_roi(1,1)+thresh_limit & voxel_lv(:,1)>vox_roi(1,1)-thresh_limit)...
            & (voxel_lv(:,2)< vox_roi(1,2)+thresh_limit & voxel_lv(:,2)>vox_roi(1,2)-thresh_limit)...
            & (voxel_lv(:,3)< vox_roi(1,3)+thresh_limit & voxel_lv(:,3)>vox_roi(1,3)-thresh_limit));
        
        if vox_roi(1,4)==max(voxel_lv(roi_idx,4))
            p=p+1;
            %fprintf(1,'Number of peaks = %.f\n',p)
            peak_voxel(p,:)=vox_roi_org;
            p_idx=[p_idx n];
        end
    end
else
    peak_voxel=[];
end
p_idx=lv_idx(p_idx);
