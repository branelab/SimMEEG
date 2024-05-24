function [SD_PSF, SD_CTF, Rm, nRm, nSD_PSF, nSD_CTF] = sm_calc_SPF_CTF_inv_soln(H,ori,wts,voxel_locs)
%   This function calculates the spatial dispersion of the point-source spread function (PSF) and cross-talk function (CTF) for each voxels within the leadfield matrix (H) for
%   specific inverse solution's weights/ImagingKernel (wts). The scalar leadfield (Hscalar) is calculated 
%
%   Hscalar = l
%
%   Rm = Resolution matrix   R = wts'*Hscalar
%
% calculations are based on those from Hauk et al., 2011 NeuroImage 54(3): 1966-1974

Hscalar = sm_get_scalar_leadfield(H,ori);

if size(wts,3)>1    % vector inv_soln - select wts with maximum orientation 
    [~,m_idx] = max(ori,[],2); 
    wts2 = zeros(size(wts,1), size(wts,2)); 
    for v=1:size(wts,2); wts2(:,v) = wts(:,v,m_idx(v)); end
    wts=wts2; clear wts2;
end

Rm = wts'*Hscalar;  % see Hauk et al., 2011 NeuroImage 54(3): 1966-1974
nRm = Rm/max(max(abs(Rm)));  % normalized resolution matrix

vx_num = size(voxel_locs,1);
% diff_locs = nan(vx_num,vx_num); 
SD_PSF = nan(vx_num,1);
SD_CTF = SD_PSF;
nSD_PSF = SD_PSF;
nSD_CTF = SD_PSF;

n_out = nargout>4;
for v=1:size(voxel_locs,1)
    %     diff_locs(v,:) = sqrt( (voxel_locs(v,1)-voxel_locs(:,1) ).^2 + (voxel_locs(v,2)-voxel_locs(:,2) ).^2 + (voxel_locs(v,3)-voxel_locs(:,3) ).^2);
    diff_locs = sqrt( (voxel_locs(v,1)-voxel_locs(:,1) ).^2 + (voxel_locs(v,2)-voxel_locs(:,2) ).^2 + (voxel_locs(v,3)-voxel_locs(:,3) ).^2);
    SD_PSF(v) = sqrt( (sum ( diff_locs .* Rm(:,v).^2)) ./ sum(Rm(:,v).^2) ) ; % spatial dispersion for point-source function (voxel's point-source spread out to other locations)
    SD_CTF(v) = sqrt( (sum ( diff_locs .* Rm(v,:)'.^2)) ./ sum(Rm(v,:)'.^2) ) ; % spatial dispersion for cross-talk function (other voxel's spread or cross-talk into this voxel's location)
end
if n_out==1
    for v=1:size(voxel_locs,1)
        nSD_PSF(v) = sqrt( (sum ( diff_locs .* nRm(:,v).^2)) ./ sum(nRm(:,v).^2) ) ; % spatial dispersion for point-source function (voxel's point-source spread out to other locations)
        nSD_CTF(v) = sqrt( (sum ( diff_locs .* nRm(v,:)'.^2)) ./ sum(nRm(v,:)'.^2) ) ; % spatial dispersion for cross-talk function (other voxel's spread or cross-talk into this voxel's location)
    end
end
