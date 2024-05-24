function Hscalar = sm_get_scalar_leadfield(H,ori)
% This function get the scalar leadfield (Hscalar) for the leadfield (H [chans x 3 orientations x voxels]) with the scalar orientations (ori).

parfor v=1:size(H,3) %1:length(H_idx); %1:size(H,3);
    lf = H(:,:,v)*ori(v,:)';     % lead field for dipole with best orientation (yields highest signal).
    lf = lf/norm(lf);              % normalizes the leadfield across channels for each voxel (Needed to get proper weights)
    Hscalar(:,v)=lf;             % lead-field for dipole with best orientation (i.e., with largest eigenvalue).
end
