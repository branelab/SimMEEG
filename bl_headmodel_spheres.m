function bl_hdm = bl_headmodel_spheres(varargin)
% Creates spherical headmodel based on mesh layers (.bnd) shells (tri, pos) to be used with "bst_eeg_sph.m" to get leadfields
% 
% INPUT
%   single input
%       hdm.bnd = mesh layers with .tri and .pos
%       hdm.cond = conductivities of each mesh layer e.g., [0.3300 0.0042 0.3300] for mesh layers in order of scalp, skull, and brain
% or
%   mesh = .tri and .pos (units = 'meters') of layers in headmodel for calculating sphere origins and radii
%   cond = conductivities of each mesh layer e.g., = [0.3300 0.0042 0.3300] for mesh layers in order of scalp, skull, and brain
%
% OUTPUT
%   hdm 



if nargin==1 
    bl_hdm = varargin{1}; 
    bl_hdm.type = 'Spheres';
elseif nargin==2
    bl_hdm.type = 'Spheres';
    bl_hdm.units = 'm';
    bl_hdm.bnd = varargin{1};
    bl_hdm.cond = varargin{2};
else
    warning('Too many inputs'); return
end

%% Following was extracted and modified from "ft_headmodel_concentricspheres.m"
% concatenate the vertices of all surfaces
pos = {bl_hdm.bnd.pos};
pos = cat(1, pos{:});

% remove double vertices
pos  = unique(pos, 'rows');
npos = size(pos, 1);

% fit a single sphere to all combined headshape points, this is used for the center of the spheres
[single_o, single_r] = fitsphere(pos);
% fprintf('initial sphere: total number of unique surface points = %d\n', npos);
% fprintf('initial sphere: center = [%.1f %.1f %.1f]\n', single_o(1), single_o(2), single_o(3));
% fprintf('initial sphere: radius = %.1f\n', single_r);

% fit the radius of a single sphere to the corresponding surface points of each mesh
for i = 1:numel(bl_hdm.bnd)
  npos     = size(bl_hdm.bnd(i).pos,1);
  dist     = sqrt(sum(((bl_hdm.bnd(i).pos - repmat(single_o, npos, 1)).^2), 2));
  bl_hdm.r(i)    = mean(dist);
  bl_hdm.cond(i) = bl_hdm.cond(i);
end

% specify the center of the spheres
bl_hdm.o = single_o;

% sort the spheres from the smallest to the largest
[dum, indx] = sort(bl_hdm.r);

% order the spheres from the smallest to the largest ('insidefirst' order)
bl_hdm.r = bl_hdm.r(indx);
bl_hdm.cond = bl_hdm.cond(indx);

