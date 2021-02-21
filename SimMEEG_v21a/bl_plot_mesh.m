function p=bl_plot_mesh(vol,opt)
% function bl_plot_mesh(vol,opt);
% This plots the vol data from field trip using 'patch'. It also removes
% dependencies on Field Trip's tools that throw errors due to% compatibility 
% issues with common built-in matlab programs. This program is tested only
% on MATLAB versions R2017b.
% 
% INPUT:
%   vol(v) = volumetric data for volumes (v)
%   vol(v).tri = Faces for 'patch.m' function
%   vol(v).pos = Vertices for 'patch.m' function
%   vol(v).img = image map with same vertices as vol(v).tri;
%
%   opt.vol_nums = indices of volumes(v) to plot (default is to plot all volumes).
%   opt.caxis = colorscale of image
%   
% optional (will be updated with defaults if not present)
%   vol(v).FaceAlpha=.25;   % use 'interp' to scale transparency with the color scale.
%   vol(v).FaceColor=[1 .1 .5]*.6;
%   vol(v).EdgeAlpha=.25;
%   vol(v).EdgeColor='none';
%
%
% written by Anthony Herdman July 20, 2018
%

% set default 'opt' options
if nargin<2 
    opt.vol_nums=1:length(vol);     % indices of volumes within vol to plot
end
% checking for options and setting default if they don't exist
for v=1:length(vol)
    if ~isfield(vol(v),'FaceAlpha') || isempty(vol(v).FaceAlpha); vol(v).FaceAlpha=.25; end
    if ~isfield(vol(v),'FaceColor') || isempty(vol(v).FaceColor); vol(v).FaceColor=[1 1 1]*.6; end
    if ~isfield(vol(v),'EdgeAlpha') || isempty(vol(v).EdgeAlpha); vol(v).EdgeAlpha=.25; end
    if ~isfield(vol(v),'EdgeColor') || isempty(vol(v).EdgeColor); vol(v).EdgeColor='none'; end
end


% read in 'opt' options
for v=opt.vol_nums
    if isfield(vol(v),'img') && ~isempty(vol(v).img)
        if size(vol(v).pos,1)==length(vol(v).img)
            p(v)=patch('Faces',vol(v).tri,'Vertices',vol(v).pos,'FaceVertexCData',vol(v).img,'FaceColor','interp','FaceAlpha',vol(v).FaceAlpha,'EdgeColor',vol(v).EdgeColor,'EdgeAlpha',vol(v).EdgeAlpha);
            
            if isfield(opt,'min_max')
                min_max=opt.min_max;
            else
                min_max=[min(vol(v).img) max(vol(v).img)];
            end
            
            if ~isnumeric(vol(v).FaceAlpha)
                p(v).FaceColor='interp'; p(v).FaceAlpha='interp';
                alpha_scale=linspace(0,1,1000);
                min_img=min_max(1)/min_max(2);
                alpha_scale(alpha_scale<min_img)=0;
                alpha('color'); alphamap(alpha_scale);
                shading interp;
                caxis(min_max)
            end
            
        else
            fprintf('ERROR! Size of img and vertices/pos do NOT match. Only plotting vol data\n');
            p(v)=patch('Faces',vol(v).tri,'Vertices',vol(v).pos,'FaceColor',vol(v).FaceColor,'FaceAlpha',vol(v).FaceAlpha,'EdgeColor',vol(v).EdgeColor,'EdgeAlpha',vol(v).EdgeAlpha);
        end
    else
        p(v)=patch('Faces',vol(v).tri,'Vertices',vol(v).pos,'FaceColor',vol(v).FaceColor,'FaceAlpha',vol(v).FaceAlpha,'EdgeColor',vol(v).EdgeColor,'EdgeAlpha',vol(v).EdgeAlpha);
    end
end
axis off;
material dull; 

