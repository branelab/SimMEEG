function [p,s1]=bl_plot_source_locs(vol,source_locs_mm,opt)
%
% see bl_plot_mesh.m for "vol" input options
% INPUT:
%   vol(v) = volumetric data for volumes (v)
%   vol(v).tri = Faces for 'patch.m' function
%   vol(v).pos = Vertices for 'patch.m' function
%   vol(v).img = image map with same vertices as vol(v).tri;
%
%   opt.axes_h = handle of axes to plot headmodel.
%   opt.vol_nums = indices of volumes(v) to plot (default is to plot all volumes).
%   opt.source_clr = source colors; (default = [0 0 1; 1 0 0; 0 0.8 0])
%   opt.source_size = size of source location marker; (default = [100 100 100])
%   
% optional (will be updated with defaults if not present)
%   vol(v).FaceAlpha=.25; 
%   vol(v).FaceColor=[1 .1 .5]*.6;
%   vol(v).EdgeAlpha=.25;
%   vol(v).EdgeColor='none';
%
%
if ~isfield(opt,'source_clr')
    opt.source_clr=[0 0 1; 1 0 0; 0 0.8 0]; % default 
end
if ~isfield(opt,'source_size')
%     opt.source_size=[100 100 100]; % default 
    opt.source_size=repmat(100,size(source_locs_mm,1),1); % default 
end
if ~isfield(opt,'axes_h')
    figure; set(gcf,'color','w'); cla; ax=gca;
    opt.axes_h=ax;
else
    if ~isvalid(opt.axes_h)
        figure; set(gcf,'color','w'); cla; ax=gca;
        opt.axes_h=ax;
    end
end

axes(opt.axes_h); cla;
for v=1:size(source_locs_mm,1)
s1(v)=scatter3(source_locs_mm(v,1),source_locs_mm(v,2),source_locs_mm(v,3),'MarkerFaceColor',opt.source_clr(v,:),'MarkerEdgeColor',opt.source_clr(v,:)*.3,'SizeData',opt.source_size(v));
end
hold on;
p=bl_plot_mesh(vol,opt);
% rotate3d;
