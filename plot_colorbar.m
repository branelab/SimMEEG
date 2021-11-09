function plot_colorbar(h,y_scale,cmap,horiz_orient)
%function plot_colorbar(h,y_scale,cmap,horiz_orient)
%   h = handle for axes
%   y_scale = scale for color bar y_scale=[max min]
%   cmap = colormap style = 'bone' 'jet' 'hot' ...
%   horiz_orient = (1) horizontal colorbar

if nargin < 4; horiz_orient=0;end

if horiz_orient==1;
img_h = image([y_scale(1) y_scale(2)],[0,0],flipud(1:size(cmap,1))); 
cb_ytick=get(h,'Xtick');
set(h,'XtickLabel',num2str(flipud(cb_ytick')),'Ytick',[]);
colormap(cmap);
else
img_h = image([0,0],[y_scale(2) y_scale(1)],flipud(1:size(cmap,1))'); 
cb_ytick=get(h,'Ytick');
set(h,'YtickLabel',num2str(flipud(cb_ytick')),'Xtick',[]);
colormap(cmap);
end
  