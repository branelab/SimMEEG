function sm_plot_headmodel(varargin)
global h; 

xidx = findobj(h.axes_anatomy.Children,'Type','Patch'); num_patches = length(xidx); % existing patches
hdm_clr = [1 0 1];
axes(h.axes_anatomy); hold on;
if num_patches>5 && h.slider_transparency_hdm.Value==0
    msgbox(sprintf('Head Model is already plotted but HDM transparency is ZERO\n\nChange HDM transparency to see Head Model'));
elseif num_patches>5 
    msgbox(sprintf('Head Model is already plotted\n\nIf you cannot see it then try "plot 3D" followed by "plot HeadModel" again'));
else
if strcmp(h.anatomy.headmodel.type,'localspheres')
    btn = questdlg(sprintf('\nPlotting HeadModel\n\nThis may take some time for MEG localspheres\n\nClick "OK" to plot'),'Plot Head Model?','OK','Cancel','Cancel');
else
    btn = 'OK';
end
switch btn
    case 'OK'
        ft_plot_headmodel(h.anatomy.headmodel,'FaceColor',hdm_clr,'FaceAlpha',.1,'EdgeColor','none','EdgeAlpha',.01);
        xidx = findobj(h.axes_anatomy.Children,'Type','Patch');
        h.hdm_patch = xidx(1:end-num_patches); % the last patch one added
        set_transparency([],[],'hdm');
    case 'Cancel'
end
end