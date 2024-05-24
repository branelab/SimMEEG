function sm_plot_leadfield_grid(varargin)
global h;

if size(h.anatomy.leadfield.voxel_pos,1)>40000
    msgbox(sprintf('Number of lead field voxels > 40,000\n\nToo many to plot\n'));
else
    xidx = findobj(h.axes_anatomy.Children,'Type','Patch'); num_patches = length(xidx); % existing patches
    lf_clr = [0 1 0];
    axes(h.axes_anatomy); hold on;
    if num_patches>6 && h.slider_transparency_lf_grids.Value==0
        msgbox(sprintf('Head Model is already plotted but LF transparency is ZERO\n\nChange HDM transparency to see Head Model'));
    elseif num_patches>6
        msgbox(sprintf('Head Model is already plotted\n\nIf you cannot see it then try "plot 3D" followed by "plot HeadModel" again'));
    else
        h.lf_grids_patch = scatter3(h.anatomy.leadfield.voxel_pos(:,1),h.anatomy.leadfield.voxel_pos(:,2),h.anatomy.leadfield.voxel_pos(:,3),'Marker','o','SizeData',2,...
            'MarkerFaceColor',lf_clr,'MarkerEdgeColor',lf_clr,'MarkerFaceAlpha',h.slider_transparency_lf_grids.Value,'MarkerEdgeAlpha',h.slider_transparency_lf_grids.Value);
        set_transparency([],[],'lf');
    end
end
