function get_patch_axes_anatomy(varargin)
global h

xidx = findobj(h.axes_anatomy.Children,'Type','Patch');
for t=1:length(xidx)
    if strcmpi(xidx(t).FaceColor,'interp') % topo plot patch
        h.topo_plot_patch = xidx(t);
    elseif sum(xidx(t).FaceColor == h.sens_clr)==3 || sum(xidx(t).FaceColor == h.chan_clr)==3   % sensor plot patch
        h.sens_plot_patch = xidx(t);
    end
end
