function set_transparency(varargin)
global h

src = varargin;
try
    switch src{3}
        case 'scalp'
            h.scalp_plot_patch.FaceAlpha = h.slider_transparency_scalp.Value;
        case 'brain'
            h.brain_plot_patch.FaceAlpha = h.slider_transparency_brain.Value;
        case 'topo'
            h.topo_plot_patch.FaceAlpha = h.slider_transparency_topo.Value;
        case 'hdm'
            for p=1:length(h.hdm_patch); h.hdm_patch(p).FaceAlpha = h.slider_transparency_hdm.Value; end
        case 'lf'
            %             for p=1:length(h.lf_grids_patch)
            h.lf_grids_patch.MarkerFaceAlpha = h.slider_transparency_lf_grids.Value;
            h.lf_grids_patch.MarkerEdgeAlpha = h.slider_transparency_lf_grids.Value;
            %             end
    end
catch
end

