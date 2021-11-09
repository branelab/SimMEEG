function listbox_peaks_found_Callback(varargin)
global h

if h.radio_3D_peak_locs.Value==1    % peaks need to be shown on map to highlight them
    if length(h.listbox_inv_solns.Value)==1
        ln_clr = h.map3D_peak_locs(1).CData; %bsxfun(@mtimes,ones(length(h.current_3D_peak_idx),3),[.9 .6 .3]*.75);  %lines(length(h4.current_3D_peak_idx));
        %         ln_clr(1:3,:) = h.src_clr;
        for v=1:size(h.current_peak_swf,2)
            h.current_inv_swf_plots(v).Color = ln_clr(v,:);
            h.current_inv_swf_plots_fft(v).Color = ln_clr(v,:);
            %     if v>3
            h.current_inv_swf_plots(v).Color(4) = h.false_positive_lineAlpha;
            h.current_inv_swf_plots_fft(v).Color(4) = h.false_positive_lineAlpha;
            %     end
        end
        
        sel_locs = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx(h.listbox_peaks_found.Value),:);
        if isfield(h,'current_peak_selected')
            if isvalid(h.current_peak_selected)
                h.current_peak_selected.XData = sel_locs(:,1); h.current_peak_selected.YData = sel_locs(:,2); h.current_peak_selected.ZData = sel_locs(:,3);
                %         h.current_inv_swf_plots(h.listbox_peaks_found.Value).Color(1:3) = [1 0 1];
                for v=1:length(h.listbox_peaks_found.Value)
                    h.current_inv_swf_plots(h.listbox_peaks_found.Value(v)).Color(4) = 1;
                    h.current_inv_swf_plots_fft(h.listbox_peaks_found.Value(v)).Color(4) = 1;
                end
            else
                h.current_peak_selected = scatter3(h.axes_3D_images,sel_locs(:,1),sel_locs(:,2),sel_locs(:,3),'o','MarkerEdgeColor',[1 0 1]*0,'SizeData',h.map3D_peak_locs(1).SizeData,'linewidth',2);
            end
        else
            h.current_peak_selected = scatter3(h.axes_3D_images,sel_locs(:,1),sel_locs(:,2),sel_locs(:,3),'o','MarkerEdgeColor',[1 0 1]*0,'SizeData',h.map3D_peak_locs(1).SizeData,'linewidth',2);
        end
    end
end
if h.radio_inv_plot_peak_tfr_connectivity.Value == 1  %&& h.menu_inv_tfr_type.Value>3  % also update TFR plot if its on
    sm_plot_tfr_connectivity; x.Type = 'none'; sm_get_inv_tfr_point(x);
    sm_plot_connectivity_graph;
end

if h.radio_inv_plot_true_tfr_connectivity.Value == 1  %&& h.menu_inv_tfr_type.Value>3  % also update TFR plot if its on
    sm_plot_true_tfr_connectivity; x.Type = 'none'; sm_get_inv_tfr_point(x);
    sm_plot_connectivity_graph;
end

h.axes_source_waves.XLim = str2num(h.edit_plot_time_int.String);