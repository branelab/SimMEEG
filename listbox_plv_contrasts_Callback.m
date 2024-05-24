function listbox_plv_contrasts_Callback(varargin)
global h

%
% fc_idx = h.listbox_plv_contrasts.Value;

% h.current_inv_TFR_waves_all

if h.radio_plot_TFR_waves.Value==1 && (h.radio_inv_plot_peak_tfr_connectivity.Value==1 || h.radio_inv_plot_true_tfr_connectivity.Value==1)
    %% plotting connectivity
    sm_plot_connectivity_graph;
    %     sm_plot_true_tfr_connectivity;
    x.Type = 'none';
    sm_get_inv_tfr_point(x);
    sm_plot_connectivity(varargin{end});
    
    %     switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
    %         case {'PLV' 'PLI' 'dPLI'}
    %             hit_plv_idx = [1 2; 1 3; 2 3];
    %             v_idx=[];
    %             for v=h.listbox_plv_contrasts.Value
    %                 v_idx = [v_idx, find(sum(hit_plv_idx==v,2)==1)'];
    %             end
    %             v_idx = unique(v_idx);
    %             v_idx = v_idx(isgraphics(h.current_inv_TFR_waves(v_idx),'line'));
    %             for v=v_idx
    %                 h.current_inv_TFR_waves(v).Color(4) = 1;
    %             end
    %         case {'Total Power' 'Evoked Power' 'Induced Power'}
    %             for v=h.listbox_plv_contrasts.Value
    %                 h.current_inv_TFR_waves(v).Color(4) = 1;
    %             end
    %     end
else
    
    if h.radio_inv_plot_peak_tfr_connectivity.Value == 1   %&& h.menu_inv_tfr_type.Value>3  % also update TFR plot if its on
        sm_plot_tfr_connectivity; x.Type = 'none'; sm_get_inv_tfr_point(x);
        sm_plot_connectivity_graph;
    end
    if h.radio_inv_plot_true_tfr_connectivity.Value == 1 %&& h.menu_inv_tfr_type.Value>3  % also update TFR plot if its on
        sm_plot_true_tfr_connectivity; x.Type = 'none';
        sm_get_inv_tfr_point(x);
        sm_plot_connectivity_graph;
    end
    
end

%% plotting highlight circle at TFR/FC locations
if ~isempty(h.current_3D_plv_contrasts)
    switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
        case {'Total Power' 'Evoked Power' 'Induced Power'}
            x_locs = h.current_3D_plv_contrasts( h.current_3D_plv_contrasts_listbox_order(h.listbox_plv_contrasts.Value));
        case {'PLV' 'PLI' 'dPLI'}
            x_locs = h.current_3D_plv_contrasts( h.current_3D_plv_contrasts_listbox_order(h.listbox_plv_contrasts.Value),: );
            if size(x_locs,1)==1
                x_locs = x_locs(~isnan(x_locs));
            else
                x_locs = x_locs(sum(~isnan(x_locs),2)==2);
            end
    end
    
    
    if all(~isnan(x_locs))
        sel_locs = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(x_locs ,:);
    else
        return
    end
else
    sel_locs = [];
end

try
    if isfield(h,'current_peak_selected')
        if isvalid(h.current_peak_selected)
            h.current_peak_selected.XData = sel_locs(:,1); h.current_peak_selected.YData = sel_locs(:,2); h.current_peak_selected.ZData = sel_locs(:,3);
            %         h.current_inv_swf_plots(h.listbox_plv_contrasts.Value).Color(1:3) = [1 0 1];
            %         for v=h.listbox_plv_contrasts.Value
            %             h.current_inv_swf_plots(v).Color(4) = 1;
            %             h.current_inv_swf_plots_fft(v).Color(4) = 1;
            %         end
        else
            h.current_peak_selected = scatter3(h.axes_3D_images,sel_locs(:,1),sel_locs(:,2),sel_locs(:,3),'o','MarkerEdgeColor',[1 0 1]*0,'SizeData',h.map3D_peak_locs(1).SizeData,'linewidth',2);
        end
    else
        h.current_peak_selected = scatter3(h.axes_3D_images,sel_locs(:,1),sel_locs(:,2),sel_locs(:,3),'o','MarkerEdgeColor',[1 0 1]*0,'SizeData',h.map3D_peak_locs(1).SizeData,'linewidth',2);
    end
catch
end