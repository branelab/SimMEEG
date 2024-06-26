function set_inv_tfr_caxis(varargin)
global h

% Peak Source TFR plot
switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
    case {'Total Power' 'Evoked Power' 'Induced Power'}
        c_lim = str2num(h.edit_inv_tfr_caxis.String);
        if c_lim(1)>=c_lim(2); c_lim(1) = -c_lim(2); h.edit_inv_tfr_caxis.String = num2str(c_lim); end
        if unique(c_lim)==0; c_lim = [-1 1]; h.edit_inv_tfr_caxis.String = num2str(c_lim); end
    case {'PLV' 'PLI' 'dPLI'}
        c_lim = str2num(h.edit_inv_plv_caxis.String);
        if c_lim(1)>=c_lim(2); c_lim(1) = -c_lim(2); h.edit_inv_plv_caxis.String = num2str(c_lim); end
        if unique(c_lim)==0; c_lim = [-1 1]; h.edit_inv_plv_caxis.String = num2str(c_lim); end
end

if h.radio_inv_plot_peak_tfr_connectivity.Value == 1
    h.axes_inv_soln_tfr.CLim = c_lim;
end

% True source TFR plot
if h.radio_inv_plot_true_tfr_connectivity.Value == 1
    h.axes_true_source_tfr.CLim = c_lim;
end

if h.radio_plot_TFR_waves.Value == 1  % update wave Ylim
    h.axes_source_waves.YLim = str2num(h.edit_inv_tfr_caxis.String);
end

%% updating connectivity graph 
if strmatch(varargin{end},'FC')
sm_plot_connectivity_graph;
end

