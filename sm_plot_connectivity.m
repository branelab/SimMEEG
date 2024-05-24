function sm_plot_connectivity(varargin)
global h
h.FC_alpha_level = str2num(h.edit_3D_FC_stats_alpha_level.String);
switch varargin{end}
    case 'Peaks'
        if h.radio_inv_plot_peak_tfr_connectivity.Value==1
            sm_plot_tfr_connectivity; %set_inv_tfr_caxis('FC');
        end
    case 'True'
        if h.radio_inv_plot_true_tfr_connectivity.Value==1
            sm_plot_true_tfr_connectivity; %set_inv_tfr_caxis('FC');
        end
        x.Type = 'none'; sm_get_inv_tfr_point(x);
    case 'Both'
        if h.radio_inv_plot_peak_tfr_connectivity.Value==1
            sm_plot_tfr_connectivity; %set_inv_tfr_caxis('FC');
        end
        if h.radio_inv_plot_true_tfr_connectivity.Value==1
            sm_plot_true_tfr_connectivity; %set_inv_tfr_caxis('FC');
        end
end





