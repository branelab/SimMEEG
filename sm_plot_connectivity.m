function sm_plot_connectivity(varargin)
global h

h.FC_alpha_level = str2num(h.edit_3D_FC_stats_alpha_level.String);
if h.radio_inv_plot_peak_tfr_connectivity.Value==1
sm_plot_tfr_connectivity;
end
if h.radio_inv_plot_true_tfr_connectivity.Value==1
sm_plot_true_tfr_connectivity;
end
x.Type = 'none'; sm_get_inv_tfr_point(x);



