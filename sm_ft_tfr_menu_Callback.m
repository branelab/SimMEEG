function sm_ft_tfr_menu_Callback(varargin)
global h

%% turning options 'off'
    h.edit_preproc_ft_tfr_t_ftimwin_txt.Visible = 'off'; h.edit_preproc_ft_tfr_t_ftimwin.Visible = 'off';
    h.edit_preproc_ft_tfr_toi_txt.Visible = 'off'; h.edit_preproc_ft_tfr_toi.Visible = 'off';
    h.edit_preproc_ft_tfr_t_ftimwin_txt.Visible = 'off'; h.edit_preproc_ft_tfr_t_ftimwin.Visible = 'off';
    h.edit_preproc_ft_tfr_toi_txt.Visible = 'off'; h.edit_preproc_ft_tfr_toi.Visible = 'off';
    h.edit_preproc_ft_tfr_wt_width_txt.Visible = 'off'; h.edit_preproc_ft_tfr_wt_width.Visible = 'off';
    h.edit_preproc_ft_tfr_gwidth_txt.Visible = 'off'; h.edit_preproc_ft_tfr_gwidth.Visible = 'off';

% if length(str2num(h.edit_preproc_ft_tfr_foilim.String))==2
%     h.edit_preproc_ft_tfr_tapsmofrq_txt.Visible = 'off'; h.edit_preproc_ft_tfr_tapsmofrq.Visible = 'off';
% else
    h.edit_preproc_ft_tfr_tapsmofrq_txt.Visible = 'on'; h.edit_preproc_ft_tfr_tapsmofrq.Visible = 'on';
% end

%% selectively turning options 'on'
switch h.menu_preproc_ft_tfr_method.String{h.menu_preproc_ft_tfr_method.Value}
    
    case {'mtmfft'}
    case {'mtmconvol' }
        h.edit_preproc_ft_tfr_t_ftimwin_txt.Visible = 'on'; h.edit_preproc_ft_tfr_t_ftimwin.Visible = 'on';
        h.edit_preproc_ft_tfr_toi_txt.Visible = 'on'; h.edit_preproc_ft_tfr_toi.Visible = 'on';
    case {'wavelet'}
        h.edit_preproc_ft_tfr_t_ftimwin_txt.Visible = 'on'; h.edit_preproc_ft_tfr_t_ftimwin.Visible = 'on';
        h.edit_preproc_ft_tfr_toi_txt.Visible = 'on'; h.edit_preproc_ft_tfr_toi.Visible = 'on';
        h.edit_preproc_ft_tfr_wt_width_txt.Visible = 'on'; h.edit_preproc_ft_tfr_wt_width.Visible = 'on';
        h.edit_preproc_ft_tfr_gwidth_txt.Visible = 'on'; h.edit_preproc_ft_tfr_gwidth.Visible = 'on';
    case {'tfr'}
         h.edit_preproc_ft_tfr_wt_width_txt.Visible = 'on'; h.edit_preproc_ft_tfr_wt_width.Visible = 'on';
        h.edit_preproc_ft_tfr_gwidth_txt.Visible = 'on'; h.edit_preproc_ft_tfr_gwidth.Visible = 'on';
         h.edit_preproc_ft_tfr_toi_txt.Visible = 'on'; h.edit_preproc_ft_tfr_toi.Visible = 'on';
  case {'mvar'}
end



