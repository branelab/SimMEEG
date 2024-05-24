function menu_inv_soln_txt_Callback(varargin)
global h
h.panel_SPA_params.Visible = 'off';
h.panel_inv_ft_params.Visible = 'off';
h.panel_inv_bst_params.Visible = 'off';
h.menu_inv_weightnorm_txt.Visible = 'off'; h.menu_inv_weightnorm.Visible='off';
h.menu_inv_scalesourcecov_txt.Visible = 'off'; h.menu_inv_scalesourcecov.Visible = 'off';
h.menu_inv_maxvectorori_txt.Visible = 'off'; h.menu_inv_maxvectorori.Visible = 'off';
h.menu_inv_datatype_txt.Visible = 'off'; h.menu_inv_datatype.Visible = 'off'; h.menu_inv_datatype.Enable = 'on';

%% BST Panel Setting - making specifics options visible 
h.edit_inv_bst_reg_lambda_snr_txt.Visible = 'on'; h.edit_inv_bst_reg_lambda_snr.Visible = 'on';
h.edit_inv_bst_reg_noise_factor_txt.Visible = 'on'; h.edit_inv_bst_reg_noise_factor.Visible = 'on';
h.radio_inv_bst_depth_weight.Visible = 'on';
h.edit_inv_bst_depth_order_txt.Visible = 'on'; h.edit_inv_bst_depth_order.Visible = 'on';
h.edit_inv_bst_depth_max_txt.Visible = 'on'; h.edit_inv_bst_depth_max.Visible = 'on';
h.menu_inv_bst_maxvectorori_txt.Visible = 'on'; h.menu_inv_bst_maxvectorori.Visible = 'on';
if h.menu_inv_bst_reg_noise_type.Value >1
    h.edit_inv_bst_reg_noise_factor_txt.Visible = 'off'; h.edit_inv_bst_reg_noise_factor.Visible = 'off'; 
end
if h.radio_inv_bst_depth_weight.Value == 0
    h.edit_inv_bst_depth_order_txt.Visible = 'off'; h.edit_inv_bst_depth_order.Visible = 'off'; 
    h.edit_inv_bst_depth_max_txt.Visible = 'off'; h.edit_inv_bst_depth_max.Visible = 'off'; 
end


%% selecting InvSoln
switch h.menu_inv_soln.String{h.menu_inv_soln.Value}
    case {'SPA'}    % single-source scalar beamformer seee Herdman et al., 2018 Brain Topography
        h.panel_SPA_params.Visible = 'on';
        h.panel_SPA_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_SPA_loc_flag.String = {'pseudo Z','event-related','reduced event-related','activity index'};
        h.edit_SPA_max_sources_txt.Visible = 'off'; h.edit_SPA_max_sources.Visible = 'off';
        
    case {'SIA' 'MIA'}      % multi-source scalar beamformers seee Herdman et al., 2018 Brain Topography
        h.panel_SPA_params.Visible = 'on';
        h.panel_SPA_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_SPA_loc_flag.String = {'MPZ','MER','RMER','MAI'};
        h.edit_SPA_max_sources_txt.Visible = 'on'; h.edit_SPA_max_sources.Visible = 'on';
        
    case {'LCMV (BST)'}
        h.panel_inv_bst_params.Visible = 'on';
        h.edit_inv_bst_reg_lambda_snr_txt.Visible = 'off'; h.edit_inv_bst_reg_lambda_snr.Visible = 'off';
        h.radio_inv_bst_depth_weight.Visible = 'off';
        h.edit_inv_bst_depth_order_txt.Visible = 'off'; h.edit_inv_bst_depth_order.Visible = 'off';
        h.edit_inv_bst_depth_max_txt.Visible = 'off'; h.edit_inv_bst_depth_max.Visible = 'off';

    case {'MNE (BST)'}
    h.panel_inv_bst_params.Visible = 'on';
    case {'sLORETA (BST)'}
    h.panel_inv_bst_params.Visible = 'on';
        h.radio_inv_bst_depth_weight.Visible = 'off';
        h.edit_inv_bst_depth_order_txt.Visible = 'off'; h.edit_inv_bst_depth_order.Visible = 'off';
        h.edit_inv_bst_depth_max_txt.Visible = 'off'; h.edit_inv_bst_depth_max.Visible = 'off';
    case {'LCMV (FT)'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_inv_weightnorm_txt.Visible = 'on'; h.menu_inv_weightnorm.Visible='on';
    case {'SAM (FT)'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'off';
%         h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
%         h.menu_inv_weightnorm_txt.Visible = 'on'; h.menu_inv_weightnorm.Visible='on';
        
    case {'MNE (FT)'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
         h.menu_inv_scalesourcecov_txt.Visible = 'on'; h.menu_inv_scalesourcecov.Visible = 'on';
       h.menu_inv_maxvectorori_txt.Visible = 'on'; h.menu_inv_maxvectorori.Visible = 'on';
    
    case {'eLORETA (FT)'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_inv_maxvectorori_txt.Visible = 'on'; h.menu_inv_maxvectorori.Visible = 'on';

    case {'sLORETA (FT)'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        
    case {'dics (FT)'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_inv_datatype_txt.Visible = 'on'; h.menu_inv_datatype.Visible = 'on';
        h.menu_inv_datatype.Value = 2; h.menu_inv_datatype.Enable = 'inactive';
           
    case {'pcc (FT)'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_inv_datatype_txt.Visible = 'on'; h.menu_inv_datatype.Visible = 'on';
         
    case {'sMCMV' 'bRAPBeam' 'TrapMUSIC'}   % Alex Moiseev's sub-space multi-source scalar beamformer and TRAP-MUSIC
        h.panel_SPA_params.Visible = 'on';
        h.panel_SPA_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_SPA_loc_flag.String = {'MPZ','MER','RMER','MAI'};
        h.edit_SPA_max_sources_txt.Visible = 'on'; h.edit_SPA_max_sources.Visible = 'on';
end


