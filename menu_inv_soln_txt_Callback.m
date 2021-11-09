function menu_inv_soln_txt_Callback(varargin)
global h
h.panel_SPA_params.Visible = 'off';
h.panel_inv_ft_params.Visible = 'off';
h.menu_inv_weightnorm_txt.Visible = 'off'; h.menu_inv_weightnorm.Visible='off';
h.menu_inv_scalesourcecov_txt.Visible = 'off'; h.menu_inv_scalesourcecov.Visible = 'off';
h.menu_inv_maxvectorori_txt.Visible = 'off'; h.menu_inv_maxvectorori.Visible = 'off';
h.menu_inv_datatype_txt.Visible = 'off'; h.menu_inv_datatype.Visible = 'off'; h.menu_inv_datatype.Enable = 'on';


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
        
    case {'LCMV'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_inv_weightnorm_txt.Visible = 'on'; h.menu_inv_weightnorm.Visible='on';
        
    case {'MNE'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
         h.menu_inv_scalesourcecov_txt.Visible = 'on'; h.menu_inv_scalesourcecov.Visible = 'on';
       h.menu_inv_maxvectorori_txt.Visible = 'on'; h.menu_inv_maxvectorori.Visible = 'on';
    
    case {'eLORETA'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_inv_maxvectorori_txt.Visible = 'on'; h.menu_inv_maxvectorori.Visible = 'on';

    case {'sLORETA'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        
    case {'dics'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_inv_datatype_txt.Visible = 'on'; h.menu_inv_datatype.Visible = 'on';
        h.menu_inv_datatype.Value = 2; h.menu_inv_datatype.Enable = 'inactive';
           
    case {'pcc'} % Field Trip's Inverse Solutions
        h.panel_inv_ft_params.Visible = 'on';
        h.panel_inv_ft_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_inv_datatype_txt.Visible = 'on'; h.menu_inv_datatype.Visible = 'on';
         
    case {'sMCMV' 'bRAPBeam' 'TrapMUSIC'}   % Alex Moiseev's sub-space multi-source scalar beamformer and TRAP-MUSIC
        h.panel_SPA_params.Visible = 'on';
        h.panel_SPA_params.Title = sprintf('%s Parameters', h.menu_inv_soln.String{h.menu_inv_soln.Value});
        h.menu_SPA_loc_flag.String = {'MPZ','MER','RMER','MAI'};
        h.edit_SPA_max_sources_txt.Visible = 'on'; h.edit_SPA_max_sources.Visible = 'on';
end


