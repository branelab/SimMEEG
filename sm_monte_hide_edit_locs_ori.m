function sm_monte_hide_edit_locs_ori(varargin)

global h



inact_clr = [1 1 1]*.9; 
on_clr = [1 1 1];

if h.radio_monte_rand_source_locs.Value == 1 % hide edit boxes for Location Change Range and Orientation Range
    for n = 1:length(h.edit_monte_source_loc_range_X)
        h.edit_monte_source_loc_range_X(n).Enable = 'inactive'; h.edit_monte_source_loc_range_X(n).BackgroundColor = inact_clr; 
        h.edit_monte_source_loc_range_Y(n).Enable = 'inactive'; h.edit_monte_source_loc_range_Y(n).BackgroundColor = inact_clr;
        h.edit_monte_source_loc_range_Z(n).Enable = 'inactive'; h.edit_monte_source_loc_range_Z(n).BackgroundColor = inact_clr;
        h.edit_monte_source_ori_range_Az(n).Enable = 'inactive'; h.edit_monte_source_ori_range_Az(n).BackgroundColor = inact_clr;
        h.edit_monte_source_ori_range_El(n).Enable = 'inactive'; h.edit_monte_source_ori_range_El(n).BackgroundColor = inact_clr;
    end
        h.edit_monte_source_loc_StdDev_X.Enable = 'inactive'; h.edit_monte_source_loc_StdDev_X.BackgroundColor = inact_clr; 
        h.edit_monte_source_loc_StdDev_Y.Enable = 'inactive'; h.edit_monte_source_loc_StdDev_Y.BackgroundColor = inact_clr; 
        h.edit_monte_source_loc_StdDev_Z.Enable = 'inactive'; h.edit_monte_source_loc_StdDev_Z.BackgroundColor = inact_clr; 
        h.edit_monte_source_ori_StdDev_Az.Enable = 'inactive'; h.edit_monte_source_ori_StdDev_Az.BackgroundColor = inact_clr; 
        h.edit_monte_source_ori_StdDev_El.Enable = 'inactive'; h.edit_monte_source_ori_StdDev_El.BackgroundColor = inact_clr; 
else                % show edit boxes for Location Change Range and Orientation Range
     for n = 1:length(h.edit_monte_source_loc_range_X)
        h.edit_monte_source_loc_range_X(n).Enable = 'on'; h.edit_monte_source_loc_range_X(n).BackgroundColor = on_clr; 
        h.edit_monte_source_loc_range_Y(n).Enable = 'on'; h.edit_monte_source_loc_range_Y(n).BackgroundColor = on_clr;
        h.edit_monte_source_loc_range_Z(n).Enable = 'on'; h.edit_monte_source_loc_range_Z(n).BackgroundColor = on_clr;
        h.edit_monte_source_ori_range_Az(n).Enable = 'on'; h.edit_monte_source_ori_range_Az(n).BackgroundColor = on_clr;
        h.edit_monte_source_ori_range_El(n).Enable = 'on'; h.edit_monte_source_ori_range_El(n).BackgroundColor = on_clr;
    end
        h.edit_monte_source_loc_StdDev_X.Enable = 'on'; h.edit_monte_source_loc_StdDev_X.BackgroundColor = on_clr; 
        h.edit_monte_source_loc_StdDev_Y.Enable = 'on'; h.edit_monte_source_loc_StdDev_Y.BackgroundColor = on_clr; 
        h.edit_monte_source_loc_StdDev_Z.Enable = 'on'; h.edit_monte_source_loc_StdDev_Z.BackgroundColor = on_clr; 
        h.edit_monte_source_ori_StdDev_Az.Enable = 'on'; h.edit_monte_source_ori_StdDev_Az.BackgroundColor = on_clr; 
        h.edit_monte_source_ori_StdDev_El.Enable = 'on'; h.edit_monte_source_ori_StdDev_El.BackgroundColor = on_clr; 
end
