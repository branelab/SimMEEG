function sm_menu_inv_analyses_CallBack(varargin)
% set menu options for the in sub-panel "Analyses" on the "Source Modeling" Tab
global h
h.panel_3D_PLV_PLI.Visible = 'off';

if isempty(h.listbox_inv_solns.String)
    hm=warndlg(sprintf('Please load perform inverse modeling using "Run Modeling"'),'No Inverse Solution');
    h.menu_inv_analyses.Value=1;
else
    %% %%%%% PLV & PLI analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if h.menu_inv_analyses.Value==2
        h.panel_3D_PLV_PLI.Visible = 'on';
        
        %% setting Seed indices
        seed_idx = [];
        if any(h.listbox_inv_plv_seed_locs.Value == 1)   % add True Source locations
            clear vx_idx;
            % finding nearest voel to true source within the leadfield because inv_soln may be calculated using different leadfields than the one that simulated the sens_data
            for v=1:3; vx_idx(v) = find_nearest_voxel(h.sim_data.cfg.source.vx_locs(v,:),h.inv_soln(h.current_inv_soln).leadfield.voxel_pos);  end  % source's voxel index from leadfield positions
            seed_idx = [seed_idx vx_idx];
        end
        if any(h.listbox_inv_plv_seed_locs.Value == 2)   % add "Hit" Peak Source locations
            seed_idx = [seed_idx h.inv_soln(h.current_inv_soln).peak_idx(1:3)'];
        end
        if any(h.listbox_inv_plv_seed_locs.Value == 3)   % add Peak Source locations
            seed_idx = [seed_idx h.inv_soln(h.current_inv_soln).peak_idx'];
        end
        
        if any(h.listbox_inv_plv_seed_locs.Value == 4)   % add Lead Field locations - all or downsampled
            if isfield(h.inv_soln(h.current_inv_soln),'plv_leadfield_grid_idx')
                sm_downsample_leadfield;
                lf_grid_idx = h.inv_soln(h.current_inv_soln).plv_leadfield_grid_idx;
            else
                lf_grid_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
            end
            
            seed_idx = [seed_idx lf_grid_idx];
        end
        h.inv_soln(h.current_inv_soln).plv_seed_idx = unique(seed_idx,'stable');
        
        %% setting Comparison indices
        comp_idx = [];
        if any(h.listbox_inv_plv_comp_locs.Value == 1)   % add True Source locations
            clear vx_idx;
            % finding nearest voel to true source within the leadfield because inv_soln may be calculated using different leadfields than the one that simulated the sens_data
            for v=1:3; vx_idx(v) = find_nearest_voxel(h.sim_data.cfg.source.vx_locs(v,:),h.inv_soln(h.current_inv_soln).leadfield.voxel_pos);  end  % source's voxel index from leadfield positions
            comp_idx  = [comp_idx vx_idx];
        end
        if any(h.listbox_inv_plv_comp_locs.Value == 2)   % add "Hit" Peak Source locations
            comp_idx = [comp_idx h.inv_soln(h.current_inv_soln).peak_idx(1:3)'];
        end
        if any(h.listbox_inv_plv_comp_locs.Value == 3)   % add Peak Source locations
            comp_idx = [comp_idx h.inv_soln(h.current_inv_soln).peak_idx'];
        end
        
        if any(h.listbox_inv_plv_comp_locs.Value == 4)   % add Lead Field locations - all or downsampled
            if isfield(h.inv_soln(h.current_inv_soln),'plv_leadfield_grid_idx')
                lf_grid_idx = h.inv_soln(h.current_inv_soln).plv_leadfield_grid_idx;
            else
                lf_grid_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
            end
            
            try comp_idx = [comp_idx lf_grid_idx]; catch; comp_idx = [comp_idx' lf_grid_idx]; end
            
        end
        h.inv_soln(h.current_inv_soln).plv_comp_idx = unique(comp_idx,'stable');
        
        if length(h.inv_soln(h.current_inv_soln).plv_seed_idx)>1
            if length(h.inv_soln(h.current_inv_soln).plv_seed_idx) == length(h.inv_soln(h.current_inv_soln).plv_comp_idx) ...
                    && isempty( setdiff(h.inv_soln(h.current_inv_soln).plv_seed_idx,h.inv_soln(h.current_inv_soln).plv_comp_idx) )  % same voxels in seed and comps thus nchoose2 calculation
                num_contrasts = nchoosek(length(h.inv_soln(h.current_inv_soln).plv_seed_idx) ,2);
            else
                num_contrasts = length(h.inv_soln(h.current_inv_soln).plv_seed_idx) * length(h.inv_soln(h.current_inv_soln).plv_comp_idx);
            end
        else
            num_contrasts = [];
        end
        
        h.edit_inv_plv_num_seeds_comps_txt.String =sprintf('Seeds = %.f\nComparisons = %.f\nFC Contrasts = %.f',...
            length(h.inv_soln(h.current_inv_soln).plv_seed_idx),length(h.inv_soln(h.current_inv_soln).plv_comp_idx),...
            num_contrasts );
        sm_plot_plv_locs;
        
        %% %%%%% OTHER ANALYSES can go under here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif h.menu_inv_analyses.Value==3
        
        
    end
    
end




