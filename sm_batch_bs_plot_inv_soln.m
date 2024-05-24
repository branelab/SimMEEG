function h = sm_batch_bs_plot_inv_soln(h)

%% Resetting current_peak...
h.current_3D_peak_idx = [];
h.current_3D_peak_voxels = [];
h.current_peak_hit_idx = [];
h.current_peak_miss_idx = [];
h.current_peak_fa_idx = [];
h.current_3D_thresh =[];
h.current_inv_soln_show_peak_idx =[];
h.current_inv_soln_hide_peak_idx =[];
h.current_norm_peak_swf =[]; h.current_peak_swf = [];

%%
%% 
sel_idx = h.listbox_peaks_found.Value;

%% map type: 'image' 'norm image' 'point spread function' ...
% h = sm_plot_replace_3Dmap(h);


%% disable slider so multiple plot attempts are not done
% h.slider_3D_image_thresh.Enable = 'inactive';
%% intializing plot settings
seed_idx = 1:length(h.cfg.source.vx_idx);  ln_wdth = 1; ln_wdth2 = 1;

min_max = h.inv_soln(h.current_inv_soln).soln.plot_min_max; %[min(h.inv_soln(h.current_inv_soln).soln.P.img) max(h.inv_soln(h.current_inv_soln).soln.P.img)];
if min_max(1)==min_max(2); min_max(2)=min_max(1)+1;
elseif min_max(1)>min_max(2); min_max(2)=min_max(1)+1;
end

null_thresh = 0; %h.inv_soln(h.current_inv_soln).soln.plot_thresh;
search_dist = str2num(h.edit_inv_peak_spread.String); % minimum distance between peak spreads for sptaially searching for peaks within the image
%% Getting peak_voxels for image in inv_soln.soln.P.img
switch h.inv_soln(h.current_inv_soln).headmodel_type
    case 'Whole Brain' % whole brain

        if h.radio_find_spatiotemp_peaks.Value == 1 % Find peaks using spatiotemporal search across active interval
            
            msgbox('Spatiotemporal mapping not implemented yet for command-line batch scripting','Inverse Modeling'); 
            return
            sm_spatiotemp_mapping;
            %             h.inv_soln(h.current_inv_soln).soln.plot_min_max = [min(h.inv_soln(h.current_inv_soln).soln.P.img) max(h.inv_soln(h.current_inv_soln).soln.P.img)];
            h.inv_soln(h.current_inv_soln).soln.plot_min_max = [-max(abs(h.inv_soln(h.current_inv_soln).soln.P.img)) max(abs(h.inv_soln(h.current_inv_soln).soln.P.img))];
            min_max = h.inv_soln(h.current_inv_soln).soln.plot_min_max;
            null_thresh = 0; %h.slider_3D_image_thresh.Value;
            %% Plot Whole Brain Soln combined spatiotemporal maps --> current_3D_peak_voxels set in sm_spatiotemp_mapping
            [~,~,h.s1,h.p1,h.func3D,h.anat3D]=bl_plot_lcmv_peak_img_FT_new(h.inv_soln(h.current_inv_soln).soln.P.img,...
                h.inv_soln(h.current_inv_soln).soln.ori,null_thresh,search_dist,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,parula(255),...
                min_max,vol,h.anatomy.sens,h.cfg.study.bl_bmf.vw_angle,pk_flag,1,h.inv_soln(h.current_inv_soln).leadfield.pos,inside_idx,h.cfg.study.bl_bmf.slice_orient,sFaceAlpha);
        else        % plot using map created by inverse solution
            
            %% Plot Whole Brain Soln
            if isfield(h.inv_soln(h.current_inv_soln).soln.P,'img_org')     % reverting back to original inv_soln image
                switch h.inv_soln(h.current_inv_soln).Type
                    case {'SIA' 'MIA'}
                        h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).soln.P.img_org;     % using standard MCMV
                    case {'SPA' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}
                        h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).soln.P.img_org;
                    case {'LCMV (FT)' 'sLORETA (FT)' 'MNE (FT)' 'dics (FT)' 'pcc  (FT)' 'SAM (FT)' 'LCMV (BST)' 'sLORETA (BST)' 'MNE (BST)'}
                        inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
                        h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).soln.P.img_org(inside_idx);
                end
            else
                switch h.inv_soln(h.current_inv_soln).Type
                    case {'SIA' 'MIA'}
                    case {'LCMV (FT)' 'sLORETA (FT)' 'MNE (FT)' 'dics (FT)' 'pcc  (FT)' 'SAM (FT)' 'LCMV (BST)' 'sLORETA (BST)' 'MNE (BST)'}
                        if size(h.inv_soln(h.current_inv_soln).soln.P.img,1)>size(h.inv_soln(h.current_inv_soln).leadfield.H,3)
                            inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
                            h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).soln.P.img(inside_idx);
                        end
                end
            end
%             [h.current_3D_peak_voxels,h.current_3D_peak_idx,h.s1,h.p1,h.func3D,h.anat3D]=bl_plot_lcmv_peak_img_FT_new(h.inv_soln(h.current_inv_soln).soln.P.img,...
%                 h.inv_soln(h.current_inv_soln).soln.ori,null_thresh,search_dist,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,parula(255),...
%                 min_max,vol,h.anatomy.sens,h.cfg.study.bl_bmf.vw_angle,pk_flag,1,h.inv_soln(h.current_inv_soln).leadfield.pos,inside_idx,h.cfg.study.bl_bmf.slice_orient,sFaceAlpha);
%                         
            voxel_vals=[h.inv_soln(h.current_inv_soln).leadfield.voxel_pos, h.inv_soln(h.current_inv_soln).soln.P.img];
            thresh_val = null_thresh;
            [h.current_3D_peak_voxels,h.current_3D_peak_idx]=BRANELab_find_peak_voxel_thresh(voxel_vals,thresh_val,search_dist);
 
        end
       
        
    case 'Cortical Surface'    % Cortical Surface
        
        if h.radio_find_spatiotemp_peaks.Value == 1 % Find peaks using spatiotemporal search across active interval
                       msgbox('Spatiotemporal mapping not implemented yet for command-line batch scripting','Inverse Modeling'); 
            return
            sm_spatiotemp_mapping;
        else        % plot using map created by inverse solution
            % Peak locations
            voxel_vals=[h.inv_soln(h.current_inv_soln).leadfield.voxel_pos, vol(3).img];
            thresh_val = null_thresh;
            [h.current_3D_peak_voxels,h.current_3D_peak_idx]=BRANELab_find_peak_voxel_thresh(voxel_vals,thresh_val,search_dist);
        end
         
end

if ~isempty(h.current_3D_peak_idx)
    h.current_3D_peak_voxels(:,5) = h.current_3D_peak_idx;


%% Do not update these once called set here - only to be set during run_source_modeling.m
if h.run_inv_soln_flag == 1  % Update only on first instance of run_source_modeling.m
    h.inv_soln(h.current_inv_soln).peak_voxels = h.current_3D_peak_voxels;  % setting these to be re-ordered in [true_idx fa_idx]
    h.inv_soln(h.current_inv_soln).peak_idx = h.current_3D_peak_idx;
end
%% setting slider value to zero
h.slider_3D_image_thresh.Value = 0;


%% Plotting Scatter Locations and Orientations after searching for Hits & Reordering peaks to be in order of those nearest to source 1, 2, 3 locations
if isempty(h.current_3D_peak_voxels) % no peaks found
    h.current_3D_peak_idx =[]; h.current_3D_peak_voxels = []; h.current_norm_peak_swf =[]; h.current_peak_swf = [];
    sm_calc_localizer_performance;
else
    
    if h.radio_find_spatiotemp_peaks.Value == 1
        % Do Nothing because nearest peaks already found
    else
        if h.run_inv_soln_flag ==1  % called during by run_source_modeling.m
            h.current_inv_soln_show_peak_idx = 1:length(h.current_3D_peak_idx);
            h = sm_batch_sm_search_for_hits(h,'initial search'); % reorders "h.current_3D_peak_voxels" after performing search based on user's selection of 'Nearest', 'Wave Error', 'Wave Correlation'
        else
            h.slider_3D_image_thresh.Value = 0; %h.inv_soln(h.current_inv_soln).soln.plot_thresh;
            h = sm_batch_sm_search_for_hits(h,'slider thresh');
        end
    end
    h = sm_batch_sm_calc_localizer_performance(h);
    
    %% peak locs
    if isempty(h.current_3D_peak_voxels) % no peaks found
        text(0,0,0,'No Peak Sources Found'); h.listbox_peaks_found.String = ''; h.listbox_peaks_found.Value = 1;
        view(h.axes_3D_images,cur_vw); axis(h.axes_3D_images,'tight');
        h.current_3D_peak_idx =[]; h.current_3D_peak_voxels = []; h.current_norm_peak_swf =[]; h.current_peak_swf = [];
        if isfield(h,'colorbar_3D')
            if isvalid(h.colorbar_3D)
                delete(h.colorbar_3D); h.colorbar_3D = colorbar(h.axes_3D_images,'Location','southoutside','Position',[.75 .35 .2 .03]);
            else
                h.colorbar_3D = colorbar(h.axes_3D_images,'Location','southoutside','Position',[.75 .35 .2 .03]);
            end
        else
            h.colorbar_3D = colorbar(h.axes_3D_images,'Location','southoutside','Position',[.75 .35 .2 .03]);
        end
    else
        
        h = sm_batch_bs_calc_errors_inv_soln(h);
        h = sm_batch_sm_calc_localizer_performance(h);
    end
end

h.slider_3D_image_thresh.Value = h.inv_soln(h.current_inv_soln).soln.plot_min_max(1);
if h.slider_3D_image_thresh.Value<h.slider_3D_image_thresh.Min
    h.slider_3D_image_thresh.Value = h.slider_3D_image_thresh.Min;
elseif h.slider_3D_image_thresh.Value>h.slider_3D_image_thresh.Max
    h.slider_3D_image_thresh.Value = h.slider_3D_image_thresh.Max;
end

if h.run_inv_soln_flag ==1  % called during by run_source_modeling.m
    h.slider_3D_image_thresh.Value = 0;
else
    h.slider_3D_image_thresh.Value = h.inv_soln(h.current_inv_soln).soln.plot_thresh;
end

%% updating Localization Performance Results Graphs
if ~isfield(h.inv_soln(h.current_inv_soln),'classifier_results')
     h = sm_batch_sm_calc_results_loc_performance(h);
else
    if isempty(h.inv_soln(h.current_inv_soln).classifier_results)
         h = sm_batch_sm_calc_results_loc_performance(h);
    end
end
h = sm_batch_sm_calc_results_loc_performance(h);

else
    h.current_3D_peak_voxels = [];
    h.inv_soln(h.current_inv_soln).classifier_results = []; 
    h.inv_soln(h.current_inv_soln).classifier_metrics = []; 
    
    if ~isfield(h.inv_soln(h.current_inv_soln),'classifier_results')
        h = sm_batch_sm_calc_results_loc_performance(h);
    else
        if isempty(h.inv_soln(h.current_inv_soln).classifier_results)
            h = sm_batch_sm_calc_results_loc_performance(h);
        end
    end

end
