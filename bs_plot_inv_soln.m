function bs_plot_inv_soln(varargin)
global h

h.panel_3D_image_plot_msg.Visible = 'on'; drawnow;
% if h.current_inv_soln>size(h.inv_soln,2); h.current_inv_soln=size(h.inv_soln,2); end
h.panel_3D_image_plot_msg_txt.String = sprintf('%s\nRegular mapping',h.inv_soln(h.current_inv_soln).Type);

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
h.func3D_midline = ''; % midline slices 

h.slider_3D_image_thresh.Value = 0;
%% 
sel_idx = h.listbox_peaks_found.Value; 
axes(h.axes_3D_images); cla; %view(h.axes_3D_images,cur_vw);

%% map type: 'image' 'norm image' 'point spread function' ...
sm_plot_replace_3Dmap();

%%
%% disable slider so multiple plot attempts are not done
% h.slider_3D_image_thresh.Enable = 'inactive';
%% intializing plot settings
cur_vw = h.axes_3D_images.View;
h.cfg.study.bl_bmf.vw_angle = cur_vw;
seed_idx = 1:length(h.cfg.source.vx_idx);  ln_wdth = 1; ln_wdth2 = 1;
if h.monte_carlo_flag == 1
    
else
    h.waitfor_panel.Visible='off';
end
if h.radio_find_spatiotemp_peaks.Value == 1
%     hm = msgbox(sprintf('Plotting Inverse Solution\n\nSpatiotemporal mapping')); %WinOnTop(hm);
h.panel_3D_image_plot_msg_txt.String = sprintf('%s\nSpatiotemporal mapping',h.inv_soln(h.current_inv_soln).Type);

else
%     hm = msgbox(sprintf('Plotting Inverse Solution\n\nRegular mapping')); %WinOnTop(hm);
end
min_max = h.inv_soln(h.current_inv_soln).soln.plot_min_max; %[min(h.inv_soln(h.current_inv_soln).soln.P.img) max(h.inv_soln(h.current_inv_soln).soln.P.img)];
if min_max(1)==min_max(2); min_max(2)=min_max(1)+1;
elseif min_max(1)>min_max(2); min_max(2)=min_max(1)+1;
end

null_thresh = 0; %h.inv_soln(h.current_inv_soln).soln.plot_thresh;
search_dist = str2num(h.edit_inv_peak_spread.String); % minimum distance between peak spreads for sptaially searching for peaks within the image
%% Getting peak_voxels for image in inv_soln.soln.P.img
switch h.inv_soln(h.current_inv_soln).headmodel_type
    case 'Whole Brain' % whole brain
        
%         axes(h.axes_3D_images); cla; %view(h.axes_3D_images,cur_vw);
        % parameters for plotting functional maps
        h.cfg.study.bl_bmf.vw_angle = h.axes_3D_images.View;

            pk_flag=4; % map slices at peaks and midline locations
        
        inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
        sFaceAlpha = h.slider_3D_transparency_func.Value;
        %         vol = h.inv_soln(h.current_inv_soln).headmodel_mesh;
        vol = h.anatomy.mesh_volumes([1 4]);
        
        
        if h.radio_find_spatiotemp_peaks.Value == 1 && h.btn_3D_plot_peak_waves.Value == 0 % Find peaks using spatiotemporal search across active interval
            sm_spatiotemp_mapping;
            %             h.inv_soln(h.current_inv_soln).soln.plot_min_max = [min(h.inv_soln(h.current_inv_soln).soln.P.img) max(h.inv_soln(h.current_inv_soln).soln.P.img)];
%             h.inv_soln(h.current_inv_soln).soln.plot_min_max = [-max(abs(h.inv_soln(h.current_inv_soln).soln.P.img)) max(abs(h.inv_soln(h.current_inv_soln).soln.P.img))];
            min_max = h.inv_soln(h.current_inv_soln).soln.plot_min_max;
            null_thresh = 0; %h.slider_3D_image_thresh.Value;
            %% Plot Whole Brain Soln combined spatiotemporal maps --> current_3D_peak_voxels set in sm_spatiotemp_mapping
            [h.current_3D_peak_voxels,h.current_3D_peak_idx,h.s1,h.p1,h.func3D,h.anat3D,h.func3D_midline]=bl_plot_lcmv_peak_img_FT_new(h.axes_3D_images, h.inv_soln(h.current_inv_soln).soln.P.img,...
                h.inv_soln(h.current_inv_soln).soln.ori,null_thresh,search_dist,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,parula(255),...
                min_max,vol,h.anatomy.sens,h.cfg.study.bl_bmf.vw_angle,pk_flag,1,h.inv_soln(h.current_inv_soln).leadfield.pos,inside_idx,h.cfg.study.bl_bmf.slice_orient,sFaceAlpha);
        elseif h.radio_find_spatiotemp_peaks.Value == 0 && h.btn_3D_plot_peak_waves.Value == 1        % plot patiotemp map at time point selected in time domain waves
            %% Plot Whole Brain Soln
            [h.current_3D_peak_voxels,h.current_3D_peak_idx,h.s1,h.p1,h.func3D,h.anat3D,h.func3D_midline]=bl_plot_lcmv_peak_img_FT_new(h.axes_3D_images, h.inv_soln(h.current_inv_soln).soln.P.img,...
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
                    case {'LCMV (FT)' 'sLORETA (FT)' 'MNE (FT)' 'dics (FT)' 'pcc  (FT)'  'SAM (FT)' 'LCMV (BST)' 'sLORETA (BST)' 'MNE (BST)'}
                        inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
                        h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).soln.P.img_org(inside_idx);
                end
            else
                switch h.inv_soln(h.current_inv_soln).Type
                    case {'SIA' 'MIA'}
                    case {'LCMV (FT)' 'sLORETA (FT)' 'MNE (FT)' 'dics (FT)' 'pcc  (FT)'  'SAM (FT)' 'LCMV (BST)' 'sLORETA (BST)' 'MNE (BST)'}
                        if size(h.inv_soln(h.current_inv_soln).soln.P.img,1)>size(h.inv_soln(h.current_inv_soln).leadfield.H,3)
                            inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
                            h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).soln.P.img(inside_idx);
                        end
                end
            end
            [h.current_3D_peak_voxels,h.current_3D_peak_idx,h.s1,h.p1,h.func3D,h.anat3D,h.func3D_midline]=bl_plot_lcmv_peak_img_FT_new(h.axes_3D_images, h.inv_soln(h.current_inv_soln).soln.P.img,...
                h.inv_soln(h.current_inv_soln).soln.ori,null_thresh,search_dist,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,parula(255),...
                min_max,vol,h.anatomy.sens,h.cfg.study.bl_bmf.vw_angle,pk_flag,1,h.inv_soln(h.current_inv_soln).leadfield.pos,inside_idx,h.cfg.study.bl_bmf.slice_orient,sFaceAlpha);
        end
        try; set(h.s1,'Visible','off'); set(h.p1,'Visible','off'); end
        
        try sm_view_map_slices(); catch; end
        alpha_gain = [.25 1 .5];
        for a=1:length(h.func3D);  h.func3D(a).FaceAlpha=h.slider_3D_transparency_func.Value; end
        for a=1:length(h.anat3D);  h.anat3D(a).FaceAlpha=h.slider_3D_transparency_anat.Value*alpha_gain(a); end
        h.anat3D(1).FaceColor = h.scalp_clr;
        title(sprintf('%s',h.inv_soln(h.current_inv_soln).Type));
        
    case 'Cortical Surface'    % Cortical Surface
        
        clear vol;
        if length(h.anatomy.mesh_volumes)>=5
            vol = h.anatomy.mesh_volumes([1 4 4]);
        else
            vol(1) = h.inv_soln(h.current_inv_soln).headmodel_mesh;
            vol(2) = h.inv_soln(h.current_inv_soln).headmodel_mesh;
            vol(3) = h.inv_soln(h.current_inv_soln).headmodel_mesh;
        end
        
        alpha_gain = [.25 1 .5];
        vol(1).FaceColor = h.scalp_clr; % scalp
        vol(2).FaceColor = h.brain_clr; % brain
        vol(1).FaceAlpha = h.slider_3D_transparency_anat.Value*alpha_gain(1);
        vol(2).FaceAlpha = h.slider_3D_transparency_anat.Value*alpha_gain(1);
        vol(3).FaceAlpha = h.slider_3D_transparency_func.Value;
        vol(3).img = h.inv_soln(h.current_inv_soln).soln.P.img;
        vol(3).img(vol(3).img<null_thresh)=nan;
        
        opt.vol_nums=1:3;
        opt.caxis = h.axes_3D_images.CLim;
        
        if h.radio_find_spatiotemp_peaks.Value == 1 % Find peaks using spatiotemporal search across active interval
            sm_spatiotemp_mapping;
        else        % plot using map created by inverse solution
            % Peak locations
            voxel_vals=[h.inv_soln(h.current_inv_soln).leadfield.voxel_pos, vol(3).img];
            thresh_val = null_thresh;
            [h.current_3D_peak_voxels,h.current_3D_peak_idx]=BRANELab_find_peak_voxel_thresh(voxel_vals,thresh_val,search_dist);
        end
        try h.current_3D_peak_voxels = [h.current_3D_peak_voxels h.current_3D_peak_idx];
        catch; h.current_3D_peak_voxels = [h.current_3D_peak_voxels h.current_3D_peak_idx'];
        end
        
        
%         axes(h.axes_3D_images); cla; %view(cur_vw);
        h.func3D=bl_plot_mesh(vol,opt);
        %         h.axes_3D_images.SortMethod='childorder';
        h.axes_3D_images.SortMethod='depth';
        title(sprintf('%s',h.inv_soln(h.current_inv_soln).Type));
end

%% Do not update these once called set here - only to be set during run_source_modeling.m
if h.run_inv_soln_flag == 1 || h.radio_update_peaks.Value==1 % Update only on first instance of run_source_modeling.m
    h.inv_soln(h.current_inv_soln).peak_voxels = h.current_3D_peak_voxels;  % setting these to be re-ordered in [true_idx fa_idx]
    h.inv_soln(h.current_inv_soln).peak_idx = h.current_3D_peak_idx;
end
%% setting slider value to zero
h.slider_3D_image_thresh.Value = 0;

%% Setting colormap limit
try
    h.axes_3D_images.CLim = str2num(h.edit_3D_min_max.String);
catch
    h.axes_3D_images.CLim = [0 1]; h.edit_3D_min_max.String='0 1';
end

h.axes_3D_images.Colormap=eval(sprintf('%s(255)',h.menu_InvMap_colormap.String{h.menu_InvMap_colormap.Value}));
hold on; h.mrk_size=150; h.mrk_size2=50;
% true locs
if isfield(h,'sim_data')
    if isfield(h.sim_data,'cfg')
        true_locs = h.sim_data.cfg.source.vx_locs;
    else
        true_locs = h.cfg.source.vx_locs;
    end
else
    true_locs = h.cfg.source.vx_locs;
end
%% Scatter plot of True locations as baskets
if isfield(h,'map3D_true_locs'); h = rmfield(h,'map3D_true_locs'); end % clear field for overwriting
for v=1:length(h.cfg.source.vx_idx)
    %     h.map3D_true_locs(1,v) = scatter3(h.sim_data.cfg.source.vx_locs(seed_idx(v),1),h.sim_data.cfg.source.vx_locs(seed_idx(v),2),h.sim_data.cfg.source.vx_locs(seed_idx(v),3),'+','MarkerEdgeColor',h.cfg.source.src_clr(v,:),'sizedata',mrk_size,'linewidth',ln_wdth,'Visible',h.radio_3D_true_locs.Value==1);
    %     h.map3D_true_locs(2,v) =scatter3(h.sim_data.cfg.source.vx_locs(seed_idx(v),1),h.sim_data.cfg.source.vx_locs(seed_idx(v),2),h.sim_data.cfg.source.vx_locs(seed_idx(v),3),'s','MarkerEdgeColor',h.cfg.source.src_clr(v,:),'sizedata',mrk_size,'linewidth',ln_wdth2,'Visible',h.radio_3D_true_locs.Value==1);
    h.map3D_true_locs(1,v) = scatter3(true_locs(seed_idx(v),1),true_locs(seed_idx(v),2),true_locs(seed_idx(v),3),'+','MarkerEdgeColor',h.cfg.source.src_clr(v,:),'sizedata',h.mrk_size,'linewidth',ln_wdth,'Visible',h.radio_3D_true_locs.Value==1);
    h.map3D_true_locs(2,v) =scatter3(true_locs(seed_idx(v),1),true_locs(seed_idx(v),2),true_locs(seed_idx(v),3),'s','MarkerEdgeColor',h.cfg.source.src_clr(v,:),'sizedata',h.mrk_size,'linewidth',ln_wdth2,'Visible',h.radio_3D_true_locs.Value==1);
end

%% Plotting Scatter Locations and Orientations after searching for Hits & Reordering peaks to be in order of those nearest to source 1, 2, 3 locations
if isempty(h.current_3D_peak_voxels) % no peaks found
    text(0,0,0,'No Peak Sources Found'); h.listbox_peaks_found.String = ''; h.listbox_peaks_found.Value = 1;
%     view(h.axes_3D_images,cur_vw); axis(h.axes_3D_images,'tight');
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
    sm_calc_localizer_performance;
else
    
    if h.radio_find_spatiotemp_peaks.Value == 1
        % Do Nothing because nearest peaks already found
            h.current_inv_soln_show_peak_idx = 1:length(h.current_3D_peak_idx);
            sm_search_for_hits('initial search'); % reorders "h.current_3D_peak_voxels" after performing search based on user's selection of 'Nearest', 'Wave Error', 'Wave Correlation'
    else
        if h.run_inv_soln_flag ==1  % called during by run_source_modeling.m
            h.current_inv_soln_show_peak_idx = 1:length(h.current_3D_peak_idx);
            sm_search_for_hits('initial search'); % reorders "h.current_3D_peak_voxels" after performing search based on user's selection of 'Nearest', 'Wave Error', 'Wave Correlation'
        else
            h.slider_3D_image_thresh.Value = 0; %h.inv_soln(h.current_inv_soln).soln.plot_thresh;
            sm_search_for_hits('slider thresh');
        end
    end
    sm_calc_localizer_performance();
    
    %% peak locs
    if isempty(h.current_3D_peak_voxels) % no peaks found
        text(0,0,0,'No Peak Sources Found'); h.listbox_peaks_found.String = ''; h.listbox_peaks_found.Value = 1;
%         view(h.axes_3D_images,cur_vw); axis(h.axes_3D_images,'tight');
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
        
        sm_plot_3Dmap_locs();
        try
            update_image_thresh_txt;
        catch
            sm_plot_results_loc_performance; update_image_thresh_txt;
        end
            
        update_listbox_peaks_found();
        bs_plot_peak_waves;
        
    end
end

h.slider_3D_image_thresh.Value = h.inv_soln(h.current_inv_soln).soln.plot_min_max(1);
if h.slider_3D_image_thresh.Value<h.slider_3D_image_thresh.Min
    h.slider_3D_image_thresh.Value = h.slider_3D_image_thresh.Min;
elseif h.slider_3D_image_thresh.Value>h.slider_3D_image_thresh.Max
    h.slider_3D_image_thresh.Value = h.slider_3D_image_thresh.Max;
end



%% update FC contrasts if they exist
if isfield(h.inv_soln(h.current_inv_soln),'plv_contrasts') && isfield(h.inv_soln(h.current_inv_soln),'plv_seed_idx') && any(~isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.Hits))
    
    if ~isempty(h.inv_soln(h.current_inv_soln).plv_seed_idx)
        sm_update_plv_contrast_order();
    else
        h.listbox_plv_contrasts.String = '';
    end
else
    h.listbox_plv_contrasts.String = '';
end

if isfield(h.inv_soln(h.current_inv_soln),'TFR_results')
    if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results) && any(~isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.Hits))
        sm_plot_tfr_connectivity;
    else
        h.axes_inv_soln_tfr.Visible = 'off'; for a=1:length(h.axes_inv_soln_tfr.Children); h.axes_inv_soln_tfr.Children(a).Visible='off'; end
        %% need to delete colorbar
        if isfield(h,'colorbar_axes_inv_soln_tfr'); if isvalid(h.colorbar_axes_inv_soln_tfr); delete(h.colorbar_axes_inv_soln_tfr); end; end
        h.radio_inv_plot_peak_tfr_connectivity.Value = 0;
    end
else
    
end

% Turning lighting back on if it was on to begin with
% h.toggle_light_OnOff.Value = light_flag;
toggle_light_OnOff_Callback;
set_3D_transparency;

if h.run_inv_soln_flag ==1  % called during by run_source_modeling.m
    h.slider_3D_image_thresh.Value = 0;
else
    h.slider_3D_image_thresh.Value = h.inv_soln(h.current_inv_soln).soln.plot_thresh;
end

set_3D_image_thresh(); % updating maps based on slider

%% updating Localization Performance Results Graphs
if ~isfield(h.inv_soln(h.current_inv_soln),'classifier_results')
    sm_calc_results_loc_performance();
else
    if isempty(h.inv_soln(h.current_inv_soln).classifier_results)
        sm_calc_results_loc_performance();
    end
end
sm_plot_results_loc_performance();

% if isvalid(hm)
%     delete(hm)
% end

%% returning selected peaks back to original
h.listbox_peaks_found.Value = sel_idx; 

axis(h.axes_3D_images,'off');
% h.slider_3D_image_thresh.Enable = 'on';
h.panel_3D_image_plot_msg.Visible = 'off';

%% View
view(h.axes_3D_images,cur_vw); axis(h.axes_3D_images,'tight');

