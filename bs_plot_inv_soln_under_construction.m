function bs_plot_inv_soln(varargin)
global h

if h.monte_carlo_flag == 1

else
    h.waitfor_panel.Visible='off';
end

if h.radio_find_spatiotemp_peaks.Value == 1
    hm = msgbox(sprintf('Plotting Inverse Solution\n\nSpatiotemporal mapping')); WinOnTop(hm);
else
    hm = msgbox(sprintf('Plotting Inverse Solution\n\nRegular mapping')); WinOnTop(hm);
end


%%
cur_vw = h.axes_3D_images.View;
seed_idx = 1:3;  ln_wdth = 1; ln_wdth2 = 1;
h.current_3D_peak_idx =[]; h.current_3D_peak_voxels = []; h.current_norm_peak_swf =[]; h.current_peak_swf = [];

min_max = h.inv_soln(h.current_inv_soln).soln.plot_min_max; %[min(h.inv_soln(h.current_inv_soln).soln.P.img) max(h.inv_soln(h.current_inv_soln).soln.P.img)];

if min_max(1)==min_max(2); min_max(2)=min_max(1)+1;
elseif min_max(1)>min_max(2); min_max(2)=min_max(1)+1;
end

null_thresh = h.inv_soln(h.current_inv_soln).soln.plot_thresh;

%% for plotting
switch h.inv_soln(h.current_inv_soln).headmodel_type
    case 'Whole Brain' % whole brain
        
        axes(h.axes_3D_images); cla; view(cur_vw);
        
        switch h.inv_soln(h.current_inv_soln).Type
            case {'SPA'}    % BRANE Lab beamformers
                search_dist = 15;
            case {'SIA','MIA'}    % BRANE Lab beamformers
                search_dist = 15;
           case {'sMCMV','bRAPBeam','TrapMUSIC'}    % Alex Moiseev's beamformers
                search_dist = 15;
            case {'LCMV' 'dics' 'pcc'}    % Field Trips inverse solutions
                search_dist = 15;
            case {'MNE','eLORETA','sLORETA'}    % Field Trips inverse solutions
                search_dist = 15;
        end
        
        % parameters for plotting functional maps
        h.cfg.study.bl_bmf.vw_angle = h.axes_3D_images.View;
        pk_flag=0; % just plotting the maps
        inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
        sFaceAlpha = h.slider_3D_transparency_func.Value;
        %         vol = h.inv_soln(h.current_inv_soln).headmodel_mesh;
        vol = h.anatomy.mesh_volumes([1 4]);
        
        
        if h.radio_find_spatiotemp_peaks.Value == 1 % Find peaks using spatiotemporal search across active interval
            sm_spatiotemp_mapping;
            h.inv_soln(h.current_inv_soln).soln.plot_min_max = [min(h.inv_soln(h.current_inv_soln).soln.P.img) max(h.inv_soln(h.current_inv_soln).soln.P.img)];
            min_max = h.inv_soln(h.current_inv_soln).soln.plot_min_max;
            null_thresh = h.slider_3D_image_thresh.Value;
            %% Plot Whole Brain Soln combined spatiotemporal maps --> current_3D_peak_voxels set in sm_spatiotemp_mapping
            [~,~,h.s1,h.p1,h.func3D,h.anat3D]=bl_plot_lcmv_peak_img_FT_new(h.inv_soln(h.current_inv_soln).soln.P.img,...
                h.inv_soln(h.current_inv_soln).soln.ori,null_thresh,search_dist,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,parula(255),...
                min_max,vol,h.anatomy.sens,h.cfg.study.bl_bmf.vw_angle,pk_flag,1,h.inv_soln(h.current_inv_soln).leadfield.pos,inside_idx,h.cfg.study.bl_bmf.slice_orient,sFaceAlpha);
        else        % plot using map created by inverse solution
            
            %% Plot Whole Brain Soln
            if isfield(h.inv_soln(h.current_inv_soln).soln.P,'img_org')     % reverting back to original inv_soln image
                switch h.inv_soln(h.current_inv_soln).Type
                    case {'SIA' 'MIA'}
%                         h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).soln.P.nulled_img;
%                         h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).soln.P.img_org;
                        % overlaying MCMV found maps onto the nulled MCMV map
                        img = h.inv_soln(h.current_inv_soln).soln.P.nulled_img; 
                        img(h.inv_soln(h.current_inv_soln).soln.MCMV_idx) = h.inv_soln(h.current_inv_soln).soln.P.img_org(h.inv_soln(h.current_inv_soln).soln.MCMV_idx);
                        h.inv_soln(h.current_inv_soln).soln.P.img = img; 
                        
                    case {'SPA' 'LCMV' 'eLORETA' 'sLORETA' 'MNE' 'dics' 'pcc' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}
                        h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).soln.P.img_org;
                end
            end
            [h.current_3D_peak_voxels,h.current_3D_peak_idx,h.s1,h.p1,h.func3D,h.anat3D]=bl_plot_lcmv_peak_img_FT_new(h.inv_soln(h.current_inv_soln).soln.P.img,...
                h.inv_soln(h.current_inv_soln).soln.ori,null_thresh,search_dist,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,parula(255),...
                min_max,vol,h.anatomy.sens,h.cfg.study.bl_bmf.vw_angle,pk_flag,1,h.inv_soln(h.current_inv_soln).leadfield.pos,inside_idx,h.cfg.study.bl_bmf.slice_orient,sFaceAlpha);
            
        end
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
            thresh_limit = 15;
            [h.current_3D_peak_voxels,h.current_3D_peak_idx]=BRANELab_find_peak_voxel_thresh(voxel_vals,thresh_val,thresh_limit);
        end
        try h.current_3D_peak_voxels = [h.current_3D_peak_voxels h.current_3D_peak_idx];
        catch; h.current_3D_peak_voxels = [h.current_3D_peak_voxels h.current_3D_peak_idx'];
        end
        
        
        axes(h.axes_3D_images); cla; view(cur_vw);
        h.func3D=bl_plot_mesh(vol,opt);
        %         h.axes_3D_images.SortMethod='childorder';
        h.axes_3D_images.SortMethod='depth';
        title(sprintf('%s',h.inv_soln(h.current_inv_soln).Type));
end

try
    h.axes_3D_images.CLim = str2num(h.edit_3D_min_max.String);
catch
    h.axes_3D_images.CLim = [0 1]; h.edit_3D_min_max.String='0 1';
end


h.axes_3D_images.Colormap=jet(255);

hold on; mrk_size=150; mrk_size2=50;
% true locs
if isfield(h,'sim_data')
    if isfield(h.sim_data,'cfg')
        vx_locs = h.sim_data.cfg.source.vx_locs;
    else
        vx_locs = h.cfg.source.vx_locs;
    end
else
    vx_locs = h.cfg.source.vx_locs;
end
%% Scatter plot of True locaitons as baskets
for v=1:3
    %     h.map3D_true_locs(1,v) = scatter3(h.sim_data.cfg.source.vx_locs(seed_idx(v),1),h.sim_data.cfg.source.vx_locs(seed_idx(v),2),h.sim_data.cfg.source.vx_locs(seed_idx(v),3),'+','MarkerEdgeColor',h.src_clr(v,:),'sizedata',mrk_size,'linewidth',ln_wdth,'Visible',h.radio_3D_true_locs.Value==1);
    %     h.map3D_true_locs(2,v) =scatter3(h.sim_data.cfg.source.vx_locs(seed_idx(v),1),h.sim_data.cfg.source.vx_locs(seed_idx(v),2),h.sim_data.cfg.source.vx_locs(seed_idx(v),3),'s','MarkerEdgeColor',h.src_clr(v,:),'sizedata',mrk_size,'linewidth',ln_wdth2,'Visible',h.radio_3D_true_locs.Value==1);
    h.map3D_true_locs(1,v) = scatter3(vx_locs(seed_idx(v),1),vx_locs(seed_idx(v),2),vx_locs(seed_idx(v),3),'+','MarkerEdgeColor',h.src_clr(v,:),'sizedata',mrk_size,'linewidth',ln_wdth,'Visible',h.radio_3D_true_locs.Value==1);
    h.map3D_true_locs(2,v) =scatter3(vx_locs(seed_idx(v),1),vx_locs(seed_idx(v),2),vx_locs(seed_idx(v),3),'s','MarkerEdgeColor',h.src_clr(v,:),'sizedata',mrk_size,'linewidth',ln_wdth2,'Visible',h.radio_3D_true_locs.Value==1);
end


%% Reordering peaks to be in order of those nearest to source 1, 2, 3 locations
if isempty(h.current_3D_peak_voxels) % no peaks found
    text(0,0,0,'No Peak Sources Found');
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
    
    if h.radio_find_spatiotemp_peaks.Value == 1
        % Do Nothing because nearest peaks already found
    else
   % commented out code here is for finding MCMV peaks first and then if rest of sources but this doesn't work if MCMV has <3 found sources
%         if isfield(h.inv_soln(h.current_inv_soln).soln,'MCMV_idx')  % if SIA or MIA, first find closest sources using found MCMV idx then use nulled image peaks to find nearest sources
% 
%             p_locs = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.inv_soln(h.current_inv_soln).soln.MCMV_idx,:);
%             [v_idx]=find_nearest_voxel(p_locs,vx_locs);    % find nearest true_source for each peak source
%             diff_locs=[];
%             for vx=1:3; diff_locs(vx,:) = ( nanmean( (p_locs-vx_locs(vx,:)).^2,2)).^.5; end
%             for vx=1:3; vvx=find(v_idx==vx); [mx_val,mx] = min( diff_locs(vx,vvx) ); min_idx(vx)=vvx(mx); end
%             
%             v_idx = min_idx;
%             diff_idx = setxor(h.current_3D_peak_voxels(:,5),h.inv_soln(h.current_inv_soln).soln.MCMV_idx(v_idx));
%             v_idx = [h.inv_soln(h.current_inv_soln).soln.MCMV_idx(v_idx) diff_idx'];
%             
%             h.current_3D_peak_idx(1:length(v_idx)) = v_idx; % MCMV_idx found and re-ordered to be nearest to true sources
%             h.current_3D_peak_voxels = [h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(v_idx,:) h.inv_soln(h.current_inv_soln).soln.P.img(v_idx) v_idx'];  % re-ordered
%         else
%         [v_idx]=find_nearest_voxel(vx_locs,h.current_3D_peak_voxels(:,1:3));
        [v_idx]=find_nearest_voxel(h.current_3D_peak_voxels(:,1:3),vx_locs);    % find nearest true_source for each peak source
        % finding closest peak source to each true source
        diff_locs=[];
        for vx=1:3; diff_locs(vx,:) = ( nanmean( (h.current_3D_peak_voxels(:,1:3)-vx_locs(vx,:)).^2,2)).^.5; end
        for vx=1:length(v_idx); vvx=find(v_idx==vx); [mx_val,mx] = min( diff_locs(vx,vvx) ); min_idx(vx)=vvx(mx); end
%         diag(diff_locs(:,min_idx))
        
        v_idx = min_idx;
        diff_idx = setxor(1:size(h.current_3D_peak_voxels,1),v_idx);
        v_idx = [v_idx diff_idx];
        h.current_3D_peak_idx = h.current_3D_peak_voxels(v_idx,5); % re-ordered
        h.current_3D_peak_voxels = h.current_3D_peak_voxels(v_idx,:);  % re-ordered

%         end
        
    end
    
    %% peak locs
     ln_clr = bsxfun(@mtimes,ones(length(h.current_3D_peak_idx),3),[.9 .6 .3]*.5);  %lines(length(h4.current_3D_peak_idx));
     ln_clr(1:3,:) = h.src_clr;
    h.map3D_peak_locs(1) = scatter3(h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx,1),...
        h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx,2),...
        h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx,3),...
        'filled','ko','sizedata',mrk_size2,'linewidth',3,'Visible',h.radio_3D_peak_locs.Value==1);
    h.map3D_peak_locs(1).CData = ln_clr; 
    h.map3D_peak_locs(1).MarkerFaceAlpha = h.false_positive_FaceAlpha; h.map3D_peak_locs(1).MarkerEdgeAlpha = h.false_positive_FaceAlpha; 
    % plotting nearest peak source
    h.map3D_peak_locs(2) = scatter3(h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx(1:3),1),...
        h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx(1:3),2),...
        h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx(1:3),3),...
        'filled','ko','sizedata',mrk_size2,'linewidth',3,'Visible',h.radio_3D_peak_locs.Value==1);
    h.map3D_peak_locs(2).CData = ln_clr(1:3,:); 
    h.map3D_peak_locs(2).MarkerFaceAlpha = 1; h.map3D_peak_locs(2).MarkerEdgeAlpha = 1; 

    
    %% plotting orientations
    lf_pos =  h.inv_soln(h.current_inv_soln).leadfield.voxel_pos;
    vx_res = 5; %max(lf_pos(2,:)-lf_pos(1,:));
    ori = h.inv_soln(h.current_inv_soln).soln.ori;
    if size(ori,2)<3 % MEG solns have 2 dipoles for eLORETA
        ori(:,3) = zeros(size(ori,1),1);
    end
    h.map3D_peak_ori =[];
    
    for v=1:length(h.current_3D_peak_idx)
        amp_gain =abs( ( h.inv_soln(h.current_inv_soln).soln.P.img(h.current_3D_peak_idx(v)) / min_max(2) ) );
        vx_pos=lf_pos(h.current_3D_peak_idx(v),:);
        ori_pos=vx_pos+(2*amp_gain*vx_res*ori(h.current_3D_peak_idx(v),:));
        h.map3D_peak_ori(v)=plot3([vx_pos(1) ori_pos(1)],[vx_pos(2) ori_pos(2)],[vx_pos(3) ori_pos(3)],'color',ln_clr(v,:),'linewidth',2);
        if v<=3 % nearest peak sources
            xh = handle(h.map3D_peak_ori(v)); xh.Color(4)=1;    % sets transparency
        else
            xh = handle(h.map3D_peak_ori(v)); xh.Color(4)=h.false_positive_FaceAlpha;    % sets transparency
        end
    end
    
    %%
    h.axes_3D_images.CLim = h.inv_soln(h.current_inv_soln).soln.plot_min_max;
    h.edit_3D_min_max.String = num2str(h.axes_3D_images.CLim); %sprintf('%.3f %.3f',h.axes_3D_images.CLim);
    
    h.slider_3D_image_thresh.Min = h.axes_3D_images.CLim(1);
    h.slider_3D_image_thresh.Max = h.axes_3D_images.CLim(2);
    h.slider_3D_image_thresh.Value = h.inv_soln(h.current_inv_soln).soln.plot_thresh;
    
    view(h.axes_3D_images,cur_vw); axis(h.axes_3D_images,'tight');
    if isfield(h,'colorbar_3D')
        if isvalid(h.colorbar_3D)
            delete(h.colorbar_3D); h.colorbar_3D = colorbar(h.axes_3D_images,'Location','southoutside','Position',[.75 .35 .2 .03]);
        else
            h.colorbar_3D = colorbar(h.axes_3D_images,'Location','southoutside','Position',[.75 .35 .2 .03]);
        end
    else
        h.colorbar_3D = colorbar(h.axes_3D_images,'Location','southoutside','Position',[.75 .35 .2 .03]);
    end
    
    update_image_thresh_txt;
    update_listbox_peaks_found();
    bs_plot_peak_waves;
    
end

if h.slider_3D_image_thresh.Value<h.slider_3D_image_thresh.Min
    h.slider_3D_image_thresh.Value = h.slider_3D_image_thresh.Min;
elseif h.slider_3D_image_thresh.Value>h.slider_3D_image_thresh.Max
    h.slider_3D_image_thresh.Value = h.slider_3D_image_thresh.Max;
end

if isfield(h.inv_soln(h.current_inv_soln),'TFR_results')
    if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results)
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

if isvalid(hm)
    delete(hm)
end


