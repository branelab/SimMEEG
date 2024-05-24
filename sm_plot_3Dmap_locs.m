function sm_plot_3Dmap_locs(varargin)
global h

%% Initializing some variables for plotting
min_max = h.inv_soln(h.current_inv_soln).soln.plot_min_max; %[min(h.inv_soln(h.current_inv_soln).soln.P.img) max(h.inv_soln(h.current_inv_soln).soln.P.img)];
if min_max(1)==min_max(2); min_max(2)=min_max(1)+1;
elseif min_max(1)>min_max(2); min_max(2)=min_max(1)+1;
end

h.ln_clr = bsxfun(@mtimes,ones(length(h.current_3D_peak_idx),3),h.FA_clr);  %lines(length(h4.current_3D_peak_idx));
h.ln_clr(1:size(h.cfg.source.src_clr,1),:) = h.cfg.source.src_clr;
cur_vw = h.axes_3D_images.View;

%% deleting existing scatter of peak_voxels 
if isfield(h,'map3D_peak_locs')
    if ~isempty(h.map3D_peak_locs)
        if any(isvalid(h.map3D_peak_locs))
        delete(h.map3D_peak_locs)
        delete(h.map3D_peak_ori)
        end
    end
end

%% Plotting locs
h.map3D_peak_locs =[];
if any(~isnan(h.current_3D_peak_idx))
    p_idx = find(~isnan(h.current_3D_peak_idx));
    for v = 1:length(p_idx)
        h.map3D_peak_locs(p_idx(v)) = scatter3(h.axes_3D_images,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx(p_idx(v)),1),...
            h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx(p_idx(v)),2),...
            h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx(p_idx(v)),3),...
            'filled','ko','sizedata',h.mrk_size2,'linewidth',3,'Visible',h.radio_3D_peak_locs.Value==1);
        set(h.map3D_peak_locs(p_idx(v)),'CData',h.ln_clr(p_idx(v),:));
        
        if v<=length(h.cfg.source.vx_idx) % hit peaks
            set(h.map3D_peak_locs(p_idx(v)),'MarkerFaceAlpha',1);
            set(h.map3D_peak_locs(p_idx(v)),'MarkerEdgeAlpha',1); % sets transparency
        else
            set(h.map3D_peak_locs(p_idx(v)),'MarkerFaceAlpha',h.false_positive_FaceAlpha);
            set(h.map3D_peak_locs(p_idx(v)),'MarkerEdgeAlpha',h.false_positive_FaceAlpha); % sets transparency
        end
        
    end
    h.map3D_peak_locs = handle(h.map3D_peak_locs);
end

%% plotting orientations
lf_pos =  h.inv_soln(h.current_inv_soln).leadfield.voxel_pos;
if isfield(h.inv_soln(h.current_inv_soln).leadfield,'voxel_res')
vx_res = h.inv_soln(h.current_inv_soln).leadfield.voxel_res; %round(min(pdist(lf_pos,'euclidean')));  %max(lf_pos(2,:)-lf_pos(1,:));
if isempty(vx_res); vx_res = 5; end % default
else
vx_res = 5;    
end
ori = h.inv_soln(h.current_inv_soln).soln.ori;
if size(ori,2)<3 % MEG solns have 2 dipoles for eLORETA
    ori(:,3) = zeros(size(ori,1),1);
end
h.map3D_peak_ori =[];

% plotting hit ori
%         if ~isempty(h.current_3D_peak_idx_hit_idx)
if any(~isnan(h.current_3D_peak_idx))
    lf_idx = h.current_3D_peak_idx(find(~isnan(h.current_3D_peak_idx)));
    p_idx = find(~isnan(h.current_3D_peak_idx));
    for v=1:length(lf_idx)
        amp_gain =abs( ( h.inv_soln(h.current_inv_soln).soln.P.img(lf_idx(v)) / min_max(2) ) );
        vx_pos=lf_pos(lf_idx(v),:);
        ori_pos=vx_pos+(4*amp_gain*vx_res*ori(lf_idx(v),:));
        h.map3D_peak_ori(p_idx(v))=plot3(h.axes_3D_images,[vx_pos(1) ori_pos(1)],[vx_pos(2) ori_pos(2)],[vx_pos(3) ori_pos(3)],'color',h.ln_clr(p_idx(v),:),'linewidth',2);
        if p_idx(v)<=length(h.sim_data.cfg.source.vx_idx) % hit peaks
            xh = handle(h.map3D_peak_ori(p_idx(v))); xh.Color(4)=1;    % sets transparency
        else
            xh = handle(h.map3D_peak_ori(p_idx(v))); xh.Color(4)=h.false_positive_FaceAlpha;    % sets transparency
        end
    end
    h.map3D_peak_ori = handle(h.map3D_peak_ori);

    if h.monte_carlo_flag==0; drawnow; end

end


%%
if diff(h.inv_soln(h.current_inv_soln).soln.plot_min_max)<=0; h.inv_soln(h.current_inv_soln).soln.plot_min_max = [0 2]; end
h.axes_3D_images.CLim = h.inv_soln(h.current_inv_soln).soln.plot_min_max;
% h.edit_3D_min_max.String = num2str(h.axes_3D_images.CLim); %sprintf('%.3f %.3f',h.axes_3D_images.CLim);
h.edit_3D_min_max.String = sprintf('%.3f %.3f',h.axes_3D_images.CLim);

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
