function [p1, s1] = bl_plot_FC_graph(fc_data, vol, vx_locs, fc_contrasts, fc_thresh, fc_scale)
% function [p1, s1] = bl_plot_FC_graph(fc_data, vol, vx_locs, fc_contrasts, fc_thresh, fc_scale);
%
% INPUT:
%   fc_data = FC data for plotting lines on mesh brain volume
%   vol     = mesh volume of brain/skull (see bl_plot_mesh.m)
%   vx_locs = voxel positions (x y z) in mm
%   fc_contrasts = voxel indices for FC contrasts 
%   fc_thresh = threshold for FC map
%   
p1=[]; s1=[];
opt.vol_nums = 1:length(vol);

fc_locs = vx_locs(fc_contrasts,:);
s1 = scatter3(fc_locs(:,1),fc_locs(:,2),fc_locs(:,3),'ko','SizeData',20,'LineWidth',1,'MarkerFaceColor','flat');
fc_locs2 = reshape(fc_locs,[size(fc_locs,1)/2 2 3]); % [plv_idx x contrast x XYZ]

% plotting new FC graph
for v = 1:length(fc_data)
    if  abs(fc_data(v)) > fc_thresh
        cx_locs = squeeze(fc_locs2(v,:,:));
        %                         ln_width = abs(fc_data(v)*15);
        ln_width = abs(fc_data(v))/max(fc_scale)*4;
        if fc_data(v)>0; ln_clr = 'r'; elseif fc_data(v)<=0; ln_clr = 'b'; end % ERS=red and ERD=blue
        p1(v) = plot3(cx_locs(:,1),cx_locs(:,2),cx_locs(:,3),'LineWidth',ln_width,'Color',ln_clr);

        %                         alpha_val=abs(fc_data(v)*2); if alpha_val>1; alpha_val=1; end
        alpha_val=abs(fc_data(v))/max(fc_scale); if alpha_val>1; alpha_val=1; end
%         if ~isempty(p1(v)); p1(v).Color(4) = alpha_val; end
    end
end




