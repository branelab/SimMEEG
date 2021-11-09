function bs_calc_errors_inv_soln(varargin)
global h

%% Location Errors
true_locs = h.cfg.source.vx_locs;

if h.menu_inv_headmodel.Value == 1 && h.menu_sens_type.Value == 1   % volume MEG
    peak_locs = h.anatomy.leadfield_meg_vol.voxel_pos(h.current_3D_peak_idx(1:3),:);
elseif h.menu_inv_headmodel.Value == 2 && h.menu_sens_type.Value == 1   % cortex MEG
    peak_locs = h.anatomy.leadfield_meg_cortex.voxel_pos(h.current_3D_peak_idx(1:3),:);
elseif h.menu_inv_headmodel.Value == 1 && h.menu_sens_type.Value == 2   % volume EEG
    peak_locs = h.anatomy.leadfield_eeg_vol.voxel_pos(h.current_3D_peak_idx(1:3),:);
elseif h.menu_inv_headmodel.Value == 2 && h.menu_sens_type.Value == 2   % cortex EEG
    peak_locs = h.anatomy.leadfield_eeg_cortex.voxel_pos(h.current_3D_peak_idx(1:3),:);
else
    peak_locs = h.anatomy.leadfield.voxel_pos(h.current_3D_peak_idx(1:3),:);
end


h.current_mse_locs = sqrt( nanmean( (true_locs-peak_locs).^2 , 2));
b1=bar(h.axes_invSoln_errors_locs,h.current_mse_locs,'FaceColor','Flat','CData',h.src_clr); title(h.axes_invSoln_errors_locs,sprintf('Location\nError (mm)'));
xtips2 = double(b1.XEndPoints); ytips2 = double(b1.YEndPoints); labels2 = string(round(b1.YData)); 
t1 = text(h.axes_invSoln_errors_locs,xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);

%% Orientation Errors
try
true_ori = h.cfg.source.vx_ori;
peak_ori = nan(size(true_ori));
peak_ori2 = h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(1:3),:);
peak_ori(1:size(peak_ori2,1),1:size(peak_ori2,2)) = peak_ori2;
h.current_mse_ori = rad2deg(sqrt( nanmean( (true_ori-peak_ori).^2 , 2)));
b2=bar(h.axes_invSoln_errors_ori,h.current_mse_ori,'FaceColor','Flat','CData',h.src_clr); title(h.axes_invSoln_errors_ori,sprintf('Orientation\nError (deg)')); 
catch
    h.axes_invSoln_errors_ori.clo; text(h.axes_invSoln_errors_ori,.1,.25, sprintf('Calculation\nERROR')); 
end
xtips2 = double(b2.XEndPoints); ytips2 = double(b2.YEndPoints); labels2 = string(round(b2.YData)); 
t1 = text(h.axes_invSoln_errors_ori,xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);


%% waveform errors - as residual variances
act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
% h.current_mse_evk_waves = sqrt( nanmean( (h.norm_true_swf(act_samps,:) - h.current_norm_peak_swf(act_samps,1:3)).^2, 1) );
% b3=bar(h.axes_invSoln_errors_waves,h.current_mse_evk_waves ,'FaceColor','Flat','CData',h.src_clr); title(h.axes_invSoln_errors_waves,sprintf('Peak Wave\nError (mse)')); 
% changing this to a percentage of explained variance - "the sum of squares of the unexplained signal"
% switch h.inv_soln(h.current_inv_soln).Type
%     case {'SPA' 'SIA' 'MIA' 'LCMV' 'sLORETA'  'dics' 'pcc' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'} % scalar beamformers
%         d1 = real(h.norm_true_swf(act_samps,:)); d2 = real(h.current_norm_peak_swf(act_samps,1:3));
%         h.current_mse_evk_waves = real(100 * (sum((d1-d2).^2) ./ sum(d1.^2) )); % same as field Trip's Calculation
%     case {'eLORETA' 'MNE'}  % vector beamformers
%         d1 = real(h.norm_true_swf(act_samps,:)); d2 = real(h.current_norm_peak_swf(act_samps,1:3));
%         switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
%             case 'RMS'
%                 h.current_mse_evk_waves = real(100 * (sum((d1-d2).^2) ./ sum(d1.^2) )); % same as field Trip's Calculation
%             case 'Max'
%                 h.current_mse_evk_waves = real(100 * (sum((d1-d2).^2) ./ sum(d1.^2) )); % same as field Trip's Calculation
%             case 'avg.pow'
%                 %                 h.current_mse_evk_waves = real(100*sum( abs(h.norm_true_swf(act_samps,:) - h.current_norm_peak_swf(act_samps,1:3))) ./ sum(abs(h.norm_true_swf(act_samps,:))));     % same as field Trip's Calculation
% %                 h.current_mse_evk_waves = real(100 * (sum((d1-(d2.^.5)).^2) ./ sum(d1.^2) )); % need to sqrt avg.pow before sum of squares
%                 h.current_mse_evk_waves = real(100 * (sum((d1-d2).^2) ./ sum(d1.^2) )); % same as field Trip's Calculation
%         end
% end

% going back to mean-squared error --> residual variance was yielding >100%, also not consistent with mean-squared-error results
 d1 = real(h.norm_true_swf(act_samps,:)); d2 = real(h.current_norm_peak_swf(act_samps,1:3));
 h.current_mse_evk_waves = nanmean((d1-d2).^2); 

b3=bar(h.axes_invSoln_errors_waves,h.current_mse_evk_waves ,'FaceColor','Flat','CData',h.src_clr); title(h.axes_invSoln_errors_waves,sprintf('Wave Error\n(MSE)')); 
xtips2 = double(b3.XEndPoints); ytips2 = double(b3.YEndPoints); labels2 = string(round(b3.YData)); 
t1 = text(h.axes_invSoln_errors_waves,xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);

%% Ylims
mx=30; if max(h.current_mse_locs)*1.15>mx; mx = max(h.current_mse_locs)*1.15; end
h.axes_invSoln_errors_locs.XLabel.String = 'Peaks'; h.axes_invSoln_errors_locs.YLim = [0 mx];
mx=30; if max(h.current_mse_ori)*1.15>mx; mx = max(h.current_mse_ori)*1.15; end
h.axes_invSoln_errors_ori.XLabel.String = 'Peaks';  h.axes_invSoln_errors_ori.YLim = [0 mx];

mx=50; if max(h.current_mse_evk_waves)*1.15>mx; mx = max(h.current_mse_evk_waves)*1.15; end
h.axes_invSoln_errors_waves.XLabel.String = 'Peaks'; h.axes_invSoln_errors_waves.YLim = [0 mx];


%% halfmax widths
h.axes_invSoln_halfmax_width.clo;
h.inv_soln(h.current_inv_soln).half_width = [];
h.inv_soln(h.current_inv_soln).xslices =[];
h.inv_soln(h.current_inv_soln).yslices =[];
h.inv_soln(h.current_inv_soln).zslices =[];

switch h.inv_soln(h.current_inv_soln).headmodel_type  % only calulating half-width when volume head model - can not calcualte spread function for cortically constrained
    case 'Whole Brain'
        h.axes_invSoln_halfmax_width.YLim = [0 30];
        switch h.inv_soln(h.current_inv_soln).Type
            case {'SPA' 'LCMV' 'eLORETA' 'dics' 'pcc' 'sLORETA' 'MNE'}
                [h.inv_soln(h.current_inv_soln).half_width,h.inv_soln(h.current_inv_soln).xslices,h.inv_soln(h.current_inv_soln).yslices,h.inv_soln(h.current_inv_soln).zslices]=bs_calc_halfmax_spread(h.inv_soln(h.current_inv_soln).soln.P.img,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,h.current_3D_peak_idx);
                b4=bar(h.axes_invSoln_halfmax_width, [h.inv_soln(h.current_inv_soln).half_width(1:3).average] ,'FaceColor','Flat','CData',h.src_clr);
                title(h.axes_invSoln_halfmax_width,sprintf('Half-Max\n Width (mm)')); h.axes_invSoln_halfmax_width.XLabel.String = 'Peaks';
                mx=30; if max([h.inv_soln(h.current_inv_soln).half_width(1:3).average])*1.15>mx; mx = max([h.inv_soln(h.current_inv_soln).half_width(1:3).average])*1.15; end
                h.axes_invSoln_halfmax_width.YLim = [0 mx];
                xtips2 = double(b4.XEndPoints); ytips2 = double(b4.YEndPoints); labels2 = string(round(b4.YData));
                t1 = text(h.axes_invSoln_halfmax_width,xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
                
            case {'SIA', 'MIA' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}
                text(h.axes_invSoln_halfmax_width,0.1,15,sprintf('Not Applicable\nfor %s',h.inv_soln(h.current_inv_soln).Type));
                h.axes_invSoln_halfmax_width.YLim = [0 30];
        end
        
    case 'Cortical Surface'
        text(h.axes_invSoln_halfmax_width,0.1,15,sprintf('Not Applicable\nfor Cortical\nSurface',h.inv_soln(h.current_inv_soln).Type));
        h.axes_invSoln_halfmax_width.YLim = [0 30];
end



%% Peaks 
h.inv_soln(h.current_inv_soln).peak_idx = h.current_3D_peak_idx;
h.inv_soln(h.current_inv_soln).peak_voxels = h.current_3D_peak_voxels;

