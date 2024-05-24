function bs_calc_errors_inv_soln(varargin)
global h

hit_idx = find(~isnan(h.current_peak_hit_lf_idx));
miss_idx = find(~isnan(h.current_peak_miss_lf_idx));
true_locs = h.sim_data.cfg.source.vx_locs;
peak_locs = nan(size(true_locs));
hit_lf_idx = nan(size(true_locs,1),1);
if max(hit_idx) > size(h.inv_soln(h.current_inv_soln).peak_voxels,1)    % less peak_voxels found than true sources
    hit_lf_idx(hit_idx) = h.current_3D_peak_voxels(hit_idx,5);
else
    hit_lf_idx(hit_idx) = h.inv_soln(h.current_inv_soln).peak_voxels(hit_idx,5);
end

if isempty(hit_idx)  % No hits found - so setting "miss" for bar plots
    norm_type =''; hmw_txt = '';
    h.current_mse_locs = nan(size(hit_lf_idx));
    h.current_mse_ori = nan(size(hit_lf_idx));
    h.current_mse_evk_waves = ones(size(hit_lf_idx,1),2)*100;
    h.inv_soln(h.current_inv_soln).half_width = repmat(struct('edges_x', nan,'val_x', nan,'edges_y', nan,'val_y', nan,'edges_z', nan,'val_z', nan,'average',nan),size(hit_lf_idx,1),1);
else        % hits found
    
    %% Location Errors
    peak_locs(hit_idx,:) = h.current_3D_peak_voxels(hit_idx,1:3);
    
    % h.current_mse_locs = nan(size(h.sim_data.cfg.source.vx_locs,1),1);
    h.current_mse_locs = sqrt( nanmean( (true_locs-peak_locs).^2 , 2));
    
    %% Orientation Errors
    true_ori = nan(size(true_locs));
    true_ori(hit_idx,:) = h.sim_data.cfg.source.vx_ori(hit_idx,:);
    peak_ori = nan(size(true_ori));
    peak_ori(hit_idx,:) = h.inv_soln(h.current_inv_soln).soln.ori(h.current_peak_hit_lf_idx(hit_idx),:);
    h.current_mse_ori = nan(size(h.sim_data.cfg.source.vx_locs,1),1);
    h.current_mse_ori(hit_idx) = rad2deg(sqrt( nanmean( (true_ori(hit_idx,:)-peak_ori(hit_idx,:)).^2 , 2)));
    h.inv_soln(h.current_inv_soln).ori_error_mse = h.current_mse_ori';
    
    %% Waveform errors - mean-squared error
    h.current_mse_evk_waves = nan(size(h.sim_data.cfg.source.vx_locs,1),1);
    act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
    ctrl_samps = h.inv_soln(h.current_inv_soln).params.ctrl_samps;
    %% finding MSE for normalized swf waves and true waves
    [pd, true_data] = sm_calc_normalized_swf(hit_lf_idx(~isnan(hit_lf_idx)));
    peak_data = zeros(size(true_data));
    peak_data(:,hit_idx) = pd;
    %% Maximum RMS of true_wave to normalize so that error is 100% for misses
    true_mse_vals_act = sqrt( nanmean( true_data(act_samps,:) .^2));
    true_mse_vals_ctrl = sqrt( nanmean( true_data(ctrl_samps,:) .^2));
    %% MSE act_samps
    if strcmp(h.inv_soln(h.current_inv_soln).maxvectorori_Type,'RMS')
        xdiff = abs(true_data) - peak_data;     % peak waveforms are already in rms so make true_data rms which is simply taking 'abs'
    else
        xdiff = true_data - peak_data;
    end
    mse_vals_act = sqrt( nanmean( xdiff(act_samps,:) .^2)) ./ true_mse_vals_act;
    mse_vals_ctrl = sqrt( nanmean( xdiff(ctrl_samps,:) .^2)) ./ true_mse_vals_ctrl;
    norm_type = ''; %'Normalized (%)';
    %% updating inv_soln
    % h.current_mse_evk_waves(hit_idx) = mse_vals;
    h.current_mse_evk_waves = [mse_vals_act; mse_vals_ctrl]'*100;
    h.inv_soln(h.current_inv_soln).wave_error_mse_norm_act = mse_vals_act*100;
    h.inv_soln(h.current_inv_soln).wave_error_mse_norm_ctrl = mse_vals_ctrl*100;
    h.inv_soln(h.current_inv_soln).wave_error_mse_norm_type = 'Wave Error (MSE) Normalized (%% max true)';
    
    %% Halfmax widths
    h.axes_invSoln_halfmax_width.clo;
    h.inv_soln(h.current_inv_soln).half_width = [];
    h.inv_soln(h.current_inv_soln).xslices =[];
    h.inv_soln(h.current_inv_soln).yslices =[];
    h.inv_soln(h.current_inv_soln).zslices =[];
    
    switch h.inv_soln(h.current_inv_soln).headmodel_type  % only calulating half-width when volume head model - can not calcualte spread function for cortically constrained
        case 'Whole Brain'
            h.axes_invSoln_halfmax_width.YLim = [0 30];
            switch h.inv_soln(h.current_inv_soln).Type
                case {'SPA' 'LCMV (FT)' 'MNE (FT)' 'sLORETA (FT)' 'eLORETA (FT)' 'SAM (FT)' 'dics (FT)' 'pcc (FT)' 'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}
                    [h.inv_soln(h.current_inv_soln).half_width,h.inv_soln(h.current_inv_soln).xslices,h.inv_soln(h.current_inv_soln).yslices,h.inv_soln(h.current_inv_soln).zslices]=bs_calc_halfmax_spread(h.inv_soln(h.current_inv_soln).soln.P.img,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,hit_lf_idx);
                    hmw_txt = '';
                case {'SIA', 'MIA' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}
                    nan_idx = nan(size(hit_lf_idx));    % setting these to nan to be able to still plot with nulls
                    [h.inv_soln(h.current_inv_soln).half_width,h.inv_soln(h.current_inv_soln).xslices,h.inv_soln(h.current_inv_soln).yslices,h.inv_soln(h.current_inv_soln).zslices]=bs_calc_halfmax_spread(h.inv_soln(h.current_inv_soln).soln.P.img,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,nan_idx);
                    hmw_txt = sprintf('Not Applicable\nfor %s',h.inv_soln(h.current_inv_soln).Type);
            end
            
        case 'Cortical Surface'
            nan_idx = nan(size(hit_lf_idx));
            [h.inv_soln(h.current_inv_soln).half_width,h.inv_soln(h.current_inv_soln).xslices,h.inv_soln(h.current_inv_soln).yslices,h.inv_soln(h.current_inv_soln).zslices]=bs_calc_halfmax_spread(h.inv_soln(h.current_inv_soln).soln.P.img,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,nan_idx);
            hmw_txt = sprintf('Not Applicable\nfor %s',h.inv_soln(h.current_inv_soln).Type);
            h.axes_invSoln_halfmax_width.YLim = [0 30];
    end
end
h.inv_soln(h.current_inv_soln).classifier_metrics.half_width = h.inv_soln(h.current_inv_soln).half_width;

%% %%%% PLOT bar graphs

if h.run_inv_soln_flag == 0
    %% Location error
    b1=bar(h.axes_invSoln_errors_locs,h.current_mse_locs,'FaceColor','Flat','CData',h.sim_data.cfg.source.src_clr); title(h.axes_invSoln_errors_locs,sprintf('Location\nError (mm)'));
    xtips2 = double(b1.XEndPoints); ytips2 = double(b1.YEndPoints);
    labels2 = string(round(b1.YData)); labels2(miss_idx) = 'Miss'; ytips2(miss_idx) = 0;
    t1 = text(h.axes_invSoln_errors_locs,xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
    
    %% Orientation Error
    b2=bar(h.axes_invSoln_errors_ori,h.current_mse_ori,'FaceColor','Flat','CData',h.sim_data.cfg.source.src_clr); title(h.axes_invSoln_errors_ori,sprintf('Orientation\nError (deg)'));
    xtips2 = double(b2.XEndPoints); ytips2 = double(b2.YEndPoints);
    labels2 = string(round(b2.YData));  labels2(miss_idx) = 'Miss'; ytips2(miss_idx) = 0;
    t1 = text(h.axes_invSoln_errors_ori,xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
    
    %% Wave Error MSE (% max(True swf))
    % plotting bar graph
    b3=bar(h.axes_invSoln_errors_waves,h.current_mse_evk_waves,'FaceColor','Flat','CData',h.sim_data.cfg.source.src_clr); title(h.axes_invSoln_errors_waves,sprintf('Best Wave Error\n(MSE)'));
    h.axes_invSoln_errors_waves.YLabel.String = norm_type;
    xtips2 = double([b3(:).XEndPoints]); ytips2 = double([b3(:).YEndPoints]);
    labels2 = string(round([b3(:).YData])); %labels2(miss_idx) = 'Miss'; %ytips2(miss_idx) = 0;
    t1 = text(h.axes_invSoln_errors_waves,xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
    h.axes_invSoln_errors_waves.XTick = [.75 1.25 1.75 2.25 2.75 3.25];
    h.axes_invSoln_errors_waves.XTickLabel = {'act' 'ctrl' 'act' 'ctrl' 'act' 'ctrl'};
    title(h.axes_invSoln_errors_waves,sprintf('Wave Error (MSE)\nNormalized (%% max true)'));
    
    %% Halfmax Width
    if ~isempty(h.inv_soln(h.current_inv_soln).half_width)
        b4=bar(h.axes_invSoln_halfmax_width, [h.inv_soln(h.current_inv_soln).half_width(1:length(h.sim_data.cfg.source.vx_idx)).average] ,'FaceColor','Flat','CData',h.sim_data.cfg.source.src_clr);
        title(h.axes_invSoln_halfmax_width,sprintf('Half-Max\n Width (mm)')); h.axes_invSoln_halfmax_width.XLabel.String = 'Peaks';
        mx=30; if max([h.inv_soln(h.current_inv_soln).half_width(:).average])*1.15>mx; mx = max([h.inv_soln(h.current_inv_soln).half_width(:).average])*1.15; end
        h.axes_invSoln_halfmax_width.YLim = [0 mx];
        xtips2 = double(b4.XEndPoints); ytips2 = double(b4.YEndPoints);
        labels2 = string(round(b4.YData)); labels2(miss_idx) = 'Miss'; ytips2(miss_idx) = 0;
        t1 = text(h.axes_invSoln_halfmax_width,xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
        t2 = text(h.axes_invSoln_halfmax_width,0.1,15,hmw_txt);
        mx=30; if max(max([h.inv_soln(h.current_inv_soln).half_width(1:length(h.sim_data.cfg.source.vx_idx)).average]))*1.15>mx; mx = max(max([h.inv_soln(h.current_inv_soln).half_width(1:length(h.sim_data.cfg.source.vx_idx)).average]))*1.15; end
        h.axes_invSoln_halfmax_width.YLim = [0 mx];
    else
        
    end

    
    %% Ylims for all plots
    mx=30; if max(h.current_mse_locs)*1.15>mx; mx = max(h.current_mse_locs)*1.15; end
    h.axes_invSoln_errors_locs.XLabel.String = 'Peaks'; h.axes_invSoln_errors_locs.YLim = [0 mx];
    mx=30; if max(h.current_mse_ori)*1.15>mx; mx = max(h.current_mse_ori)*1.15; end
    h.axes_invSoln_errors_ori.XLabel.String = 'Peaks';  h.axes_invSoln_errors_ori.YLim = [0 mx];
    
    mx=100; if max(max(h.current_mse_evk_waves))*1.15>mx; mx = max(max(h.current_mse_evk_waves))*1.15; end
    h.axes_invSoln_errors_waves.XLabel.String = 'Peaks'; h.axes_invSoln_errors_waves.YLim = [0 mx];
    
end

