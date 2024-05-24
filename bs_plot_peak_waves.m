function bs_plot_peak_waves(varargin)
global h

xdata = [];
% axes(h.axes_source_waves); cla; hold on; legend off;
h.axes_source_waves.NextPlot = 'replace';
if isfield(h, 'current_inv_swf_plots')
    h = rmfield(h, 'current_inv_swf_plots');
end

bs = round( (h.sim_data.cfg.study.base_int-h.sim_data.cfg.study.lat_sim(1))*h.sim_data.cfg.study.srate); bs(bs==0)=1;
h.sim_data.cfg.study.base_samps = bs(1):bs(2);
avg_data = nan(size(h.sim_data.sig_final,1),size(h.sim_data.sig_final,2));

if ~isfield(h.sim_data.cfg.study,'bl_bmf')
    h.sim_data.cfg.study.bl_bmf=[];
    if ~isfield(h.sim_data.cfg.study.bl_bmf,'inside_idx')
        h.sim_data.cfg.study.bl_bmf.inside_idx = find(h.inv_soln(h.current_inv_soln).leadfield.inside==1);
    end
end

% Plot active interval as a box

%% Plot TRUE source waves
if h.radio_3D_true_locs.Value==1
    h.norm_true_swf=[];
    for v=1:length(h.cfg.source.vx_idx)
        data = squeeze(h.sim_data.sig_final(:,v,:)); data = bsxfun(@minus,data,nanmean(data(h.sim_data.cfg.study.base_samps,:,:)))*h.sim_data.cfg.source.vx_amp(v);
        avg_data(:,v) = squeeze(nanmean(data,2));
    end
    if h.radio_normalize_swf.Value == 1
        %         mu = repmat(nanmean(avg_data(h.sim_data.cfg.study.base_samps,v)),size(avg_data,1),1);
        %         sigma = repmat(nanstd(avg_data(h.sim_data.cfg.study.base_samps,v)),size(avg_data,1),1);
        %
        %% New normalization procedure --> Normalize to +/1 max(abs) within act_samps
        act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
        max_p = nanmax(nanmax(abs(avg_data(act_samps,:))));
        avg_data = 100*(avg_data ./ max_p);
    end
    h.axes_source_waves.NextPlot = 'replace';
    for v=1:length(h.cfg.source.vx_idx)
%         if h.radio_normalize_swf.Value == 1 % plot normalized waves
%             %             h.norm_true_swf(:,v) = (avg_data(:,v)-mu) ./ sigma;
            h.norm_true_swf(:,v) = avg_data(:,v);
%         else
%             h.norm_true_swf(:,v) = avg_data(:,v);
%         end
        
        switch h.inv_soln(h.current_inv_soln).Type
            case {'SPA' 'SIA' 'MIA' 'LCMV (FT)' 'sLORETA (FT)'  'dics (FT)' 'pcc (FT)' 'SAM (FT)' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}     % scalar
            case {'eLORETA (FT)' 'MNE (FT)' 'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}   % vector
                switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
                    case {'RMS' 'rms'}
                        h.norm_true_swf(:,v) = abs(h.norm_true_swf(:,v));     % taking abs to match calculated swf
                    case {'max' 'Max'}    % vector dipole orientation with maximum swf power between active and control interval
                        h.norm_true_swf(:,v) = h.norm_true_swf(:,v);     % taking abs to match calculated swf
                    case 'avg.pow'
                        h.norm_true_swf(:,v) = abs(h.norm_true_swf(:,v));     % taking abs to match calculated swf
                end
        end
        h.norm_true_swf = bsxfun(@minus,h.norm_true_swf, nanmean(h.norm_true_swf(h.cfg.study.base_samps,:)));
        
        
        h.current_true_swf_plots(v) = plot(h.axes_source_waves,h.sim_data.cfg.study.lat_sim, h.norm_true_swf(:,v),'--','color',h.cfg.source.src_clr(v,:),'linewidth',1);
        h.axes_source_waves.NextPlot = 'add';

    end
end

%% Plot Peak Locations
% try
% if h.radio_3D_peak_locs.Value==1
% calculating source waveforms 'swf' for peaks in the 3D maps
p_idx = ~isnan(h.current_3D_peak_idx);
h.current_peak_swf = nan(size(h.sim_data.sens_final,1),length(p_idx));

switch h.inv_soln(h.current_inv_soln).Type
    case {'Dipole' 'SPA' 'LCMV (FT)' 'sLORETA (FT)'  'dics (FT)' 'pcc (FT)'  'SAM (FT)' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}    % BRANE Lab beamformers
        h.current_peak_swf(:,p_idx) = [h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx(p_idx))'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';
    case {'SIA' 'MIA'}    % BRANE Lab beamformers
        %             h.current_peak_swf = [h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx)'*squeeze(nanmean(h.sim_data.sens_final,3))']';
        %         h.current_peak_swf = [h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx)'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';
        h.current_peak_swf(:,p_idx) = [h.inv_soln(h.current_inv_soln).soln.residual_wts(:,h.current_3D_peak_idx(p_idx))'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';
        %         h.current_peak_swf = [h.inv_soln(h.current_inv_soln).soln.nulled_wts(:,h.current_3D_peak_idx)'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';
        
    case {'eLORETA (FT)' 'MNE (FT)'}    % Field Trips Vector inverse solutions
        act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
        ctrl_samps = h.inv_soln(h.current_inv_soln).params.ctrl_samps;
        clear swf swf_pwr
        for ox=1:3
            swf(:,ox,:)=squeeze(nanmean(h.sim_data.sens_final,3))*squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx(p_idx),ox));
            swf_pwr(ox,:)=rms(swf(act_samps,ox,:),1)./rms(swf(ctrl_samps,ox,:),1);
        end
        %        h.current_peak_swf=[];
        switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
            case 'RMS'
                h.current_peak_swf(:,p_idx) = squeeze(rms(swf,2));
            case 'Max'    % vector dipole orientation with maximum swf power between active and control interval
                [~,max_ori]=max(swf_pwr);   % maximum orientation
                for v=1:size(swf,3); h.current_peak_swf(:,v) = squeeze(swf(:,max_ori(v),v)); end
            case 'avg.pow'
                switch h.inv_soln(h.current_inv_soln).Type
                    case 'eLORETA (FT)'
                        h.current_peak_swf(:,find(p_idx==1)) = squeeze(nanmean(abs(swf),2));
                    case 'MNE (FT)'
                        switch h.inv_soln(h.current_inv_soln).headmodel_type
                            case 'Whole Brain'
                                h.current_peak_swf(:,find(p_idx==1)) = h.inv_soln(h.current_inv_soln).soln.avg.pow(h.inv_soln(h.current_inv_soln).params.inside_idx(find(p_idx==1)),:)';
                            case 'Cortical Surface'
                                h.current_peak_swf(:,find(p_idx==1)) = h.inv_soln(h.current_inv_soln).soln.avg.pow(h.current_3D_peak_idx,:)';
                        end
                end
                
        end
    case {'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}    % Brainstorm's vector inverse solutions
        swf = h.inv_soln(h.current_inv_soln).soln.ImagingKernel * squeeze(nanmean(h.sim_data.sens_final,3))';
        switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
            case 'RMS'
                iVertSource = 1:size(swf,1);
                idx1 = 1:3:length(iVertSource); idx2 = 2:3:length(iVertSource); idx3 = 3:3:length(iVertSource);
                swf_rms = squeeze(sqrt(swf(idx1,:).^2 + swf(idx2,:).^2 + swf(idx3,:).^2));
                h.current_peak_swf(:,p_idx) = swf_rms(h.current_3D_peak_idx(p_idx),:)';
            case 'Max'    % selecting waveform for maximum ori "max_ori" vector already calculatd for the inv_soln
                dims = size(swf);
                swf = permute(reshape(swf,[3 dims(1)/3 dims(2)]),[3 2 1]); % [3dipole_vecotrs x grid locs]
                vx_idx = h.current_3D_peak_idx(p_idx);
                swf2 = [];
                for v=1:length(vx_idx); swf2(:,v) = swf(:,vx_idx(v),h.inv_soln(h.current_inv_soln).soln.max_ori(vx_idx(v))); end
                %                 iVertSource = 1:size(swf,1);
                %                 idx1 = 1:3:length(iVertSource); idx2 = 2:3:length(iVertSource); idx3 = 3:3:length(iVertSource);
                %                 [swf_max, max_ori] = max(abs(squeeze(cat(3, swf(idx1,:), swf(idx2,:), swf(idx3,:)))), [], 3);
                %                 h.current_peak_swf(:,p_idx) = swf_max(h.current_3D_peak_idx(p_idx),:)';
                h.current_peak_swf(:,p_idx) = swf2;
        end
end


% baselining
%     h.current_peak_swf = bsxfun(@minus,h.current_peak_swf,nanmean(h.current_peak_swf(h.sim_data.cfg.study.base_samps,:)));
h.current_peak_swf = bsxfun(@minus,h.current_peak_swf,nanmean(h.current_peak_swf(h.cfg.study.base_samps,:)));
% Flip wave --> flipping wave if orientation is in oppposite quadrant (i.e., close to 180deg out of phase with true orientation);

%     try
if any(~isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.Hits))   %~isempty(h.current_3D_peak_voxels)
    hit_idx = h.inv_soln(h.current_inv_soln).classifier_metrics.Hits;
    for v=1:length(hit_idx)    % only flipping the peaks associated with true sources that have been reordered from 1:true_sources
        lf_idx = h.inv_soln(h.current_inv_soln).classifier_metrics.Hits(v);
        if isnan(lf_idx)   % Miss
        else    % Hit
            true_ori = h.sim_data.cfg.source.vx_ori(v,:);
            peak_ori = h.inv_soln(h.current_inv_soln).soln.ori(lf_idx,:);
            % negative correlations means that the source wave needs to be flipped
            xr = corr(h.current_peak_swf(h.cfg.study.bl_bmf.act_samps,v),h.norm_true_swf(h.cfg.study.bl_bmf.act_samps,v));
            
            if xr<0 % nansum(abs(true_ori-peak_ori)) > nansum(abs(true_ori+peak_ori)) % flip swf
                h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(v),:) = h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(v),:)*-1;
                
                % flip orientation if the asbolute difference is larger that the sum
                peak_ori = nan(size(true_ori));
                peak_ori2 = h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(v),:);
                peak_ori(1:size(peak_ori2,1),1:size(peak_ori2,2)) = peak_ori2;
                
                if nansum(abs(true_ori-peak_ori)) > nansum(abs(true_ori+peak_ori))
                    h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(v),:) = h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(v),:)*-1;
                end
                
                switch h.inv_soln(h.current_inv_soln).Type
                    case {'SPA' 'SIA' 'MIA' 'LCMV (FT)' 'sLORETA (FT)'  'dics (FT)' 'pcc (FT)' 'SAM (FT)' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}     % scalar
                        h.current_peak_swf(:,v) = h.current_peak_swf(:,v)*-1;
                        h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx(v)) = h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx(v))*-1;
                    case {'eLORETA (FT)' 'MNE (FT)' 'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'  }   % vector
                        switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
                            case {'max' 'Max'}   % allow flipping because waveform is not absolute
                                h.current_peak_swf(:,v) = h.current_peak_swf(:,v)*-1;
                                h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx(v)) = h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx(v))*-1;
                        end
                end
                
                %         vx_pos = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx(v),:);
                %         clear ori_pos
                %         hv = handle(h.map3D_peak_ori(v));
                %         ori_pos(1,:) = hv.XData-vx_pos(1);   % center ori at vox_pos then flip it
                %         ori_pos(2,:) = hv.YData-vx_pos(2);
                %         ori_pos(3,:) = hv.ZData-vx_pos(3);
                %
                %         hv.XData = vx_pos(1) + ori_pos(1,:)*-1;   % flipping orientation on image
                %         hv.YData = vx_pos(2) + ori_pos(2,:)*-1;   % flipping orientation on image
                %         hv.ZData = vx_pos(3) + ori_pos(3,:)*-1;   % flipping orientation on image
                %
            end
        end
    end
    %     catch
    %         fprintf('ERROR in automatically flipping orientations for %s\n', h.inv_soln(h.current_inv_soln).Type)
    %     end
    
    if h.radio_normalize_swf.Value == 1 % plot normalized waves
        %         %% normalizing z-transform relative to standard deviation of the baseline
        %         mu_base = repmat(nanmean(h.current_peak_swf(h.sim_data.cfg.study.base_samps,:)),size(h.current_peak_swf,1),1);
        %         sigma_base = repmat(nanstd(h.current_peak_swf(h.sim_data.cfg.study.base_samps,:)),size(h.current_peak_swf,1),1);
        %         norm_swf = (h.current_peak_swf-mu_base) ./ sigma_base;    % Z-score
        
        %% New normalization procedure --> Normalize to +/1 max(abs) within act_samps
        if ~isempty(find(p_idx==1, 1))
            avg_data = h.current_peak_swf(:,p_idx);   % getting data just for found peaks
            act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
            max_p = nanmax(nanmax(abs(avg_data(act_samps,:))));
            norm_swf = 100* (h.current_peak_swf ./ max_p);
            h.current_norm_peak_swf = norm_swf;  % peak waves are all z-score normalized now to baseline
            h.current_norm_peak_swf = bsxfun(@minus,h.current_norm_peak_swf,nanmean(h.current_norm_peak_swf(h.cfg.study.base_samps,:)));
        else
             avg_data = h.current_peak_swf;   % getting data for any found peaks
            act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
            max_p = nanmax(nanmax(abs(avg_data(act_samps,:))));
            norm_swf = 100* (h.current_peak_swf ./ max_p);
            h.current_norm_peak_swf = norm_swf;  % peak waves are all z-score normalized now to baseline
            h.current_norm_peak_swf = bsxfun(@minus,h.current_norm_peak_swf,nanmean(h.current_norm_peak_swf(h.cfg.study.base_samps,:)));
        end
        
         %     h.axes_source_waves.YLabel.String = 'Normalized to Baseline';
        h.axes_source_waves.YLabel.String = 'Normalized (Max Hits)';
    else
        h.current_norm_peak_swf = h.current_peak_swf;  % peak waves are all z-score normalized now to baseline
        h.current_norm_peak_swf = bsxfun(@minus,h.current_norm_peak_swf,nanmean(h.current_norm_peak_swf(h.cfg.study.base_samps,:)));
        h.axes_source_waves.YLabel.String = 'Inv Solution Output';
    end
    
    
    
    if h.radio_3D_peak_locs.Value==1
        %% plotting normalized inv peaks waves
%         h.current_inv_swf_plots = [];
       h.ln_clr = bsxfun(@mtimes,ones(length(h.current_3D_peak_idx),3),h.FA_clr);  %lines(length(h4.current_3D_peak_idx));
        h.ln_clr(1:size(h.cfg.source.src_clr,1),:) = h.cfg.source.src_clr;
  
        for v=1:size(h.current_peak_swf,2)
            h.current_inv_swf_plots(v) = plot(h.axes_source_waves,h.sim_data.cfg.study.lat_sim,h.current_norm_peak_swf(:,v),'color',h.ln_clr(v,:),'linewidth',2);
            h.current_inv_swf_plots(v).Color(4) = h.false_positive_lineAlpha;
            h.axes_source_waves.NextPlot = 'add';
        end
        
        % else
        % end
    end
%     axes(h.axes_source_waves); legend off;
%     try
        if h.radio_3D_true_locs.Value==1 && h.radio_3D_peak_locs.Value==0 && ~isempty(h.current_3D_peak_idx) % only true source waves
            xdata = h.norm_true_swf;
            lg_txt=[]; for v=1:size(xdata,2); lg_txt{v} = sprintf('True Source %.f',v); end
            legend(h.axes_source_waves,h.current_true_swf_plots(1:size(xdata,2)),lg_txt);
        elseif h.radio_3D_true_locs.Value==0 && h.radio_3D_peak_locs.Value==1  && ~isempty(h.current_3D_peak_idx) % only peak waves
            xdata = h.current_norm_peak_swf;
            sel_v = h.listbox_peaks_found.Value;
            lg_txt = num2str(h.inv_soln(h.current_inv_soln).peak_idx(sel_v)');
            legend(h.axes_source_waves, h.current_inv_swf_plots(sel_v),lg_txt);
        elseif h.radio_3D_true_locs.Value==1 && h.radio_3D_peak_locs.Value==1  && ~isempty(h.current_3D_peak_idx) % both
            xdata = cat(2,h.current_norm_peak_swf,h.norm_true_swf);
            xdata = h.norm_true_swf;
            lg_txt1=[]; for v=1:size(xdata,2); lg_txt1{v} = sprintf('True Source %.f',v); end
            sel_v = h.listbox_peaks_found.Value;
%             lg_txt2 = cellstr(num2str(h.inv_soln(h.current_inv_soln).peak_idx(sel_v)'))';
            lg_txt2 = {};
            for v=1:length(sel_v)
                lg_txt2(v) = { num2str(h.inv_soln(h.current_inv_soln).peak_idx(sel_v(v)))};
            end
%             legend(h.axes_source_waves,[h.current_true_swf_plots(1:size(xdata,2)) h.current_inv_swf_plots(sel_v)],[lg_txt1 lg_txt2]);
            legend(h.axes_source_waves,[h.current_true_swf_plots(1:size(xdata,2)) h.current_inv_swf_plots(sel_v)],[lg_txt1(1:size(xdata,2)) lg_txt2{1:length(sel_v)}]);
        else
            %      xdata = h.norm_true_swf;
            %    lg_txt=[]; for v=1:size(xdata,2); lg_txt{v} = sprintf('True Source %.f',v); end
            %     legend(h.axes_source_waves,h.current_true_swf_plots(1:size(xdata,2)),lg_txt);
        end
        
        h.axes_source_waves.Legend.Position = [0.4875    0.2826    0.1350    0.1233];
        h.axes_source_waves.Legend.FontSize = 6.5;
        min_max = [-1.1*max(max(abs(xdata))) 1.1*max(max(abs(xdata))) ]; 
        if (~isempty(xdata) && any(min_max~=0)) && all(~isnan(min_max))
            h.axes_source_waves.YLim = min_max;
        else
           try; set(h.current_inv_swf_plots,'Visible','off'); end
             try; set(h.current_true_swf_plots,'Visible','off'); end
        end
        
        
%     catch
%     end
    h.axes_source_waves.XLabel.String = 'Time (sec)';
    %     h.axes_source_waves.YLabel.String = 'Normalized Amplitude';
    h.axes_source_waves.Visible = true;
    h.axes_source_waves.Title.String = 'Peak Waves for Inverse Solution';
    
    %% plot red time line
    if h.btn_3D_plot_peak_waves.Value==1
        if isfield(h,'current_swf_time_plot')
            if isvalid(h.current_swf_time_plot)
                h.current_swf_time_plot.XData = [h.current_swf_time(1) h.current_swf_time(1)];
                h.current_swf_time_plot.YData = h.axes_source_waves.YLim;
            else
                h.current_swf_time_plot = plot(h.axes_source_waves,[h.current_swf_time(1) h.current_swf_time(1)], h.axes_source_waves.YLim,'r--');
            end
        else
            h.current_swf_time_plot = plot(h.axes_source_waves,[h.current_swf_time(1) h.current_swf_time(1)], h.axes_source_waves.YLim,'r--');
        end
    end
    
    
    if ~isempty(h.current_peak_swf)
        bs_calc_fft;
    end
    
    % catch me
    %     fprintf('ERROR: Source weights do no match number of sensors\nLikely using weights calculated using EEG but current sensors are MEG, or vice versa\n\n');
    %     h.current_peak_swf = [];
    %     h.axes_source_fft.clo;
    %     h.axes_source_waves.clo;
    %     h.axes_invSoln_errors_locs.clo;
    %     h.axes_invSoln_errors_ori.clo;
    %     h.axes_invSoln_errors_waves.clo;
    % end
    disableDefaultInteractivity(h.axes_source_waves); h.axes_source_waves.Toolbar.Visible='off';
    %     h.axes_source_waves.XLim = str2num(h.edit_plot_time_int.String);
end
bs_calc_errors_inv_soln;
sm_set_source_wave_scales('Both');  % updating X & Y scales
sm_calc_localizer_performance;

%% reinitialize spatiotemp
    if h.btn_3D_plot_peak_waves.Value == 1
        h.axes_source_waves.ButtonDownFcn = @bs_plot_map_time;
    else
        h.axes_source_waves.ButtonDownFcn = [];
    end
    axis(h.axes_3D_images,'off');

% %% hiding peak waves below threshold
% if ~isempty(h.current_inv_soln_hide_peak_idx)   % Hiding peaks and slices - need to do first because some peaks share slices
%     sm_show_peaks(h.current_inv_soln_hide_peak_idx,'off');  % shows ('on') or hides ('off') selected peaks in 3D image
% end

