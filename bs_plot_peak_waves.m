function bs_plot_peak_waves(varargin)
global h
axes(h.axes_source_waves); cla; hold on; legend off;
try
    if ~isempty(h.map3D_peak_locs(1).CData)
        ln_clr =  h.map3D_peak_locs(1).CData;
    end
catch
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

%% Plot TRUE source locations
if h.radio_3D_true_locs.Value==1
    h.norm_true_swf=[];
    for v=1:3
        data = squeeze(h.sim_data.sig_final(:,v,:)); data = bsxfun(@minus,data,nanmean(data(h.sim_data.cfg.study.base_samps,:,:)))*str2num(h.edit_source_amp(v).String);
        avg_data(:,v) = squeeze(nanmean(data,2));
        mu = repmat(nanmean(avg_data(h.sim_data.cfg.study.base_samps,v)),size(avg_data,1),1);
        sigma = repmat(nanstd(avg_data(h.sim_data.cfg.study.base_samps,v)),size(avg_data,1),1);
        if h.radio_normalize_swf.Value == 1 % plot normalized waves
            h.norm_true_swf(:,v) = (avg_data(:,v)-mu) ./ sigma;
        else
            h.norm_true_swf(:,v) = avg_data(:,v);
        end
        
        switch h.inv_soln(h.current_inv_soln).Type
            case {'SPA' 'SIA' 'MIA' 'LCMV'  'dics' 'pcc' 'sLORETA' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}     % scalar
            case {'eLORETA' 'MNE'}   % vector
                switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
                    case 'RMS'
                        h.norm_true_swf(:,v) = abs(h.norm_true_swf(:,v));     % taking abs to match calculated swf
                    case 'Max'    % vector dipole orientation with maximum swf power between active and control interval
                        h.norm_true_swf(:,v) = h.norm_true_swf(:,v);     % taking abs to match calculated swf
                    case 'avg.pow'
                        h.norm_true_swf(:,v) = abs(h.norm_true_swf(:,v));     % taking abs to match calculated swf
                end
        end
        h.norm_true_swf = bsxfun(@minus,h.norm_true_swf,nanmean(h.norm_true_swf(h.cfg.study.base_samps,:)));

        
        h.current_true_swf_plots(v) = plot(h.axes_source_waves,h.sim_data.cfg.study.lat_sim, h.norm_true_swf(:,v),'-','color',h.src_clr(v,:),'linewidth',1);
    end
end

%% Plot Peak Locations
% try
% if h.radio_3D_peak_locs.Value==1
% calculating source waveforms 'swf' for peaks in the 3D maps
switch h.inv_soln(h.current_inv_soln).Type
    case {'SPA' 'LCMV' 'sLORETA'  'dics' 'pcc' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}    % BRANE Lab beamformers
        h.current_peak_swf = [h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx)'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';
    case {'SIA' 'MIA'}    % BRANE Lab beamformers
        %             h.current_peak_swf = [h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx)'*squeeze(nanmean(h.sim_data.sens_final,3))']';
%         h.current_peak_swf = [h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx)'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';
        h.current_peak_swf = [h.inv_soln(h.current_inv_soln).soln.residual_wts(:,h.current_3D_peak_idx)'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';
%         h.current_peak_swf = [h.inv_soln(h.current_inv_soln).soln.nulled_wts(:,h.current_3D_peak_idx)'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';

    case {'eLORETA' 'MNE'}    % Field Trips Vector inverse solutions
       act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
       ctrl_samps = h.inv_soln(h.current_inv_soln).params.ctrl_samps;
       clear swf swf_pwr
       for ox=1:3
           swf(:,ox,:)=squeeze(nanmean(h.sim_data.sens_final,3))*squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx,ox));
           swf_pwr(ox,:)=rms(swf(act_samps,ox,:),1)./rms(swf(ctrl_samps,ox,:),1);
       end
       h.current_peak_swf=[];
       switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
           case 'RMS'
               h.current_peak_swf = squeeze(rms(swf,2));
           case 'Max'    % vector dipole orientation with maximum swf power between active and control interval
               [~,max_ori]=max(swf_pwr);   % maximum orientation
               for v=1:size(swf,3); h.current_peak_swf(:,v) = squeeze(swf(:,max_ori(v),v)); end
           case 'avg.pow'
               switch h.inv_soln(h.current_inv_soln).Type
                   case 'eLORETA'
                       h.current_peak_swf = squeeze(nanmean(abs(swf),2));
                   case 'MNE'
                       switch h.inv_soln(h.current_inv_soln).headmodel_type
                           case 'Whole Brain'
                               h.current_peak_swf = h.inv_soln(h.current_inv_soln).soln.avg.pow(h.inv_soln(h.current_inv_soln).params.inside_idx(h.current_3D_peak_idx),:)';
                           case 'Cortical Surface'
                               h.current_peak_swf = h.inv_soln(h.current_inv_soln).soln.avg.pow(h.current_3D_peak_idx,:)';
                       end
               end
               
       end
        
end


% baselining
%     h.current_peak_swf = bsxfun(@minus,h.current_peak_swf,nanmean(h.current_peak_swf(h.sim_data.cfg.study.base_samps,:)));
h.current_peak_swf = bsxfun(@minus,h.current_peak_swf,nanmean(h.current_peak_swf(h.cfg.study.base_samps,:)));
% Flip wave --> flipping wave if orientation is in oppposite quadrant (i.e., close to 180deg out of phase with true orientation);

%     try
if ~isempty(h.current_3D_peak_voxels)
    
[v_idx]=find_nearest_voxel(h.current_3D_peak_voxels(:,1:3),h.cfg.source.vx_locs);
vx_idx = unique(v_idx);
for vvx=1:length(vx_idx)    % only flipping the first 3 that have been reordered
    v = vx_idx(vvx);
    true_ori = h.sim_data.cfg.source.vx_ori(v,:);
    peak_ori = nan(size(true_ori));
    peak_ori2 = h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(vvx),:);
    peak_ori(1:size(peak_ori2,1),1:size(peak_ori2,2)) = peak_ori2;
    %         peak_ori = h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(v),:);

    % negative correlations means that the source wave needs to be flipped
    xr = corr(h.current_peak_swf(h.cfg.study.bl_bmf.act_samps,vvx),h.norm_true_swf(h.cfg.study.bl_bmf.act_samps,v));
    
    if xr<0 % nansum(abs(true_ori-peak_ori)) > nansum(abs(true_ori+peak_ori)) % flip swf
        h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(vvx),:) = h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(vvx),:)*-1;
        
        % flip orientation if the asbolute difference is larger that the sum
        
        peak_ori = nan(size(true_ori));
        peak_ori2 = h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(vvx),:);
        peak_ori(1:size(peak_ori2,1),1:size(peak_ori2,2)) = peak_ori2;
        if nansum(abs(true_ori-peak_ori)) > nansum(abs(true_ori+peak_ori))
            h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(vvx),:) = h.inv_soln(h.current_inv_soln).soln.ori(h.current_3D_peak_idx(vvx),:)*-1;
        end
                
        switch h.inv_soln(h.current_inv_soln).Type
            case {'SPA' 'SIA' 'MIA' 'LCMV' 'sLORETA'  'dics' 'pcc' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}     % scalar
                h.current_peak_swf(:,v) = h.current_peak_swf(:,v)*-1;
                h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx(vvx)) = h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx(vvx))*-1;
            case {'eLORETA' 'MNE'}   % vector
                if strcmp(h.inv_soln(h.current_inv_soln).maxvectorori_Type,'Max')   % allow flipping because waveform is not absolute
                    h.current_peak_swf(:,v) = h.current_peak_swf(:,v)*-1;
                    h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx(vvx)) = h.inv_soln(h.current_inv_soln).soln.wts(:,h.current_3D_peak_idx(vvx))*-1;
                end
                
        end
        
        vx_pos = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(h.current_3D_peak_idx(vvx),:);
        clear ori_pos
        hv = handle(h.map3D_peak_ori(v));
        ori_pos(1,:) = hv.XData-vx_pos(1);   % center ori at vox_pos then flip it
        ori_pos(2,:) = hv.YData-vx_pos(2);
        ori_pos(3,:) = hv.ZData-vx_pos(3);
        
        hv.XData = vx_pos(1) + ori_pos(1,:)*-1;   % flipping orientation on image
        hv.YData = vx_pos(2) + ori_pos(2,:)*-1;   % flipping orientation on image
        hv.ZData = vx_pos(3) + ori_pos(3,:)*-1;   % flipping orientation on image
       
    end
end
%     catch
%         fprintf('ERROR in automatically flipping orientations for %s\n', h.inv_soln(h.current_inv_soln).Type)
%     end

if h.radio_normalize_swf.Value == 1 % plot normalized waves
    %% normalizing z-transform relative to standard deviation of the baseline
    mu_base = repmat(nanmean(h.current_peak_swf(h.sim_data.cfg.study.base_samps,:)),size(h.current_peak_swf,1),1);
    sigma_base = repmat(nanstd(h.current_peak_swf(h.sim_data.cfg.study.base_samps,:)),size(h.current_peak_swf,1),1);
    
    norm_swf = (h.current_peak_swf-mu_base) ./ sigma_base;    % Z-score
    h.current_norm_peak_swf = norm_swf;  % peak waves are all z-score normalized now to baseline
    h.axes_source_waves.YLabel.String = 'Normalized to Baseline';
else
    h.current_norm_peak_swf = h.current_peak_swf;  % peak waves are all z-score normalized now to baseline
    h.axes_source_waves.YLabel.String = 'Inv Solution Output';
end
h.current_norm_peak_swf = bsxfun(@minus,h.current_norm_peak_swf,nanmean(h.current_norm_peak_swf(h.cfg.study.base_samps,:)));

if h.radio_3D_peak_locs.Value==1
    %% plotting normalized inv peaks waves
    for v=fliplr(1:size(h.current_peak_swf,2))
        h.current_inv_swf_plots(v) = plot(h.axes_source_waves,h.sim_data.cfg.study.lat_sim,h.current_norm_peak_swf(:,v),'color',ln_clr(v,:),'linewidth',2);
        %         if v>3
        h.current_inv_swf_plots(v).Color(4) = h.false_positive_lineAlpha;
        %         end
        
    end
    
    % else
    % end
end
axes(h.axes_source_waves); legend off;
if h.radio_3D_true_locs.Value==1 && h.radio_3D_peak_locs.Value==0 && ~isempty(h.current_3D_peak_idx) % only true source waves
    legend(h.axes_source_waves,h.current_true_swf_plots,{'True Source 1' 'True Source 2' 'True Source 3'});
    xdata = h.norm_true_swf;
elseif h.radio_3D_true_locs.Value==0 && h.radio_3D_peak_locs.Value==1  && ~isempty(h.current_3D_peak_idx) % only peak waves
    legend(h.axes_source_waves,h.current_inv_swf_plots(1),{'Peak Waves'});
    xdata = h.current_norm_peak_swf;
elseif h.radio_3D_true_locs.Value==1 && h.radio_3D_peak_locs.Value==1  && ~isempty(h.current_3D_peak_idx) % both
    legend(h.axes_source_waves,{'True Source 1' 'True Source 2' 'True Source 3' 'Peak Waves'});
    xdata = cat(2,h.current_norm_peak_swf,h.norm_true_swf);
else
    legend(h.axes_source_waves,{'True Source 1' 'True Source 2' 'True Source 3'});
    xdata = h.norm_true_swf;
end
h.axes_source_waves.Legend.Position = [0.51    0.2795    0.135    0.075];
h.axes_source_waves.Legend.FontSize = 6.5;

h.axes_source_waves.YLim = [-1.1*max(max(abs(xdata))) 1.1*max(max(abs(xdata))) ];
h.axes_source_waves.XLim = str2num(h.edit_plot_time_int.String);
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
bs_calc_errors_inv_soln;

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
h.axes_source_waves.XLim = str2num(h.edit_plot_time_int.String);
end
sm_calc_localizer_performance;


