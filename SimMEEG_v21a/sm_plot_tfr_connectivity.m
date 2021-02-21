function sm_plot_tfr_connectivity(varargin)
global h

if h.menu_inv_tfr_type.Value<=3 && h.menu_inv_tfr_data_type.Value>=3; h.menu_inv_tfr_data_type.Value=1; end % Power data currently only have "Data" not "Noise" or "Surrogate" data yet
h.FC_alpha_level = str2num(h.edit_3D_FC_stats_alpha_level.String);

% turning things off
disableDefaultInteractivity(h.axes_inv_soln_tfr); h.axes_inv_soln_tfr.Toolbar.Visible='off'; pause(.1);
h.axes_inv_soln_tfr.Visible = 'off'; for a=1:length(h.axes_inv_soln_tfr.Children); h.axes_inv_soln_tfr.Children(a).Visible='off'; end
if isfield(h,'current_plv_plots'); if any(isvalid(h.current_plv_plots)); delete(h.current_plv_plots); end; end
if h.radio_inv_plot_peak_tfr_connectivity.Value==0 && h.radio_inv_plot_true_tfr_connectivity.Value==0
    h.axes_source_fft.Visible = 'on'; for a=1:length(h.axes_source_fft.Children); h.axes_source_fft.Children(a).Visible='on'; end
    h.menu_inv_tfr_type.Visible = 'off';
    h.menu_inv_tfr_data_type.Visible = 'off';
    h.edit_inv_tfr_caxis_txt.Visible = 'off';
    h.edit_inv_tfr_caxis.Visible = 'off';
    h.radio_inv_plot_connectivity_graph.Visible ='off';
    if isfield(h,'current_inv_tfr_point_lines'); if any(isvalid(h.current_inv_tfr_point_lines)); delete(h.current_inv_tfr_point_lines); end; end
    h.edit_inv_plv_thresh_txt.Visible = 'off'; h.edit_inv_plv_thresh.Visible = 'off';
end

if ~isfield(h.inv_soln(h.current_inv_soln),'TFR_results')
    msgbox(sprintf('Please perform Time-Frequency Response (TFR) analyses\n\nClick on "Calc Peak Connectivity"\n'));
    h.radio_inv_plot_peak_tfr_connectivity.Value = 0;
else
    if isempty(h.inv_soln(h.current_inv_soln).TFR_results)
        msgbox(sprintf('Please perform Time-Frequency Response (TFR) analyses\n\nClick on "Calc Peak Connectivity"\n'));
        h.radio_inv_plot_peak_tfr_connectivity.Value = 0;
    else
        %% need to delete colorbar
        if isfield(h,'colorbar_axes_inv_soln_tfr')
            if isvalid(h.colorbar_axes_inv_soln_tfr)
                delete(h.colorbar_axes_inv_soln_tfr);
            end
        end
        
        if h.radio_inv_plot_peak_tfr_connectivity.Value == 0
        else
            h.axes_inv_soln_tfr.Visible = 'on'; for a=1:length(h.axes_inv_soln_tfr.Children); h.axes_inv_soln_tfr.Children(a).Visible='on'; end
            h.axes_source_fft.Visible = 'off'; for a=1:length(h.axes_source_fft.Children); h.axes_source_fft.Children(a).Visible='off'; end
            h.menu_inv_tfr_type.Visible = 'on';
            h.menu_inv_tfr_data_type.Visible = 'on';
            h.edit_inv_tfr_caxis_txt.Visible = 'on';
            h.edit_inv_tfr_caxis.Visible = 'on';
            h.radio_inv_plot_connectivity_graph.Visible ='on';
            h.edit_inv_plv_thresh_txt.Visible = 'on'; h.edit_inv_plv_thresh.Visible = 'on';
            
            %% start plotting
            % plot tfr for selected source in list box - if multiple selected then average across all
            
            v_idx = h.listbox_peaks_found.Value;
            % find if v_idx is/are in seed_idx or comp_idx
            
            seed_true = ismember(h.inv_soln(h.current_inv_soln).plv_seed_idx,h.current_3D_peak_idx(v_idx));
            comp_true = ismember(h.inv_soln(h.current_inv_soln).plv_comp_idx,h.current_3D_peak_idx(v_idx));
            
            if all(seed_true==0) && all(comp_true==0)
                h.axes_inv_soln_tfr.clo;
                text(h.axes_inv_soln_tfr,h.axes_inv_soln_tfr.XLim(1)+range(h.axes_inv_soln_tfr.XLim)/8,range(h.axes_inv_soln_tfr.YLim/2),sprintf('Time-Frequency not calculated\nfor selected source'));
                return
            elseif any(seed_true==0) && any(comp_true==0)    % removing redundant sources from comp indices
                comp_rm = ismember(h.inv_soln(h.current_inv_soln).plv_comp_idx(comp_true),h.inv_soln(h.current_inv_soln).plv_seed_idx(seed_true));
                comp_true(comp_rm)=false;
            end
            
            %% find all plv of selected peak sources
            plv_idx = h.inv_soln(h.current_inv_soln).plv_contrast_idx;
            plv_true = any(ismember(plv_idx,v_idx),2);
            h.current_plv_idx = find(plv_true==1);
            
            lat = h.sim_data.cfg.study.lat_sim;
            lat_coi = lat;
            %% TFR plot
            tfr_seed = []; tfr_comp = [];
            switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
                case 'Total Power'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'
                            tfr_seed = h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt(:,:,seed_true);
                            tfr_comp = h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt(:,:,comp_true);
                        case 'Noise'
                            tfr_seed = h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_seed_wt(:,:,seed_true);
                            tfr_comp = h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_comp_wt(:,:,comp_true);
                    end
                case 'Evoked Power'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'
                            tfr_seed = h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_evk(:,:,seed_true);
                            tfr_comp = h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_evk(:,:,comp_true);
                        case 'Noise'
                            tfr_seed = h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_seed_wt_evk(:,:,seed_true);
                            tfr_comp = h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_comp_wt_evk(:,:,comp_true);
                    end
                case 'Induced Power'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'
                            tfr_seed = h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_ind(:,:,seed_true);
                            tfr_comp = h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_ind(:,:,comp_true);
                        case 'Noise'
                            tfr_seed = h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_seed_wt_ind(:,:,seed_true);
                            tfr_comp = h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_comp_wt_ind(:,:,comp_true);
                    end
                case 'PLV'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'; tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.plv_based(:,plv_true,:),[1 3 2]);
                        case 'Noise'; tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.Noise.plv_based(:,plv_true,:),[1 3 2]);
                        case 'Surrogate'; tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.plv_surg_based_mean(:,plv_true,:),[1 3 2]);
                        case 'Data (Surrogate thresholded)'
                            tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.plv_based(:,:,:),[1 3 2]);
                            mu = permute(h.inv_soln(h.current_inv_soln).TFR_results.plv_surg_based_mean(:,:,:),[1 3 2]);
                            sigma = permute(h.inv_soln(h.current_inv_soln).TFR_results.plv_surg_based_std(:,:,:),[1 3 2]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and surrogate
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            
                            % creating significance mask
                            x_mask = nan(size(tfr_seed)); x_mask(tfr_seed>surg_thresh_pos)=1; x_mask(tfr_seed<surg_thresh_neg)=1;
                            tfr_seed = tfr_seed.*x_mask;
                            tfr_seed(isnan(tfr_seed))=nan;
                            h.current_inv_peak_fc_data = permute(tfr_seed,[1 3 2]);
                            tfr_seed = squeeze(nanmean(tfr_seed(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_seed)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                        case 'Data (Noise thresholded)'
                            tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.plv_based(:,:,:),[1 3 2]);
                            mu = permute(h.inv_soln(h.current_inv_soln).TFR_results.Noise.plv_based(:,:,:),[1 3 2]);
                            sigma = repmat( nanstd(permute(h.inv_soln(h.current_inv_soln).TFR_results.Noise.plv_based(:,:,:),[1 3 2]),[],2) ,[1 size(mu,2) 1]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and noise
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            % creating significance mask
                            x_mask = nan(size(tfr_seed)); x_mask(tfr_seed>surg_thresh_pos)=1; x_mask(tfr_seed<surg_thresh_neg)=1;
                            tfr_seed = tfr_seed.*x_mask;
                            tfr_seed(isnan(tfr_seed))=nan;
                            h.current_inv_peak_fc_data = permute(tfr_seed,[1 3 2]);
                            tfr_seed = squeeze(nanmean(tfr_seed(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_seed)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                    end
                case 'PLI'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'; tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.pli_based(:,plv_true,:),[1 3 2]);
                        case 'Noise'; tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.Noise.pli_based(:,plv_true,:),[1 3 2]);
                        case 'Surrogate'; tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.pli_surg_based_mean(:,plv_true,:),[1 3 2]);
                        case 'Data (Surrogate thresholded)'
                            tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.pli_based(:,:,:),[1 3 2]);
                            mu = permute(h.inv_soln(h.current_inv_soln).TFR_results.pli_surg_based_mean(:,:,:),[1 3 2]);
                            sigma = permute(h.inv_soln(h.current_inv_soln).TFR_results.pli_surg_based_std(:,:,:),[1 3 2]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and surrogate
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            
                            % creating significance mask
                            x_mask = nan(size(tfr_seed)); x_mask(tfr_seed>surg_thresh_pos)=1; x_mask(tfr_seed<surg_thresh_neg)=1;
                            tfr_seed = tfr_seed.*x_mask;
                            tfr_seed(isnan(tfr_seed))=nan;
                            h.current_inv_peak_fc_data = permute(tfr_seed,[1 3 2]);
                            tfr_seed = squeeze(nanmean(tfr_seed(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_seed)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                        case 'Data (Noise thresholded)'
                            tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.pli_based(:,:,:),[1 3 2]);
                            mu = permute(h.inv_soln(h.current_inv_soln).TFR_results.Noise.pli_based(:,:,:),[1 3 2]);
                            sigma = repmat( nanstd(permute(h.inv_soln(h.current_inv_soln).TFR_results.Noise.pli_based(:,:,:),[1 3 2]),[],2) ,[1 size(mu,2) 1]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and noise
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            % creating significance mask
                            x_mask = nan(size(tfr_seed)); x_mask(tfr_seed>surg_thresh_pos)=1; x_mask(tfr_seed<surg_thresh_neg)=1;
                            tfr_seed = tfr_seed.*x_mask;
                            tfr_seed(isnan(tfr_seed))=nan;
                            h.current_inv_peak_fc_data = permute(tfr_seed,[1 3 2]);
                            tfr_seed = squeeze(nanmean(tfr_seed(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_seed)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                    end
                    lat = h.inv_soln(h.current_inv_soln).TFR_results.pli_lat;
                case 'dPLI'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'; tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.dpli_based(:,plv_true,:),[1 3 2]);
                        case 'Noise'; tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.Noise.dpli_based(:,plv_true,:),[1 3 2]);
                        case 'Surrogate'; tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.dpli_surg_based_mean(:,plv_true,:),[1 3 2]);
                        case 'Data (Surrogate thresholded)'
                            tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.dpli_based(:,:,:),[1 3 2]);
                            mu = permute(h.inv_soln(h.current_inv_soln).TFR_results.dpli_surg_based_mean(:,:,:),[1 3 2]);
                            sigma = permute(h.inv_soln(h.current_inv_soln).TFR_results.dpli_surg_based_std(:,:,:),[1 3 2]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and surrogate
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            
                            % creating significance mask
                            x_mask = nan(size(tfr_seed)); x_mask(tfr_seed>surg_thresh_pos)=1; x_mask(tfr_seed<surg_thresh_neg)=1;
                            tfr_seed = tfr_seed.*x_mask;
                            tfr_seed(isnan(tfr_seed))=nan;
                            h.current_inv_peak_fc_data = permute(tfr_seed,[1 3 2]);
                            tfr_seed = squeeze(nanmean(tfr_seed(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_seed)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                        case 'Data (Noise thresholded)'
                            tfr_seed = permute(h.inv_soln(h.current_inv_soln).TFR_results.dpli_based(:,:,:),[1 3 2]);
                            mu = permute(h.inv_soln(h.current_inv_soln).TFR_results.Noise.dpli_based(:,:,:),[1 3 2]);
                            sigma = repmat( nanstd(permute(h.inv_soln(h.current_inv_soln).TFR_results.Noise.dpli_based(:,:,:),[1 3 2]),[],2) ,[1 size(mu,2) 1]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and noise
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            % creating significance mask
                            x_mask = nan(size(tfr_seed)); x_mask(tfr_seed>surg_thresh_pos)=1; x_mask(tfr_seed<surg_thresh_neg)=1;
                            tfr_seed = tfr_seed.*x_mask;
                            tfr_seed(isnan(tfr_seed))=nan;
                            h.current_inv_peak_fc_data = permute(tfr_seed,[1 3 2]);
                            tfr_seed = squeeze(nanmean(tfr_seed(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_seed)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                    end
                    lat = h.inv_soln(h.current_inv_soln).TFR_results.pli_lat;
                case 'none'
                    tfr_seed = []; tfr_comp = []; h.axes_inv_soln_tfr.clo; h.axes_inv_soln_tfr.Title.String = 'No Results'; return;
            end
            tfr_seed(isnan(tfr_seed))=0; tfr_comp(isnan(tfr_comp))=0;
            tfr_data = cat(3,tfr_seed,tfr_comp);
            tfr_data = nanmean(tfr_data,3); % averaging across all selected peak sources
            
            tfr_freqs = h.inv_soln(h.current_inv_soln).TFR_results.TFR_freqs;
            
            tfr_caxis = str2num(h.edit_inv_tfr_caxis.String);
            h.axes_inv_soln_tfr.NextPlot='replace';
            h.current_inv_tfr_plot = surf(h.axes_inv_soln_tfr,lat,tfr_freqs,tfr_data); axis(h.axes_inv_soln_tfr,'tight');
            %     h.axes_inv_soln_tfr.Position(3)=.28;
            view(h.axes_inv_soln_tfr,0,90); shading(h.axes_inv_soln_tfr,'interp'); h.axes_inv_soln_tfr.Colormap=jet;
            h.axes_inv_soln_tfr.NextPlot='add';
            
            h.axes_inv_soln_tfr.XLim = str2num(h.edit_plot_time_int.String); h.axes_inv_soln_tfr.XLabel.String = 'Time (sec)';
            h.axes_inv_soln_tfr.YLim = str2num(h.edit_plot_freq_int.String); h.axes_inv_soln_tfr.YLabel.String = 'Frequency (Hz)';
            h.axes_inv_soln_tfr.CLim = tfr_caxis;
            axes(h.axes_inv_soln_tfr)
            h.colorbar_axes_inv_soln_tfr = colorbar('Location','eastoutside');
            
            plot3(h.axes_inv_soln_tfr,lat_coi,h.inv_soln(h.current_inv_soln).TFR_results.coi_wt2,ones(size(h.inv_soln(h.current_inv_soln).TFR_results.coi_wt2))*tfr_caxis(2),'color',[1 1 1]*.7,'linewidth',2);
            plot3(h.axes_inv_soln_tfr,[0 0],[h.axes_inv_soln_tfr.YLim],[1 1],'k--');
            %             x.Type = 'none'; sm_get_inv_tfr_point(x);
            
        end
        disableDefaultInteractivity(h.axes_inv_soln_tfr); h.axes_inv_soln_tfr.Toolbar.Visible='off';
        try h.axes_inv_soln_tfr.ButtonDownFcn = @sm_get_inv_tfr_point; end
        try h.current_inv_tfr_plot.ButtonDownFcn = @sm_get_inv_tfr_point; end
        S.Type = 'none'; sm_get_inv_tfr_point(S);
        h.axes_inv_soln_tfr.Position(3)=.3;
        
    end
end

