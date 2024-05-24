function sm_plot_true_tfr_connectivity(varargin)
global h

%% initiate variables
true_idx = 1; plv_true = 1; 
%%
update_listbox_plv_contrasts;
switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
    case {'Total Power' 'Evoked Power' 'Induced Power'}
        tfr_caxis = str2num(h.edit_inv_tfr_caxis.String);
        %% find all FCs associated with this TFR location
        true_idx = h.listbox_true_plv_contrasts.Value;
        splv_true =  find(sum(ismember(h.sim_data.cfg.source.plv_contrast_idx,true_idx),2)>0);
    case {'PLV' 'PLI' 'dPLI'}
        tfr_caxis = str2num(h.edit_inv_plv_caxis.String);
        % find only true indices of associated PLV
%                 v_idx = h.current_3D_plv_contrasts(h.current_3D_plv_contrasts_listbox_order(h.listbox_true_plv_contrasts.Value),:);
%                 plv_true =  find(ismember(h.current_3D_plv_contrasts(h.current_3D_plv_contrasts_seed_idx,:),v_idx,'rows')==1);
        true_idx = h.listbox_true_plv_contrasts.Value; %h.current_3D_plv_contrasts_listbox_order(h.listbox_true_plv_contrasts.Value);
        plv_true = true_idx;
end
%         v_idx = h.listbox_true_plv_contrasts.Value; %h.current_3D_plv_contrasts_listbox_order(h.listbox_true_plv_contrasts.Value);
%         true_idx = ismember(1:length(h.sim_data.cfg.source.vx_idx),v_idx);
%         h.current_true_plv_idx = find(true_idx==1);
%         plv_true = h.current_true_plv_idx;
% true_idx = v_idx; plv_true = v_idx; 
h.current_true_plv_idx = plv_true;

if h.menu_inv_tfr_type.Value<=3 && h.menu_inv_tfr_data_type.Value>=3; h.menu_inv_tfr_data_type.Value=1; end % Power data currently only have "Data" not "Noise" or "Surrogate" data yet
h.FC_alpha_level = str2num(h.edit_3D_FC_stats_alpha_level.String);

% turning things off
h.axes_true_source_tfr.Visible = 'off'; for a=1:length(h.axes_true_source_tfr.Children); h.axes_true_source_tfr.Children(a).Visible='off'; end
if isfield(h,'current_true_plv_plots'); if any(isvalid(h.current_true_plv_plots)); delete(h.current_true_plv_plots); end; end
if h.radio_inv_plot_peak_tfr_connectivity.Value==0 && h.radio_inv_plot_true_tfr_connectivity.Value==0
    h.axes_source_fft.Visible = 'on'; for a=1:length(h.axes_source_fft.Children); h.axes_source_fft.Children(a).Visible='on'; end
    h.menu_inv_tfr_type.Visible = 'off';
    h.menu_inv_tfr_data_type.Visible = 'off';
    h.edit_inv_tfr_caxis_txt.Visible = 'off';
    h.edit_inv_tfr_caxis.Visible = 'off';
    h.edit_inv_plv_caxis_txt.Visible = 'off';
    h.edit_inv_plv_caxis.Visible = 'off';
    h.radio_inv_plot_connectivity_graph.Visible ='off';
    if isfield(h,'current_inv_tfr_point_lines'); if any(isvalid(h.current_inv_tfr_point_lines)); delete(h.current_inv_tfr_point_lines); end; end
    if isfield(h,'current_inv_tfr_point_lines'); if any(isvalid(h.current_inv_tfr_point_lines)); delete(h.current_inv_tfr_point_lines); end; end
    h.edit_inv_plv_thresh_txt.Visible = 'off'; h.edit_inv_plv_thresh.Visible = 'off';
end

if ~isfield(h.sim_data.cfg.source,'TFR_results')
    if h.monte_carlo_flag==0
        msgbox(sprintf('Please perform Time-Frequency Response (TFR) analyses\n\nClick on "Calc True Connectivity"\n'));
    end
    h.radio_inv_plot_true_tfr_connectivity.Value = 0;
else
    if isempty(h.sim_data.cfg.source.TFR_results) && h.radio_inv_plot_true_tfr_connectivity.Value ==1
        if h.monte_carlo_flag==0
            msgbox(sprintf('Please perform Time-Frequency Response (TFR) analyses\n\nClick on "Calc True Connectivity"\n'));
        end
        h.radio_inv_plot_true_tfr_connectivity.Value = 0;
    else
        %% need to delete colorbar
        if isfield(h,'colorbar_axes_true_tfr'); if isvalid(h.colorbar_axes_true_tfr); delete(h.colorbar_axes_true_tfr); end; end
        
        if h.radio_inv_plot_true_tfr_connectivity.Value == 0
        else
            h.axes_true_source_tfr.Visible = 'on'; for a=1:length(h.axes_true_source_tfr.Children); h.axes_true_source_tfr.Children(a).Visible='on'; end
            h.axes_source_fft.Visible = 'off'; for a=1:length(h.axes_source_fft.Children); h.axes_source_fft.Children(a).Visible='off'; end
            if isfield(h,'current_true_plv_plots'); if any(isvalid(h.current_true_plv_plots)); delete(h.current_true_plv_plots); end; end
            h.menu_inv_tfr_type.Visible = 'on';
            h.menu_inv_tfr_data_type.Visible = 'on';
            h.edit_inv_tfr_caxis_txt.Visible = 'on';
            h.edit_inv_tfr_caxis.Visible = 'on';
            h.edit_inv_plv_caxis_txt.Visible = 'on';
            h.edit_inv_plv_caxis.Visible = 'on';
            h.radio_inv_plot_connectivity_graph.Visible ='on';
            h.edit_inv_plv_thresh_txt.Visible = 'on'; h.edit_inv_plv_thresh.Visible = 'on';
            
            %% start plotting
            % plot tfr for selected source in list box - if multiple selected then average across all
            
            if all(true_idx==0)
                h.axes_true_source_tfr.clo;
                text(h.axes_true_source_tfr,h.axes_true_source_tfr.XLim(1)+range(h.axes_true_source_tfr.XLim)/8,range(h.axes_true_source_tfr.YLim/2),sprintf('Time-Frequency not calculated\nfor selected source'));
                return
            end
            
            %% find all plv of selected peak sources
            %             plv_idx = nchoose2(1:3);
            %             plv_true = any(ismember(plv_idx,v_idx),2);
            %             h.current_true_plv_idx = find(plv_true==1);
            
            
            
            lat = h.sim_data.cfg.study.lat_sim;
            lat_coi = lat;
            %% TFR plot
            tfr_data = [];
            switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
                case 'Total Power'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'; tfr_data = h.sim_data.cfg.source.TFR_results.avg_true_wt(:,:,true_idx); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                        case 'Noise'; tfr_data = h.sim_data.cfg.source.TFR_results.Noise.avg_true_wt(:,:,true_idx); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                    end
                case 'Evoked Power'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'; tfr_data = h.sim_data.cfg.source.TFR_results.avg_true_wt_evk(:,:,true_idx); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                        case 'Noise'; tfr_data = h.sim_data.cfg.source.TFR_results.Noise.avg_true_wt_evk(:,:,true_idx); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                    end
                case 'Induced Power'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'; tfr_data = h.sim_data.cfg.source.TFR_results.avg_true_wt_ind(:,:,true_idx); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                        case 'Noise'; tfr_data = h.sim_data.cfg.source.TFR_results.Noise.avg_true_wt_ind(:,:,true_idx); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                    end
                case 'PLV'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'; tfr_data = permute(h.sim_data.cfg.source.TFR_results.plv_based(:,plv_true,:),[1 3 2]); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                        case 'Noise'; tfr_data = permute(h.sim_data.cfg.source.TFR_results.Noise.plv_based(:,plv_true,:),[1 3 2]); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                        case 'Surrogate'
                            if ~isempty(h.sim_data.cfg.source.TFR_results.plv_surg_based_mean)
                                tfr_data = permute(h.sim_data.cfg.source.TFR_results.plv_surg_based_mean(:,plv_true,:),[1 3 2]); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                            else; tfr_data = [];
                            end
                        case 'Data (Surrogate thresholded)'
                            tfr_data = permute(h.sim_data.cfg.source.TFR_results.plv_based(:,:,:),[1 3 2]);
                            mu = permute(h.sim_data.cfg.source.TFR_results.plv_surg_based_mean(:,:,:),[1 3 2]);
                            sigma = permute(h.sim_data.cfg.source.TFR_results.plv_surg_based_std(:,:,:),[1 3 2]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and surrogate
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            
                            % creating significance mask
                            x_mask = nan(size(tfr_data)); x_mask(tfr_data>surg_thresh_pos)=1; x_mask(tfr_data<surg_thresh_neg)=1;
                            tfr_data = tfr_data.*x_mask;
                            tfr_data(isnan(tfr_data))=0;
                            h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                            tfr_data = squeeze(nanmean(tfr_data(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_data)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                        case 'Data (Noise thresholded)'
                            tfr_data = permute(h.sim_data.cfg.source.TFR_results.plv_based(:,:,:),[1 3 2]);
                            mu = permute(h.sim_data.cfg.source.TFR_results.Noise.plv_based(:,:,:),[1 3 2]);
                            sigma = repmat( nanstd(permute(h.sim_data.cfg.source.TFR_results.Noise.plv_based(:,:,:),[1 3 2]),[],2) ,[1 size(mu,2) 1]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and noise
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            % creating significance mask
                            x_mask = nan(size(tfr_data)); x_mask(tfr_data>surg_thresh_pos)=1; x_mask(tfr_data<surg_thresh_neg)=1;
                            tfr_data = tfr_data.*x_mask;
                            tfr_data(isnan(tfr_data))=0;
                            h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                            tfr_data = squeeze(nanmean(tfr_data(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_data)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                    end
                    
                case 'PLI'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'; tfr_data = permute(h.sim_data.cfg.source.TFR_results.pli_based(:,plv_true,:),[1 3 2]); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                        case 'Noise'; tfr_data = permute(h.sim_data.cfg.source.TFR_results.Noise.pli_based(:,plv_true,:),[1 3 2]); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                        case 'Surrogate'; tfr_data = permute(h.sim_data.cfg.source.TFR_results.pli_surg_based_mean(:,plv_true,:),[1 3 2]); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                        case 'Data (Surrogate thresholded)'
                            tfr_data = permute(h.sim_data.cfg.source.TFR_results.pli_based(:,:,:),[1 3 2]);
                            mu = permute(h.sim_data.cfg.source.TFR_results.pli_surg_based_mean(:,:,:),[1 3 2]);
                            sigma = permute(h.sim_data.cfg.source.TFR_results.pli_surg_based_std(:,:,:),[1 3 2]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and surrogate
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            
                            % creating significance mask
                            x_mask = nan(size(tfr_data)); x_mask(tfr_data>surg_thresh_pos)=1; x_mask(tfr_data<surg_thresh_neg)=1;
                            tfr_data = tfr_data.*x_mask;
                            tfr_data(isnan(tfr_data))=0;
                            h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                            tfr_data = squeeze(nanmean(tfr_data(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_data)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                        case 'Data (Noise thresholded)'
                            tfr_data = permute(h.sim_data.cfg.source.TFR_results.pli_based(:,:,:),[1 3 2]);
                            mu = permute(h.sim_data.cfg.source.TFR_results.Noise.pli_based(:,:,:),[1 3 2]);
                            sigma = repmat( nanstd(permute(h.sim_data.cfg.source.TFR_results.Noise.pli_based(:,:,:),[1 3 2]),[],2) ,[1 size(mu,2) 1]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and noise
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            % creating significance mask
                            x_mask = nan(size(tfr_data)); x_mask(tfr_data>surg_thresh_pos)=1; x_mask(tfr_data<surg_thresh_neg)=1;
                            tfr_data = tfr_data.*x_mask;
                            tfr_data(isnan(tfr_data))=0;
                            h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                            tfr_data = squeeze(nanmean(tfr_data(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_data)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                    end
                    lat = h.sim_data.cfg.source.TFR_results.pli_lat;
                case 'dPLI'
                    switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                        case 'Data'; tfr_data = permute(h.sim_data.cfg.source.TFR_results.dpli_based(:,plv_true,:),[1 3 2]); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                        case 'Noise';  tfr_data = permute(h.sim_data.cfg.source.TFR_results.Noise.dpli_based(:,plv_true,:),[1 3 2]); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                        case 'Surrogate'; tfr_data = permute(h.sim_data.cfg.source.TFR_results.dpli_surg_based_mean(:,plv_true,:),[1 3 2]); h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                        case 'Data (Surrogate thresholded)'
                            tfr_data = permute(h.sim_data.cfg.source.TFR_results.dpli_based(:,:,:),[1 3 2]);
                            mu = permute(h.sim_data.cfg.source.TFR_results.dpli_surg_based_mean(:,:,:),[1 3 2]);
                            sigma = permute(h.sim_data.cfg.source.TFR_results.dpli_surg_based_std(:,:,:),[1 3 2]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and surrogate
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            
                            % creating significance mask
                            x_mask = nan(size(tfr_data)); x_mask(tfr_data>surg_thresh_pos)=1; x_mask(tfr_data<surg_thresh_neg)=1;
                            tfr_data = tfr_data.*x_mask;
                            tfr_data(isnan(tfr_data))=0;
                            h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                            tfr_data = squeeze(nanmean(tfr_data(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_data)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                        case 'Data (Noise thresholded)'
                            tfr_data = permute(h.sim_data.cfg.source.TFR_results.dpli_based(:,:,:),[1 3 2]);
                            mu = permute(h.sim_data.cfg.source.TFR_results.Noise.dpli_based(:,:,:),[1 3 2]);
                            sigma = repmat( nanstd(permute(h.sim_data.cfg.source.TFR_results.Noise.dpli_based(:,:,:),[1 3 2]),[],2) ,[1 size(mu,2) 1]);
                            x_std =tinv(1-h.FC_alpha_level,2); % using df = 2 because contrast is between data and noise
                            surg_thresh_pos = mu + (x_std*sigma);   % upper confidence interval
                            surg_thresh_neg = mu - (x_std*sigma);   % lower confidence interval
                            % creating significance mask
                            x_mask = nan(size(tfr_data)); x_mask(tfr_data>surg_thresh_pos)=1; x_mask(tfr_data<surg_thresh_neg)=1;
                            tfr_data = tfr_data.*x_mask;
                            tfr_data(isnan(tfr_data))=0;
                            h.current_inv_true_fc_data = permute(tfr_data,[1 3 2]);
                            tfr_data = squeeze(nanmean(tfr_data(:,:,plv_true),3));
                            %                             figure(2); clf; surf(squeeze(tfr_data)); view(0,90); shading interp; caxis([-.3 .3]); axis tight; colorbar; colormap(jet);
                    end
                    lat = h.sim_data.cfg.source.TFR_results.pli_lat;
                case 'none'
                    tfr_data = []; tfr_comp = []; h.axes_true_source_tfr.clo; h.axes_true_source_tfr.Title.String = 'No Results';  return;
            end
            tfr_data = squeeze(nanmean(tfr_data,3)); % averaging across selected sources
            tfr_freqs = h.sim_data.cfg.source.TFR_results.TFR_freqs;
            
            if ~isempty(tfr_data)
                %                 tfr_caxis = str2num(h.edit_inv_tfr_caxis.String);
                h.axes_true_source_tfr.NextPlot='replace';
                h.current_true_tfr_plot = surf(h.axes_true_source_tfr,lat,tfr_freqs,tfr_data); axis(h.axes_true_source_tfr,'tight');
                %     h.axes_true_source_tfr.Position(3)=.28;
                view(h.axes_true_source_tfr,0,90); shading(h.axes_true_source_tfr,'interp'); h.axes_true_source_tfr.Colormap=jet;
                h.axes_true_source_tfr.NextPlot='add';
                h.axes_true_source_tfr.XLim = str2num(h.edit_plot_time_int.String); h.axes_true_source_tfr.XLabel.String = 'Time (sec)';
                h.axes_true_source_tfr.YLim = str2num(h.edit_plot_freq_int.String); h.axes_true_source_tfr.YLabel.String = 'Frequency (Hz)';
                h.axes_true_source_tfr.CLim = tfr_caxis;
%                 axes(h.axes_true_source_tfr);
                h.colorbar_axes_true_tfr = colorbar(h.axes_true_source_tfr, 'Location','eastoutside');
                
                plot3(h.axes_true_source_tfr,lat_coi,h.sim_data.cfg.source.TFR_results.coi_wt2,ones(size(h.sim_data.cfg.source.TFR_results.coi_wt2))*tfr_caxis(2),'color',[1 1 1]*.7,'linewidth',2);
                plot3(h.axes_true_source_tfr,[0 0],[h.axes_true_source_tfr.YLim],[1 1],'k--');
                %             x.Type = 'none'; sm_get_inv_tfr_point(x);
            else
                h.axes_true_source_tfr.clo;
                text(h.axes_true_source_tfr,h.axes_true_source_tfr.XLim(1)+range(h.axes_true_source_tfr.XLim)/8,range(h.axes_true_source_tfr.YLim/2),sprintf('Time-Frequency not calculated\nfor selected source'));
            end
        end
        disableDefaultInteractivity(h.axes_true_source_tfr); h.axes_true_source_tfr.Toolbar.Visible='off';
        try  h.axes_true_source_tfr.ButtonDownFcn = @sm_get_inv_tfr_point; catch; end
        try  h.current_true_tfr_plot.ButtonDownFcn = @sm_get_inv_tfr_point; catch; end
        S.Type = 'none'; sm_get_inv_tfr_point(S);
        h.axes_true_source_tfr.Position(3)=.3;
        
    end
end
