function sm_plot_TFR_waves(varargin)
global h
% try
if ~isfield(h,'current_inv_tfr_freq_point'); h.current_inv_tfr_freq_point = h.sim_data.cfg.source.TFR_results.TFR_freqs(1); end
if isempty(h.current_inv_tfr_freq_point); h.current_inv_tfr_freq_point = h.sim_data.cfg.source.TFR_results.TFR_freqs(1); end

try
    ss=find(h.sim_data.cfg.source.TFR_results.TFR_freqs<=h.current_inv_tfr_freq_point);
catch
    ss=find(h.inv_soln(h.current_inv_soln).TFR_results.TFR_freqs<=h.current_inv_tfr_freq_point);
end

f_samp = ss(end);
% t_samp = h.current_inv_tfr_time_samp;
lat = h.cfg.study.lat_sim;

vx=1; tvx = 1; % initiating for legend
switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
    case {'Total Power' 'Evoked Power' 'Induced Power'}
        tfr_caxis = str2num(h.edit_inv_tfr_caxis.String);
    case {'PLV' 'PLI' 'dPLI'}
        tfr_caxis = str2num(h.edit_inv_plv_caxis.String);
end

true_clr = h.true_clr; % true color
peak_clr = h.peak_clr; % true color
nonsel_alpha = .2; % slpah transparency for non-selected waveforms

% h.axes_source_waves.clo;
h.axes_source_waves.NextPlot = 'replace';
% initializing variables for legends
p1=plot(h.axes_source_waves,1,1); delete(p1); p2=plot(h.axes_source_waves,1,1); delete(p2); p3=plot(h.axes_source_waves,1,1); delete(p3);
h.axes_source_waves.NextPlot = 'add';

other_lgnd ='Other Connections';
hit_lgnd =''; true_lgnd='';
% h = rmfield(h,'current_inv_TFR_waves');

if h.radio_plot_TFR_waves.Value==1
    %% %%%%% plotting inv_soln TFR waves
    %% find all plv of selected peak sources
    v_idx = h.current_3D_plv_contrasts_listbox_order(h.listbox_plv_contrasts.Value);
    %     plv_idx = h.inv_soln(h.current_inv_soln).plv_contrast_idx;
    %     plv_true = any(ismember(plv_idx,v_idx),2);
    %     h.current_plv_idx = find(plv_true==1);
    %     plv_true = h.current_plv_idx;
    
    % plotting zero axis lines
    plot(h.axes_source_waves,[0 0],h.axes_source_waves.YLim,'k-','linewidth',1);
    plot(h.axes_source_waves,h.axes_source_waves.XLim,[0 0],'k-','linewidth',1);
    
    if h.radio_inv_plot_peak_tfr_connectivity.Value == 1   % show the connectivity graph on the 3D maps
        switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
            case 'PLV'
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data';  
                        fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.plv_based(f_samp,:,:));
                    case 'Noise'
                        if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'Noise')
                            if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.Noise.plv_based(f_samp,:,:)); fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.Noise.plv_based(f_samp,:,:)); else; true_data = []; end
                        end
                    case 'Surrogate';   fc_data = squeeze(h.inv_soln(h.current_inv_soln).TFR_results.plv_surg_based_mean(f_samp,:,:));
                    case 'Data (Surrogate thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,:,:);
                    case 'Data (Noise thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,:,:);
                end
            case 'PLI'
                lat = h.inv_soln(h.current_inv_soln).TFR_results.pli_lat;
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data';  fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.pli_based(f_samp,:,:));
                    case 'Noise'
                        if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'Noise')
                            if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.Noise.pli_based(f_samp,:,:)); fc_data = squeeze(h.inv_soln(h.current_inv_soln).TFR_results.Noise.pli_based(f_samp,:,:)); else; true_data = []; end
                        else
                            fc_data = [];
                        end
                    case 'Surrogate';   fc_data = squeeze(h.inv_soln(h.current_inv_soln).TFR_results.pli_surg_based_mean(f_samp,:,:));
                    case 'Data (Surrogate thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,:,:);
                    case 'Data (Noise thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,:,:);
                end
            case 'dPLI'
                lat = h.inv_soln(h.current_inv_soln).TFR_results.pli_lat;
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data';  fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.dpli_based(f_samp,:,:));
                    case 'Noise'
                        if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'Noise')
                            if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.Noise.dpli_based(f_samp,:,:)); fc_data = squeeze(h.inv_soln(h.current_inv_soln).TFR_results.Noise.dpli_based(f_samp,:,:)); else; true_data = []; end
                        else
                            fc_data = [];
                        end
                    case 'Surrogate';   fc_data = squeeze(h.inv_soln(h.current_inv_soln).TFR_results.dpli_surg_based_mean(f_samp,:,:));
                    case 'Data (Noise thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,:,:);
                    case 'Data (Surrogate thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,:,:);
                end
            case 'Total Power'
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data'
                        %                         fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt(f_samp,:,:)));
                        fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt(f_samp,:,:));
                    case 'Noise'
                        if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'Noise')
                            %                             fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_seed_wt(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_comp_wt(f_samp,:,:)));
                            fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_comp_wt(f_samp,:,:));
                        else
                            fc_data =[];
                        end
                    case 'Surrogate'
                        %                         fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt(f_samp,:,:)));
                        fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt(f_samp,:,:));
                    case 'Data (Surrogate thresholded)'
                        %                         fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt(f_samp,:,:)));
                        fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt(f_samp,:,:));
                    case 'Data (Noise thresholded)'
                        %                         fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt(f_samp,:,:)));
                        fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt(f_samp,:,:));
                end
            case 'Evoked Power'
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data' ; fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_evk(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_evk(f_samp,:,:)));
                    case 'Noise'; fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_seed_wt_evk(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_comp_wt_evk(f_samp,:,:)));
                    case 'Surrogate'; fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_evk(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_evk(f_samp,:,:)));
                    case 'Data (Surrogate thresholded)'; fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_evk(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_evk(f_samp,:,:)));
                    case 'Data (Noise thresholded)'; fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_evk(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_evk(f_samp,:,:)));
                end
            case 'Induced Power'
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data' ; fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_ind(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_ind(f_samp,:,:)));
                    case 'Noise'; fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_seed_wt_ind(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.Noise.avg_comp_wt_ind(f_samp,:,:)));
                    case 'Surrogate'; fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_ind(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_ind(f_samp,:,:)));
                    case 'Data (Surrogate thresholded)'; fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_ind(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_ind(f_samp,:,:)));
                    case 'Data (Noise thresholded)'; fc_data =  cat(2,squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_seed_wt_ind(f_samp,:,:)),squeeze(h.inv_soln(h.current_inv_soln).TFR_results.avg_comp_wt_ind(f_samp,:,:)));
                end
            case {'none'}
                return
        end
        fc_data = squeeze(fc_data);
        fc_data(abs(fc_data)<str2num(h.edit_inv_plv_thresh.String))=nan;     % thresholding data
         if ~isempty(fc_data)
            if size(fc_data,1)~=length(lat); fc_data = fc_data'; end
            % plotting all PLV/TFR
            p1=plot(h.axes_source_waves,lat,fc_data,'Color','b','LineWidth',1);
            for v=1:length(p1); p1(v).Color = peak_clr(v,:); p1(v).Color(4) = nonsel_alpha; end
            h.current_inv_TFR_waves = p1;
            sel_idx = h.listbox_plv_contrasts.Value;
            for vx = 1:length(sel_idx)
                %                 h.current_inv_TFR_waves(vx)=plot(h.axes_source_waves,lat,fc_data(:,vx),'Color',peak_clr(vx,:),'LineWidth',2);
                set(h.current_inv_TFR_waves(h.current_3D_plv_contrasts_listbox_order(sel_idx(vx))),'Color',peak_clr(sel_idx(vx),:),'LineWidth',2);
            end
            switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
                case {'PLV' 'PLI' 'dPLI'}
                    h.axes_source_waves.YLabel.String = sprintf('%s',h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value});
                case {'Total Power' 'Evoked Power' 'Induced Power'}
                    h.axes_source_waves.YLabel.String = sprintf('%s (dB re:baseline)',h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value});
            end
            try h.axes_source_waves.Legend.Visible='off'; catch; end
        end
        
    end
    
    %% %%%%% plotting true TFR waves
    if h.radio_inv_plot_true_tfr_connectivity.Value == 1 % && isvalid(h.current_true_tfr_plot)   % also plot true source connectivity
        switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
            case 'PLV'
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data'
                        if ~isempty(h.sim_data.cfg.source.TFR_results.plv_based); true_data = squeeze(h.sim_data.cfg.source.TFR_results.plv_based(f_samp,:,:)); else; true_data = []; end
                    case 'Noise'
                        if isfield(h.sim_data.cfg.source.TFR_results,'Noise')
                            if ~isempty(h.sim_data.cfg.source.TFR_results.Noise.plv_based); true_data = squeeze(h.sim_data.cfg.source.TFR_results.Noise.plv_based(f_samp,:,:)); else; true_data = []; end
                        else
                            true_data = [];
                        end
                    case 'Surrogate'
                        if ~isempty(h.sim_data.cfg.source.TFR_results.plv_surg_based_mean); true_data = squeeze(h.sim_data.cfg.source.TFR_results.plv_surg_based_mean(f_samp,:,:)); else; true_data = []; end
                    case 'Data (Surrogate thresholded)'
                        try
                            if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,:,:); else; true_data = []; end
                        catch
                            sm_plot_true_tfr_connectivity;
                            if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,:,:); else; true_data = []; end
                        end
                    case 'Data (Noise thresholded)'
                        try
                            if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,:,:); else; true_data = []; end
                        catch
                            sm_plot_true_tfr_connectivity;
                            if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,:,:); else; true_data = []; end
                        end
                end
            case 'PLI'
                lat = h.inv_soln(h.current_inv_soln).TFR_results.pli_lat;
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.pli_based(f_samp,:,:));
                    case 'Noise'
                        if isfield(h.sim_data.cfg.source.TFR_results,'Noise')
                            if ~isempty(h.sim_data.cfg.source.TFR_results.Noise.pli_based); true_data = squeeze(h.sim_data.cfg.source.TFR_results.Noise.pli_based(f_samp,:,:)); else; true_data = []; end
                        else
                            true_data = [];
                        end
                    case 'Surrogate';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.pli_surg_based_mean(f_samp,:,:));
                    case 'Data (Surrogate thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,:,:);
                    case 'Data (Noise thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,:,:);
                end
            case 'dPLI'
                lat = h.inv_soln(h.current_inv_soln).TFR_results.pli_lat;
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.dpli_based(f_samp,:,:));
                    case 'Noise'
                        if isfield(h.sim_data.cfg.source.TFR_results,'Noise')
                            if ~isempty(h.sim_data.cfg.source.TFR_results.Noise.dpli_based); true_data = squeeze(h.sim_data.cfg.source.TFR_results.Noise.dpli_based(f_samp,:,:)); else; true_data = []; end
                        else
                            true_data = [];
                        end
                    case 'Surrogate';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.dpli_surg_based_mean(f_samp,:,:));
                    case 'Data (Surrogate thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,:,:);
                    case 'Data (Noise thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,:,:);
                end
            case 'Total Power'
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt(f_samp,:,:));
                    case 'Noise'
                        if isfield(h.sim_data.cfg.source.TFR_results,'Noise')
                            true_data = squeeze(h.sim_data.cfg.source.TFR_results.Noise.avg_true_wt(f_samp,:,:));
                        else
                            true_data = [];
                        end
                    case 'Surrogate';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt(f_samp,:,:));
                    case 'Data (Surrogate thresholded)'; true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt(f_samp,:,:));
                    case 'Data (Noise thresholded)';    true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt(f_samp,:,:));
                end
            case 'Evoked Power'
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt_evk(f_samp,:,:));
                    case 'Noise'
                        if isfield(h.sim_data.cfg.source.TFR_results,'Noise')
                            true_data = squeeze(h.sim_data.cfg.source.TFR_results.Noise.avg_true_wt_evk(f_samp,:,:));
                        else
                            true_data = [];
                        end
                        
                    case 'Surrogate';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt_evk(f_samp,:,:));
                    case 'Data (Surrogate thresholded)'; true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt_evk(f_samp,:,:));
                    case 'Data (Noise thresholded)';    true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt_evk(f_samp,:,:));
                end
            case 'Induced Power'
                switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                    case 'Data';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt_ind(f_samp,:,:));
                    case 'Noise'
                        if isfield(h.sim_data.cfg.source.TFR_results,'Noise')
                            true_data = squeeze(h.sim_data.cfg.source.TFR_results.Noise.avg_true_wt_ind(f_samp,:,:));
                        else
                            true_data = [];
                        end
                        
                    case 'Surrogate';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt_ind(f_samp,:,:));
                    case 'Data (Surrogate thresholded)'; true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt_ind(f_samp,:,:));
                    case 'Data (Noise thresholded)';    true_data = squeeze(h.sim_data.cfg.source.TFR_results.avg_true_wt_ind(f_samp,:,:));
                end
            case {'none'}
                return
        end
        
        %% Plotting
        ln_style = ':';
        true_data = squeeze(true_data);
        true_data(true_data==0)=nan;
        true_data(abs(true_data)<str2num(h.edit_inv_plv_thresh.String))=nan;
        tp1=plot(h.axes_source_waves,lat,true_data,'Color','b','LineWidth',2,'LineStyle',ln_style);
        for v=1:length(tp1); tp1(v).Color = true_clr(v,:); tp1(v).Color(4) = nonsel_alpha; end
        h.current_inv_true_TFR_waves = tp1;
        sel_idx = h.listbox_true_plv_contrasts.Value;
        for tvx = 1:length(sel_idx)
            %             h.current_inv_true_TFR_waves(tvx)=plot(h.axes_source_waves,lat,true_data(:,tvx),'Color',true_clr(tvx,:),'LineWidth',2,'LineStyle',ln_style);
            set(h.current_inv_true_TFR_waves(sel_idx(tvx)),'Color',true_clr(sel_idx(tvx),:),'LineWidth',2,'LineStyle',ln_style);
        end
        
        switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
            case {'PLV' 'PLI' 'dPLI'}
                h.axes_source_waves.YLabel.String = sprintf('%s',h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value});
            case {'Total Power' 'Evoked Power' 'Induced Power'}
                h.axes_source_waves.YLabel.String = sprintf('%s (dB re:baseline)',h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value});
        end
        try h.axes_source_waves.Legend.Visible='off'; catch; end
    end
    % plot time point line
    plot(h.axes_source_waves,[h.current_inv_tfr_time_point h.current_inv_tfr_time_point],h.axes_source_waves.YLim,'k:','linewidth',2);
    
    %     legend([p1(1) tp1(1)],{'Peak' 'True'});
    try
        if isvalid(h.current_inv_TFR_waves(vx)) && isvalid(h.current_inv_true_TFR_waves(tvx))
            legend([h.current_inv_TFR_waves(vx) h.current_inv_true_TFR_waves(tvx)],{'Peak' 'True'});
        elseif isvalid(h.current_inv_TFR_waves(vx)) && ~isvalid(h.current_inv_true_TFR_waves(tvx))
            legend([h.current_inv_TFR_waves(vx)],{'Peak'});
        elseif ~isvalid(h.current_inv_TFR_waves(vx)) && isvalid(h.current_inv_true_TFR_waves(tvx))
            legend([h.current_inv_true_TFR_waves(tvx)],{'True'});
        end
    catch
    end
    
    %% set axes limts
    h.axes_source_waves.YLim = tfr_caxis;
    h.axes_source_waves.XLim = str2num(h.edit_inv_source_waves_xscale.String);
else
    bs_plot_peak_waves;
end
% catch
%     msgbox(sprintf('Please perform Time-Frequency Response (TFR) analyses\n\nClick on "Calc Peak Connectivity"\n'));
% h.radio_plot_TFR_waves.Value = 0;
% end