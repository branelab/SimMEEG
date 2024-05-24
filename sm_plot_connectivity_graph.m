function sm_plot_connectivity_graph(varargin)
global h

if h.radio_inv_plot_connectivity_graph.Value==0     % delete FC plot
    % deleting FC graph
    if isfield(h,'current_plv_plots'); if any(isvalid(h.current_plv_plots)); delete(h.current_plv_plots); end; end
    if isfield(h,'current_true_plv_plots'); if any(isvalid(h.current_true_plv_plots)); delete(h.current_true_plv_plots); end; end
    
elseif h.radio_inv_plot_connectivity_graph.Value==1  % replot FC graph
    if isfield(h,'current_plv_plots'); if any(isvalid(h.current_plv_plots)); delete(h.current_plv_plots); end; end
    if isfield(h,'current_true_plv_plots'); if any(isvalid(h.current_true_plv_plots)); delete(h.current_true_plv_plots); end; end
    if ~isfield(h,'current_inv_tfr_freq_samp')    % set current_time_freq_point
        hm = msgbox(sprintf('Please select a Time-Frequency point by clicking on the TFR plot'),'Select Time-Frequency Point');
    else
        
        if h.radio_inv_plot_connectivity_graph.Value == 1   % show the connectivity graph on the 3D maps
            
            %% Plotting Peak Source Connection graph
            if h.radio_inv_plot_peak_tfr_connectivity.Value == 1
                
                % get connectivity data
                switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
                    case {'Total Power' 'Evoked Power' 'Induced Power'}
                        %% find all FCs associated with this TFR location
                        v_idx = h.current_3D_plv_contrasts(h.current_3D_plv_contrasts_listbox_order(h.listbox_plv_contrasts.Value));
                        plv_true =  find(sum(ismember(h.inv_soln(h.current_inv_soln).plv_contrasts,v_idx),2)>0);
%                         v_idx = h.listbox_plv_contrasts.Value; %h.current_3D_plv_contrasts_listbox_order(h.listbox_plv_contrasts.Value);
%                         plv_true = v_idx;
                    case {'PLV' 'PLI' 'dPLI'}
%                         plv_true = h.current_3D_plv_contrasts_listbox_order(h.listbox_plv_contrasts.Value);
                        v_idx = h.listbox_plv_contrasts.Value; %h.current_3D_plv_contrasts_listbox_order(h.listbox_plv_contrasts.Value);
                        plv_true = h.current_3D_plv_contrasts_listbox_order(v_idx);
                end
                
                f_samp = h.current_inv_tfr_freq_samp;
                t_samp = h.current_inv_tfr_time_samp;
                
                switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
                    case 'PLV'
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data';  fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.plv_based(f_samp,plv_true,t_samp));
                            case 'Noise'; fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.Noise.plv_based(f_samp,plv_true,t_samp));
                            case 'Surrogate';   fc_data = squeeze(h.inv_soln(h.current_inv_soln).TFR_results.plv_surg_based_mean(f_samp,plv_true,t_samp));
                            case 'Data (Surrogate thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,plv_true,t_samp);
                             case 'Data (Noise thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,plv_true,t_samp);
                       end
                    case 'PLI'
                        ss = find(h.inv_soln(h.current_inv_soln).TFR_results.pli_lat<h.current_inv_tfr_time_point); t_samp=ss(end);
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data';  fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.pli_based(f_samp,plv_true,t_samp));
                            case 'Noise'; fc_data = [];
                            case 'Surrogate';   fc_data = squeeze(h.inv_soln(h.current_inv_soln).TFR_results.pli_surg_based_mean(f_samp,plv_true,t_samp));
                            case 'Data (Surrogate thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,plv_true,t_samp);
                            case 'Data (Noise thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,plv_true,t_samp);
                        end
                    case 'dPLI'
                        ss = find(h.inv_soln(h.current_inv_soln).TFR_results.pli_lat<h.current_inv_tfr_time_point); t_samp=ss(end);
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data';  fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.dpli_based(f_samp,plv_true,t_samp));
                            case 'Noise'; fc_data = [];
                            case 'Surrogate';   fc_data = squeeze(h.inv_soln(h.current_inv_soln).TFR_results.dpli_surg_based_mean(f_samp,plv_true,t_samp));
                            case 'Data (Noise thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,plv_true,t_samp);
                            case 'Data (Surrogate thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,plv_true,t_samp);
                        end
                    case {'Total Power' 'Evoked Power' 'Induced Power'}
                        % defaulting to PLV
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data';  fc_data =  squeeze(h.inv_soln(h.current_inv_soln).TFR_results.plv_based(f_samp,plv_true,t_samp));
                            case 'Noise'; fc_data = [];
                            case 'Surrogate';   fc_data = squeeze(h.inv_soln(h.current_inv_soln).TFR_results.plv_surg_based_mean(f_samp,plv_true,t_samp));
                            case 'Data (Surrogate thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,plv_true,t_samp);
                            case 'Data (Noise thresholded)'; fc_data = h.current_inv_peak_fc_data(f_samp,plv_true,t_samp);
                        end
                    case {'none'}
                        return
                end
                
                %% get fc_locs
                vx_cont_idx = h.inv_soln(h.current_inv_soln).plv_contrasts(plv_true,:);
%                 vx_cont_idx = vx_cont_idx(sum(~isnan(vx_cont_idx),2)==2,:);
                valid_idx = ~isnan(vx_cont_idx);
                fc_locs = nan(size(plv_true,1)*2,3);
                fc_locs(valid_idx,:) = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(vx_cont_idx(valid_idx),:);
                fc_locs = reshape(fc_locs,[size(fc_locs,1)/2 2 3]); % [plv_idx x contrast x XYZ]
                
                %% plotting connectivity graph
                % deleting old FC graph if it exists
                if isfield(h,'current_plv_plots'); if isvalid(h.current_plv_plots); delete(h.current_plv_plots); end; end
                
                % plotting new FC graph
                for v = 1:size(fc_data,2)
                    if  abs(fc_data(v)) > str2num(h.edit_inv_plv_thresh.String)
                        cx_locs = squeeze(fc_locs(v,:,:));
%                         ln_width = abs(fc_data(v)*15);
                        ln_width = abs(fc_data(v))/max(str2num(h.edit_inv_plv_caxis.String))*4;
                        if fc_data(v)>0; ln_clr = 'r'; elseif fc_data(v)<=0; ln_clr = 'b'; end % ERS=red and ERD=blue
                        h.current_plv_plots(v) = plot3(h.axes_3D_images,cx_locs(:,1),cx_locs(:,2),cx_locs(:,3),'LineWidth',ln_width,'Color',ln_clr);
                        
%                         alpha_val=abs(fc_data(v)*2); if alpha_val>1; alpha_val=1; end
                        alpha_val=abs(fc_data(v))/max(str2num(h.edit_inv_plv_caxis.String)); if alpha_val>1; alpha_val=1; end
                        h.current_plv_plots(v).Color(4) = alpha_val;
                    end
                end
            end
            
             %% Plotting True Source Connection graph
           if h.radio_inv_plot_true_tfr_connectivity.Value == 1    % also plot true source connectivity
                ss=find(h.sim_data.cfg.source.TFR_results.TFR_freqs<=h.current_inv_tfr_freq_point);
                f_samp = ss(end);
%                 f_samp = h.current_inv_tfr_freq_samp;
                t_samp = h.current_inv_tfr_time_samp;
                
                true_idx = h.current_true_plv_idx; %h.listbox_true_plv_contrasts.Value; 
                
                switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
                    case 'PLV'
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data'  
                                if ~isempty(h.sim_data.cfg.source.TFR_results.plv_based); true_data = squeeze(h.sim_data.cfg.source.TFR_results.plv_based(f_samp,true_idx,t_samp)); else; true_data = []; end 
                            case 'Noise'; true_data = [];
                            case 'Surrogate'
                                if ~isempty(h.sim_data.cfg.source.TFR_results.plv_surg_based_mean); true_data = squeeze(h.sim_data.cfg.source.TFR_results.plv_surg_based_mean(f_samp,true_idx,t_samp)); else; true_data = []; end
                            case 'Data (Surrogate thresholded)'
                                try
                                    if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,true_idx,t_samp); else; true_data = []; end
                                catch
                                    sm_plot_true_tfr_connectivity;
                                    if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,true_idx,t_samp); else; true_data = []; end
                                end
                            case 'Data (Noise thresholded)'
                                try 
                                    if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,true_idx,t_samp); else; true_data = []; end 
                                catch
                                   sm_plot_true_tfr_connectivity;
                                   if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,true_idx,t_samp); else; true_data = []; end
                                end
                        end
                    case 'PLI'
                        ss = find(h.sim_data.cfg.source.TFR_results.pli_lat<h.current_inv_tfr_time_point); t_samp=ss(end);
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.pli_based(f_samp,true_idx,t_samp));
                            case 'Noise'; true_data = [];
                            case 'Surrogate';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.pli_surg_based_mean(f_samp,true_idx,t_samp));
                            case 'Data (Surrogate thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,true_idx,t_samp);
                            case 'Data (Noise thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,true_idx,t_samp);
                        end
                    case 'dPLI'
                        ss = find(h.sim_data.cfg.source.TFR_results.pli_lat<h.current_inv_tfr_time_point); t_samp=ss(end);
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.dpli_based(f_samp,true_idx,t_samp));
                            case 'Noise'; true_data = [];
                            case 'Surrogate';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.dpli_surg_based_mean(f_samp,true_idx,t_samp));
                            case 'Data (Surrogate thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,true_idx,t_samp);
                            case 'Data (Noise thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,true_idx,t_samp);
                        end
                    case {'Total Power' 'Evoked Power' 'Induced Power'}
                        % defaulting to PLV
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.plv_based(f_samp,true_idx,t_samp));
                            case 'Noise'; true_data = [];
                            case 'Surrogate';  true_data = squeeze(h.sim_data.cfg.source.TFR_results.plv_surg_based_mean(f_samp,true_idx,t_samp));
                            case 'Data (Surrogate thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,true_idx,t_samp);
                            case 'Data (Noise thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,true_idx,t_samp);
                        end
                end
                
                for vc=1:length(true_idx) %1:h.current_3D_plv_contrasts_seed_idx %% length(h.current_true_plv_idx)  % looping through true_data --> plv indices
                  
                    if ~isempty(true_data)
                        if  abs(true_data(vc)) > str2num(h.edit_inv_plv_thresh.String)
                            vx = h.sim_data.cfg.source.plv_contrast_idx(true_idx(vc),:);    % voxel contrasts
                            cx_locs = squeeze(h.cfg.source.vx_locs(vx,:));
%                             ln_width = abs(true_data(vc)*2);
                            ln_width = abs(true_data(vc))/max(str2num(h.edit_inv_plv_caxis.String))*4;
                            
                            if true_data(vc)>0; ln_style = '-'; elseif true_data(vc)<=0; ln_style = ':'; end % ERS=red and ERD=blue
                            h.current_true_plv_plots(vc) = plot3(h.axes_3D_images,cx_locs(:,1),cx_locs(:,2),cx_locs(:,3),'LineStyle',ln_style,'LineWidth',ln_width,'Color','k'); %ln_clr);
%                             if true_data(vc)>0; ln_clr = 'k'; ln_style = '-'; elseif true_data(vc)<=0; ln_clr = 'k'; ln_style = '-.'; end % ERS=red and ERD=blue
%                             h.current_true_plv_plots(vc) = plot3(h.axes_3D_images,cx_locs(:,1),cx_locs(:,2),cx_locs(:,3),'LineStyle',ln_style,'LineWidth',ln_width,'Color',ln_clr);
                            alpha_val=abs(true_data(vc)*2); if alpha_val>1; alpha_val=1; end
                            h.current_true_plv_plots(vc).Color(4) = alpha_val;
                        end
                    end
                end
            end
        end
    end
end

if h.radio_plot_TFR_waves.Value ==1
    sm_plot_TFR_waves;
end
