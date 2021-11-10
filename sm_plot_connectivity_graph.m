function sm_plot_connectivity_graph(varargin)
global h

if h.radio_inv_plot_connectivity_graph.Value==0     % delete FC plot
    % deleting FC graph
    if isfield(h,'current_plv_plots'); if any(isvalid(h.current_plv_plots)); delete(h.current_plv_plots); end; end
    if isfield(h,'current_true_plv_plots'); if any(isvalid(h.current_true_plv_plots)); delete(h.current_true_plv_plots); end; end
    
elseif h.radio_inv_plot_connectivity_graph.Value==1  % replot FC grpah
    if isfield(h,'current_plv_plots'); if any(isvalid(h.current_plv_plots)); delete(h.current_plv_plots); end; end
    if isfield(h,'current_true_plv_plots'); if any(isvalid(h.current_true_plv_plots)); delete(h.current_true_plv_plots); end; end
    if ~isfield(h,'current_inv_tfr_freq_samp')    % set current_time_freq_point
        hm = msgbox(sprintf('Please select a Time-Frequency point by clicking on the TFR plot'),'Select Time-Frequency Point');
    else
        
        if h.radio_inv_plot_connectivity_graph.Value == 1   % show the connectivity graph on the 3D maps
            
            %% PLotting Peak Source Connection graph
            if h.radio_inv_plot_peak_tfr_connectivity.Value == 1
                
                % get connectivity data
                plv_true = h.current_plv_idx;
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
                fc_locs = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(vx_cont_idx,:);
                fc_locs = reshape(fc_locs,[size(fc_locs,1)/2 2 3]); % [plv_idx x contrast x XYZ]
                
                %% plotting connectivity graph
                % deleting old FC graph if it exists
                if isfield(h,'current_plv_plots'); if isvalid(h.current_plv_plots); delete(h.current_plv_plots); end; end
                
                % plotting new FC graph
                for v = 1:length(fc_data)
                    if  abs(fc_data(v)) > str2num(h.edit_inv_plv_thresh.String)
                        cx_locs = squeeze(fc_locs(v,:,:));
                        ln_width = abs(fc_data(v)*15);
                        if fc_data(v)>0; ln_clr = 'r'; elseif fc_data(v)<=0; ln_clr = 'b'; end % ERS=red and ERD=blue
                        h.current_plv_plots(v) = plot3(h.axes_3D_images,cx_locs(:,1),cx_locs(:,2),cx_locs(:,3),'LineWidth',ln_width,'Color',ln_clr);
                        
                        alpha_val=abs(fc_data(v)*2); if alpha_val>1; alpha_val=1; end
                        h.current_plv_plots(v).Color(4) = alpha_val;
                    end
                end
            end
            
            if h.radio_inv_plot_true_tfr_connectivity.Value == 1    % also plot true source connectivity
                ss=find(h.cfg.source.TFR_results.TFR_freqs<=h.current_inv_tfr_freq_point);
                f_samp = ss(end);
%                 f_samp = h.current_inv_tfr_freq_samp;
                t_samp = h.current_inv_tfr_time_samp;
                
                switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
                    case 'PLV'
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data'  
                                if ~isempty(h.cfg.source.TFR_results.plv_based); true_data = squeeze(h.cfg.source.TFR_results.plv_based(f_samp,:,t_samp)); else; true_data = []; end 
                            case 'Noise'; true_data = [];
                            case 'Surrogate'
                                if ~isempty(h.cfg.source.TFR_results.plv_surg_based_mean); true_data = squeeze(h.cfg.source.TFR_results.plv_surg_based_mean(f_samp,:,t_samp)); else; true_data = []; end
                            case 'Data (Surrogate thresholded)'
                                try
                                    if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,:,t_samp); else; true_data = []; end
                                catch
                                    sm_plot_true_tfr_connectivity;
                                    if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,:,t_samp); else; true_data = []; end
                                end
                            case 'Data (Noise thresholded)'
                                try 
                                    if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,:,t_samp); else; true_data = []; end 
                                catch
                                   sm_plot_true_tfr_connectivity;
                                   if ~isempty(h.current_inv_true_fc_data); true_data = h.current_inv_true_fc_data(f_samp,:,t_samp); else; true_data = []; end
                                end
                        end
                    case 'PLI'
                        ss = find(h.cfg.source.TFR_results.pli_lat<h.current_inv_tfr_time_point); t_samp=ss(end);
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data';  true_data = squeeze(h.cfg.source.TFR_results.pli_based(f_samp,:,t_samp));
                            case 'Noise'; true_data = [];
                            case 'Surrogate';  true_data = squeeze(h.cfg.source.TFR_results.pli_surg_based_mean(f_samp,:,t_samp));
                            case 'Data (Surrogate thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,:,t_samp);
                            case 'Data (Noise thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,:,t_samp);
                        end
                    case 'dPLI'
                        ss = find(h.cfg.source.TFR_results.pli_lat<h.current_inv_tfr_time_point); t_samp=ss(end);
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data';  true_data = squeeze(h.cfg.source.TFR_results.dpli_based(f_samp,:,t_samp));
                            case 'Noise'; true_data = [];
                            case 'Surrogate';  true_data = squeeze(h.cfg.source.TFR_results.dpli_surg_based_mean(f_samp,:,t_samp));
                            case 'Data (Surrogate thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,:,t_samp);
                            case 'Data (Noise thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,:,t_samp);
                        end
                    case {'Total Power' 'Evoked Power' 'Induced Power'}
                        % defaulting to PLV
                        switch h.menu_inv_tfr_data_type.String{h.menu_inv_tfr_data_type.Value}
                            case 'Data';  true_data = squeeze(h.cfg.source.TFR_results.plv_based(f_samp,:,t_samp));
                            case 'Noise'; true_data = [];
                            case 'Surrogate';  true_data = squeeze(h.cfg.source.TFR_results.plv_surg_based_mean(f_samp,:,t_samp));
                            case 'Data (Surrogate thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,:,t_samp);
                            case 'Data (Noise thresholded)'; true_data = h.current_inv_true_fc_data(f_samp,:,t_samp);
                        end
                end
                true_cont_idx = nchoose2(1:3);
                for vc=1:length(h.current_true_plv_idx)  % looping through true_data --> plv indices
                    v=h.current_true_plv_idx(vc);
                    vx = true_cont_idx(h.current_true_plv_idx(vc),:);    % voxel contrasts
                    if ~isempty(true_data)
                        if  abs(true_data(v)) > str2num(h.edit_inv_plv_thresh.String)
                            cx_locs = squeeze(h.cfg.source.vx_locs(vx,:));
                            ln_width = abs(true_data(v)*15);
                            if true_data(v)>0; ln_clr = 'r'; elseif true_data(v)<=0; ln_clr = 'b'; end % ERS=red and ERD=blue
                            h.current_true_plv_plots(v) = plot3(h.axes_3D_images,cx_locs(:,1),cx_locs(:,2),cx_locs(:,3),'-.','LineWidth',ln_width,'Color',ln_clr);
                            alpha_val=abs(true_data(v)*2); if alpha_val>1; alpha_val=1; end
                            h.current_true_plv_plots(v).Color(4) = alpha_val;
                        end
                    end
                end
            end
        end
    end
end


