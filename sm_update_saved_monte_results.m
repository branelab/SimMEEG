function sm_update_saved_monte_results(varargin)
global h
% update listboxes with fieldnames

%% sim_data
h.listbox_monte_saved_sim_data.String = fieldnames(h.sim_data);
%% True Source
h.listbox_monte_saved_true_source.String = fieldnames(h.sim_data.cfg.source);

if isfield(h.sim_data.cfg.source(1),'TFR_results')
    if ~isempty(h.sim_data.cfg.source(1).TFR_results)
        h.listbox_monte_saved_true_source_TFR_results.String = fieldnames(h.sim_data.cfg.source(1).TFR_results);
    end
end


%% inv_soln
h.listbox_monte_saved_inv_soln.String = fieldnames(h.inv_soln);
if isfield(h,'inv_soln')
   if ~isempty(h.inv_soln)
    fnames = [];
    if isfield(h.inv_soln,'TFR_results') && ~isempty(h.inv_soln)
        if ~isempty(h.inv_soln(1).TFR_results)
            fnames = fieldnames(h.inv_soln(1).TFR_results);
        end
    end
    h.listbox_monte_saved_inv_soln_TFR_results.String = fnames;
    
    fnames = [];
    for a=1:length(h.inv_soln)
        fnames = [fnames fieldnames(h.inv_soln(a).soln)'];
    end
    fnames = unique(fnames);
    
    h.listbox_monte_saved_inv_soln_soln.String = fnames;
   else
      h.listbox_monte_saved_inv_soln.String = [];
      h.listbox_monte_saved_inv_soln_soln.String = [];
    h.listbox_monte_saved_inv_soln_TFR_results.String = [];
   end     
end

switch varargin{end}
    case 'default' % reset to selecting all fieldnames
        %% select deafult fieldnames - basic metrics --> loc error ...
        %% sim_data
        def_names = {'cfg', 'sens_final' 'sens_noise_final' 'sig_final' 'source_waveform_type' 'intersource_correlations' 'brain_noise' 'noise_gain'};
        sel_idx = ismember(h.listbox_monte_saved_sim_data.String,def_names);
        h.listbox_monte_saved_sim_data.Value = find(sel_idx==1);
        %% true source data
        def_names = {'sig_freqs' 'vx_ori' 'vx_idx' 'vx_amp' 'src_clr' 'vx_locs' 'true_evk_fft_data' 'true_evk_fft_freqs' 'plv_contrast_idx' 'plv_contrasts' 'TFR_results'};
        sel_idx = ismember(h.listbox_monte_saved_true_source.String,def_names);
        h.listbox_monte_saved_true_source.Value = find(sel_idx==1);
        h.listbox_monte_saved_true_source_TFR_results.Visible = 'off'; h.listbox_monte_saved_true_source_TFR_results_txt.Visible = 'off';
        def_names = {'plv_surg_based_mean' 'plv_surg_based_std'  'pli_surg_based_mean' 'pli_surg_based_std'  'dpli_surg_based_mean' 'dpli_surg_based_std' 'plv_based' 'pli_based' 'dpli_based' 'pli_lat' 'TFR_freqs' 'vx_amp' 'avg_true_wt' 'PLV_freqs' 'Noise'};
        sel_idx = ismember(h.listbox_monte_saved_true_source_TFR_results.String,def_names);
        h.listbox_monte_saved_true_source_TFR_results.Value = find(sel_idx==1);
        sm_update_saved_monte_results('true_source');
        
        %% inv_soln
        def_names = {'sens' 'headmodel_type' 'maxvectorori_Type' 'params' 'Type' 'soln' 'ListBox_name' 'org_img' 'SD_PSF' 'SD_CTF' 'peak_voxels' 'peak_idx' 'half_width' 'peak_voxels' 'classifier_metrics' 'classifier_performance' 'ori_error_mse' 'wave_error_mse_norm_act' 'wave_error_mse_norm_ctrl' 'wave_error_mse_norm_type' 'half_width' 'classifier_results' 'TFR_results' 'swf_evk_fft_data' 'swf_evk_fft_freqs' 'plv_leadfield_grid_idx' 'plv_seed_idx' 'plv_comp_idx' 'plv_contrast_idx' 'plv_contrasts' 'plv_seed_idx_in_plv_contrast_idx'};
        sel_idx = ismember(h.listbox_monte_saved_inv_soln.String,def_names);
        h.listbox_monte_saved_inv_soln.Value = find(sel_idx==1);
        
        h.listbox_monte_saved_inv_soln_TFR_results.Visible = 'off'; h.listbox_monte_saved_inv_soln_TFR_results_txt.Visible = 'off';
        def_names = {'plv_surg_based_mean' 'plv_surg_based_std'  'pli_surg_based_mean' 'pli_surg_based_std'  'dpli_surg_based_mean' 'dpli_surg_based_std' 'avg_seed_wt' 'avg_comp_wt' 'Noise' 'plv_based' 'pli_based' 'dpli_based' 'pli_lat' 'TFR_freqs' 'PLV_freqs' 'TotPwr_ROI_FC_stepwise_metrics' 'EvkPwr_ROI_FC_stepwise_metrics' 'IndPwr_ROI_FC_stepwise_metrics' 'PLV_ROI_FC_stepwise_metrics' 'PLV_ROI_FC_surg_std_metrics' 'PLI_ROI_FC_stepwise_metrics' 'PLI_ROI_FC_surg_std_metrics'};
        sel_idx = ismember(h.listbox_monte_saved_inv_soln_TFR_results.String,def_names);
        h.listbox_monte_saved_inv_soln_TFR_results.Value = find(sel_idx==1);
        
        h.listbox_monte_saved_inv_soln_soln.Visible = 'on'; h.listbox_monte_saved_inv_soln_soln_txt.Visible = 'on';
        def_names = {'Options' 'MCMV_idx' 'MCMV_ori' 'mHref' 'max_ori' 'method' 'null_thresh' 'nulled_Hscalar' 'nulled_ori' 'nulled_thresh' 'nulled_wts'...
            'peak_idx_found' 'time' 'true_source_idx' 'vectori_type' 'ori' 'ori_org' 'P' 'P_ctrl' 'cfg' 'error_locs' 'error_ori' 'error_waves' 'inside' 'mehtod' 'RNcov' 'null_thresh' 'compute_time' 'plot_min_max' 'plot_thresh' 'wts' 'wts_org'};
        sel_idx = ismember( h.listbox_monte_saved_inv_soln_soln.String,def_names);
        h.listbox_monte_saved_inv_soln_soln.Value = find(sel_idx==1);
        sm_update_saved_monte_results('inv_soln');

        
    case 'all' % reset to selecting all fieldnames
        h.listbox_monte_saved_sim_data.Value = 1:length(h.listbox_monte_saved_sim_data.String);
        h.listbox_monte_saved_true_source.Value = 1:length(h.listbox_monte_saved_true_source.String);
        h.listbox_monte_saved_true_source_TFR_results.Value = 1:length(h.listbox_monte_saved_true_source_TFR_results.String);
        h.listbox_monte_saved_inv_soln.Value = 1:length(h.listbox_monte_saved_inv_soln.String);
        h.listbox_monte_saved_inv_soln_TFR_results.Value = 1:length(h.listbox_monte_saved_inv_soln_TFR_results.String);
        h.listbox_monte_saved_inv_soln_soln.Value = 1:length(h.listbox_monte_saved_inv_soln_soln.String);
    case 'true_source'  % hiding soln & TFR results listboxes if they are not selected
        sel_names = h.listbox_monte_saved_true_source.String(h.listbox_monte_saved_true_source.Value);
        if any(ismember(sel_names,'TFR_results'))
            h.listbox_monte_saved_true_source_TFR_results.Visible = 'on'; h.listbox_monte_saved_true_source_TFR_results_txt.Visible = 'on';
        else
            h.listbox_monte_saved_true_source_TFR_results.Visible = 'off'; h.listbox_monte_saved_true_source_TFR_results_txt.Visible = 'off';
        end
        
    case 'inv_soln'
        sel_names = h.listbox_monte_saved_inv_soln.String(h.listbox_monte_saved_inv_soln.Value);
        if any(ismember(sel_names,'soln'))
            h.listbox_monte_saved_inv_soln_soln.Visible = 'on'; h.listbox_monte_saved_inv_soln_soln_txt.Visible = 'on';
        else
            h.listbox_monte_saved_inv_soln_soln.Visible = 'off'; h.listbox_monte_saved_inv_soln_soln_txt.Visible = 'off';
        end
        
        if any(ismember(sel_names,'TFR_results'))
            h.listbox_monte_saved_inv_soln_TFR_results.Visible = 'on'; h.listbox_monte_saved_inv_soln_TFR_results_txt.Visible = 'on';
        else
            h.listbox_monte_saved_inv_soln_TFR_results.Visible = 'off'; h.listbox_monte_saved_inv_soln_TFR_results_txt.Visible = 'off';
        end
        
end

