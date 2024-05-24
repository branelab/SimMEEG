function update_image_thresh_txt(varargin)
global h

if h.load_study_flag~=1
    
    pos_ratio = ( h.slider_3D_image_thresh.Value-h.slider_3D_image_thresh.Min) / range([h.slider_3D_image_thresh.Min h.slider_3D_image_thresh.Max]);
    if isnan(pos_ratio)==1; pos_ratio=0; end
    s_pos = h.slider_3D_image_thresh.Position;
    h.slider_3D_image_thresh_text_val.Position(1) = s_pos(1) - h.slider_3D_image_thresh_text_val.Position(3) - .01;
    h.slider_3D_image_thresh_text_val.Position(2) = (range(h.slider_3D_image_thresh.Position([2 4]))*pos_ratio*.625) + h.slider_3D_image_thresh.Position(2);
    
    if h.slider_3D_image_thresh.Value<.001
        txt = compose("%.1e", h.slider_3D_image_thresh.Value);  % scientific notation
        h.slider_3D_image_thresh_text_val.Position([1 3]) = [.8 .115];
    else
        txt = compose("%.3f", h.slider_3D_image_thresh.Value);  % reg notation
        h.slider_3D_image_thresh_text_val.Position([1 3]) = [.86 .055];
    end
    
    min_max = [h.slider_3D_image_thresh.Min h.slider_3D_image_thresh.Max]; %str2num(h.edit_3D_min_max.String);
    
    if min_max(1)<.001; txt_min = compose("%.1e",min_max(1)); else txt_min = compose("%.3f",min_max(1)); end
    if min_max(2)<.001; txt_max = compose("%.1e",min_max(2)); else txt_max= compose("%.3f",min_max(2)); end
    if min_max(1)==0; txt_min='0'; end
    txt_minmax = [char(txt_min) ' ' char(txt_max)];
    
    h.slider_3D_image_thresh_text_val.String = txt;
    h.slider_3D_image_thresh_text_max.String = txt_max;
    
    h.edit_3D_min_max.String = txt_minmax;
    
    % updating Perforamnce Results red hashed threshold line based on slider_3D  value
    if h.run_inv_soln_flag == 0 % only change and plot after run_source_modeling.m has been run
        
        try
            h.line_thresh_invSoln_perf_ROC_stepwise(1).XData = [h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value];
            h.line_thresh_invSoln_perf_ROC_stepwise(2).XData = [h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value];
            h.line_thresh_invSoln_perf_ROC_std(1).XData = [h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value];
            h.line_thresh_invSoln_perf_ROC_std(2).XData = [h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value];
        end
        
        % text for x-axis point
        xpos = h.slider_3D_image_thresh.Value+(range(h.axes_invSoln_perf_ROC_stepwise.XLim)*.05);
        xpos_add = 0;
        if xpos>range(h.axes_invSoln_perf_ROC_stepwise.XLim)/2  % flip to left side of line
            xpos = h.slider_3D_image_thresh.Value-(range(h.axes_invSoln_perf_ROC_stepwise.XLim)*.35);
            xpos_add = -(range(h.axes_invSoln_perf_ROC_stepwise.XLim)*.1);
        end
        ypos_add = range(h.axes_invSoln_perf_ROC_stepwise.YLim)*.1;
        
        
        %% text for # of Hits & FAs
        ss = find(h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals<=h.slider_3D_image_thresh.Value);
        num_hits = h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_hits(ss(end));
        num_fa = h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_fa(ss(end));
        h.line_thresh_invSoln_perf_ROC_stepwise_txt_hits(1).Position(1) = xpos;
        h.line_thresh_invSoln_perf_ROC_stepwise_txt_hits(1).Position(2) = num_hits - ypos_add;
        h.line_thresh_invSoln_perf_ROC_stepwise_txt_hits(1).String = sprintf('Hits = %.f',num_hits);
        h.line_thresh_invSoln_perf_ROC_stepwise_txt_fa(1).Position(1)   = xpos;
        h.line_thresh_invSoln_perf_ROC_stepwise_txt_fa(1).Position(2)   = num_fa + ypos_add;
        h.line_thresh_invSoln_perf_ROC_stepwise_txt_fa(1).String = sprintf('FA = %.f',num_fa);
        %% text for MCC
        mcc = h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_MCC(ss(end));
        h.line_thresh_invSoln_perf_ROC_stepwise_txt_MCC(1).Position(1)   = xpos + xpos_add;
        h.line_thresh_invSoln_perf_ROC_stepwise_txt_MCC(1).String = sprintf('MCC=%.2f',mcc);
        %% text for Loc Error
        loc_error = nanmean(h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_loc_error(ss(end),:),2);
        if isnan(loc_error); lc_txt = 'Miss'; else; lc_txt = sprintf('%.1f',loc_error); end
        h.line_thresh_invSoln_perf_ROC_std_txt(1).Position(1) = xpos;
        h.line_thresh_invSoln_perf_ROC_std_txt(1).String = sprintf('Avg=%s',lc_txt);
        %% text for wave error
        xdata = h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_act;
        wave_error = nanmean(xdata(ss(end),:),2);
        if isnan(wave_error); act_txt = 'Miss'; else; act_txt = sprintf('%.f',wave_error); end
        h.line_thresh_invSoln_perf_ROC_std_txt_MCC(1).Position(1) = xpos;
        h.line_thresh_invSoln_perf_ROC_std_txt_MCC(1).String = sprintf('Avg=%s',act_txt);
    end
end
