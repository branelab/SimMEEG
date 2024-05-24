function sm_plot_results_loc_performance(varargin)
global h;

if ~isempty(h.inv_soln(h.current_inv_soln).classifier_results)
try
h = rmfield(h,'line_thresh_invSoln_perf_ROC_stepwise');
h = rmfield(h,'line_thresh_invSoln_perf_ROC_std');
catch
end
%% Plotting results 
p_clr = [.8 0 .8];
%% Hits/FAs and MCC
xdata = [h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_hits; h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_fa]; 
if max(max(xdata))<=0; xdata = 1e-8; end
h.axes_invSoln_perf_ROC_stepwise.NextPlot = 'replace';
p1 = plot(h.axes_invSoln_perf_ROC_stepwise, h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals, xdata, 'linewidth',2); axis(h.axes_invSoln_perf_ROC_stepwise,[h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals([1 end]) 0 max(max(xdata))*1.2])
p1(1).Color = [0 .6 1]; p1(2).Color = [1 .6 0]; %legend(h.axes_invSoln_perf_ROC_stepwise,{'Hits' 'FA'}); 
h.axes_invSoln_perf_ROC_stepwise.NextPlot = 'add';
h.line_thresh_invSoln_perf_ROC_stepwise = plot(h.axes_invSoln_perf_ROC_stepwise,[h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value],h.axes_invSoln_perf_ROC_stepwise.YLim,'k--');
% h.line_thresh_invSoln_perf_ROC_stepwise_txt = text(h.axes_invSoln_perf_ROC_stepwise,h.slider_3D_image_thresh.Value*1.05,h.axes_invSoln_perf_ROC_stepwise.YLim(2)*.1,sprintf('%.1e',h.slider_3D_image_thresh.Value),'Color','r');
% text added for # of Hits & FA 
xpos = h.slider_3D_image_thresh.Value+(range([h.slider_3D_image_thresh.Min h.slider_3D_image_thresh.Max])*.05); 
ypos_add = double(range(h.axes_invSoln_perf_ROC_stepwise.YLim)*.1);
ss = find(h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals<=h.slider_3D_image_thresh.Value); 
num_hits = h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_hits(ss(end));
h.line_thresh_invSoln_perf_ROC_stepwise_txt_hits = text(h.axes_invSoln_perf_ROC_stepwise,xpos,xdata(1,1)+ypos_add,sprintf('Hits = %.f',num_hits),'Color',p1(1).Color,'FontSize',8);
num_fa = h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_fa(ss(end));
h.line_thresh_invSoln_perf_ROC_stepwise_txt_fa = text(h.axes_invSoln_perf_ROC_stepwise,xpos,xdata(2,1)+ypos_add,sprintf('FA = %.f',num_fa),'Color',p1(2).Color,'FontSize',8);
axis(h.axes_invSoln_perf_ROC_stepwise,[0 h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals(end)*1.1 0 max(max(xdata))*1.2]);

%% MCC
h.axes_invSoln_perf_ROC_stepwise_MCC.NextPlot = 'replace';
plot(h.axes_invSoln_perf_ROC_stepwise_MCC, h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals, h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_MCC,'Color', p_clr,'linewidth',2); axis(h.axes_invSoln_perf_ROC_stepwise_MCC,[h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals([1 end]) -.01 1.2])
h.axes_invSoln_perf_ROC_stepwise_MCC.NextPlot = 'add';
h.line_thresh_invSoln_perf_ROC_stepwise(2) = plot(h.axes_invSoln_perf_ROC_stepwise_MCC,[h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value],h.axes_invSoln_perf_ROC_stepwise.YLim,'k--');
% h.line_thresh_invSoln_perf_ROC_stepwise_txt(2) = text(h.axes_invSoln_perf_ROC_stepwise_MCC,h.slider_3D_image_thresh.Value*1.05,h.axes_invSoln_perf_ROC_stepwise_MCC.YLim(2)*.9,sprintf('%.1e',h.slider_3D_image_thresh.Value),'Color','r');
h.line_thresh_invSoln_perf_ROC_stepwise_txt_MCC = text(h.axes_invSoln_perf_ROC_stepwise_MCC,h.slider_3D_image_thresh.Value*1.05,double(h.axes_invSoln_perf_ROC_stepwise_MCC.YLim(2))*.9,sprintf('%.2f',h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_MCC(ss(end))),'Color',p_clr,'FontSize',8);
axis(h.axes_invSoln_perf_ROC_stepwise_MCC,[0 h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals(end)*1.1 0 1.2]);

%% Loc Error
xdata = [h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_loc_error]; 
if max(max(xdata))<=0; xdata = 1e-8; end
h.axes_invSoln_perf_ROC_std.NextPlot = 'replace';
p1 = plot(h.axes_invSoln_perf_ROC_std, h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals, xdata, 'linewidth',2); 
yscale = [0 max(max(xdata))*1.2]; if any(isnan(yscale)); yscale = [0 20]; end
axis(h.axes_invSoln_perf_ROC_std,[h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals([1 end]) yscale])

p1(1).Color = h.src_clr(1,:); p1(2).Color = h.src_clr(2,:); p1(3).Color = h.src_clr(3,:); 
h.axes_invSoln_perf_ROC_std.NextPlot = 'add';
h.line_thresh_invSoln_perf_ROC_std = plot(h.axes_invSoln_perf_ROC_std,[h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value],h.axes_invSoln_perf_ROC_std.YLim,'k--');
loc_error = nanmean(h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_loc_error(ss(end),:),2);
if isnan(loc_error); lc_txt = 'Miss'; else; lc_txt = sprintf('%.1f mm',loc_error); end
h.line_thresh_invSoln_perf_ROC_std_txt = text(h.axes_invSoln_perf_ROC_std,h.slider_3D_image_thresh.Value*1.05,max(max(xdata))+(range(h.axes_invSoln_perf_ROC_std.YLim)*.1),sprintf('%s',lc_txt),'Color',[0 0 0],'FontSize',8);
% legend(h.axes_invSoln_perf_ROC_std,p1,{'Hits' 'FA'}); 

%% Wave Error
xdata = h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_wave_error_mse_norm_act; 
if max(max(xdata))<=0; xdata = 1e-8; end
h.axes_invSoln_perf_ROC_std_MCC.NextPlot = 'replace';
p1 = plot(h.axes_invSoln_perf_ROC_std_MCC, h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals, xdata, 'linewidth',2); axis(h.axes_invSoln_perf_ROC_std_MCC,[h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals([1 end]) 0 max(max(xdata))*1.2])
p1(1).Color = h.src_clr(1,:); p1(2).Color = h.src_clr(2,:); p1(3).Color = h.src_clr(3,:); 
h.axes_invSoln_perf_ROC_std_MCC.NextPlot = 'add';
h.line_thresh_invSoln_perf_ROC_std(2) = plot(h.axes_invSoln_perf_ROC_std_MCC,[h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value],h.axes_invSoln_perf_ROC_std_MCC.YLim,'k--');
wave_error = nanmean(xdata(ss(end),:),2);
if isnan(wave_error); act_txt = 'Miss'; else; act_txt = sprintf('%.f',wave_error); end
h.line_thresh_invSoln_perf_ROC_std_txt_MCC = text(h.axes_invSoln_perf_ROC_std_MCC,h.slider_3D_image_thresh.Value*1.05,max(max(xdata))+(range(h.axes_invSoln_perf_ROC_std_MCC.YLim)*.1),sprintf('Avg=%s',wave_error),'Color',[0 0 0],'FontSize',8);
% legend(h.axes_invSoln_perf_ROC_std_MCC,p1,{'Hits' 'FA'}); 
axis(h.axes_invSoln_perf_ROC_std_MCC,[0 h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals(end)*1.1 0 max(max(xdata))*1.2]);


%% Titles and scaling
title(h.axes_invSoln_perf_ROC_stepwise,sprintf('Hits & FA (Step-Wise)')); ylabel(h.axes_invSoln_perf_ROC_stepwise,'Number Sources'); xlabel(h.axes_invSoln_perf_ROC_stepwise,'Threshold Value');
title(h.axes_invSoln_perf_ROC_std,sprintf('Localization Error')); ylabel(h.axes_invSoln_perf_ROC_std,'Error (mm)'); xlabel(h.axes_invSoln_perf_ROC_std,'Threshold Value');
title(h.axes_invSoln_perf_ROC_stepwise_MCC,sprintf('MCC')); ylabel(h.axes_invSoln_perf_ROC_stepwise_MCC,'MCC'); xlabel(h.axes_invSoln_perf_ROC_stepwise_MCC,'Threshold Value');
title(h.axes_invSoln_perf_ROC_std_MCC,sprintf('Wave Error (MSE)')); ylabel(h.axes_invSoln_perf_ROC_std_MCC,'MSE'); xlabel(h.axes_invSoln_perf_ROC_std_MCC,'Threshold Value');

end


