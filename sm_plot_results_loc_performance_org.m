function sm_plot_results_loc_performance(varargin)
global h;

%% Plotting results 
p_clr = [.8 0 .8];
%% Step-Wise plots 
xdata = [h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_hits; h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_fa]; 
if max(max(xdata))<=0; xdata = 1e-8; end
h.axes_invSoln_perf_ROC_stepwise.NextPlot = 'replace';
p1 = plot(h.axes_invSoln_perf_ROC_stepwise, h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals, xdata, 'linewidth',2); axis(h.axes_invSoln_perf_ROC_stepwise,[h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals([1 end]) 0 max(max(xdata))*1.2])
p1(1).Color = [0 .6 1]; p1(2).Color = [1 .6 0]; %legend(h.axes_invSoln_perf_ROC_stepwise,{'Hits' 'FA'}); 
h.axes_invSoln_perf_ROC_stepwise.NextPlot = 'add';
h.line_thresh_invSoln_perf_ROC_stepwise(1) = plot(h.axes_invSoln_perf_ROC_stepwise,[h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value],h.axes_invSoln_perf_ROC_stepwise.YLim,'r--');
% h.line_thresh_invSoln_perf_ROC_stepwise_txt(1) = text(h.axes_invSoln_perf_ROC_stepwise,h.slider_3D_image_thresh.Value*1.05,h.axes_invSoln_perf_ROC_stepwise.YLim(2)*.1,sprintf('%.1e',h.slider_3D_image_thresh.Value),'Color','r');
%% text added for # of Hits & FA 
xpos = h.slider_3D_image_thresh.Value+(range([h.slider_3D_image_thresh.Min h.slider_3D_image_thresh.Max])*.05); 
ss = find(h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals<=h.slider_3D_image_thresh.Value); 
num_hits = h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_hits(ss(end));
h.line_thresh_invSoln_perf_ROC_stepwise_txt_hits(1) = text(h.axes_invSoln_perf_ROC_stepwise,xpos,h.axes_invSoln_perf_ROC_stepwise.YLim(2)*.75,sprintf('Hits = %.f',num_hits),'Color',p1(1).Color);
num_fa = h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_num_fa(ss(end));
h.line_thresh_invSoln_perf_ROC_stepwise_txt_fa(1) = text(h.axes_invSoln_perf_ROC_stepwise,xpos,h.axes_invSoln_perf_ROC_stepwise.YLim(2)*.6,sprintf('FA = %.f',num_fa),'Color',p1(2).Color);

h.axes_invSoln_perf_ROC_stepwise_MCC.NextPlot = 'replace';
plot(h.axes_invSoln_perf_ROC_stepwise_MCC, h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals, h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_MCC,'Color', p_clr,'linewidth',2); axis(h.axes_invSoln_perf_ROC_stepwise_MCC,[h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresh_vals([1 end]) -.01 1.2])
h.axes_invSoln_perf_ROC_stepwise_MCC.NextPlot = 'add';
h.line_thresh_invSoln_perf_ROC_stepwise(2) = plot(h.axes_invSoln_perf_ROC_stepwise_MCC,[h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value],h.axes_invSoln_perf_ROC_stepwise.YLim,'r--');
% h.line_thresh_invSoln_perf_ROC_stepwise_txt(2) = text(h.axes_invSoln_perf_ROC_stepwise_MCC,h.slider_3D_image_thresh.Value*1.05,h.axes_invSoln_perf_ROC_stepwise_MCC.YLim(2)*.9,sprintf('%.1e',h.slider_3D_image_thresh.Value),'Color','r');
h.line_thresh_invSoln_perf_ROC_stepwise_txt_MCC(1) = text(h.axes_invSoln_perf_ROC_stepwise_MCC,h.slider_3D_image_thresh.Value*1.05,h.axes_invSoln_perf_ROC_stepwise_MCC.YLim(2)*.9,sprintf('%.2f',h.inv_soln(h.current_inv_soln).classifier_results.stepwise_thresholded_MCC(ss(end))),'Color',p_clr,'FontSize',8);

%% StdDev Plots 
xdata = [h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_num_hits; h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_num_fa]; 
if max(max(xdata))<=0; xdata = 1e-8; end
h.axes_invSoln_perf_ROC_std.NextPlot = 'replace';
p1 = plot(h.axes_invSoln_perf_ROC_std, h.inv_soln(h.current_inv_soln).classifier_results.std_thresh_vals, xdata, 'linewidth',2); axis(h.axes_invSoln_perf_ROC_std,[h.inv_soln(h.current_inv_soln).classifier_results.std_thresh_vals([1 end]) 0 max(max(xdata))*1.2])
p1(1).Color = [0 .6 1]; p1(2).Color = [1 .6 0]; 
h.axes_invSoln_perf_ROC_std.NextPlot = 'add';
h.line_thresh_invSoln_perf_ROC_std(1) = plot(h.axes_invSoln_perf_ROC_std,[h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value]/h.inv_soln(h.current_inv_soln).classifier_results.stddev_img,h.axes_invSoln_perf_ROC_stepwise.YLim,'r--');
std_val = h.slider_3D_image_thresh.Value/h.inv_soln(h.current_inv_soln).classifier_results.stddev_img;
h.line_thresh_invSoln_perf_ROC_std_txt(1) = text(h.axes_invSoln_perf_ROC_std,std_val+.5,h.axes_invSoln_perf_ROC_stepwise.YLim(2)*.925,sprintf('StdDev = %.1f',std_val),'Color','r');
% legend(h.axes_invSoln_perf_ROC_std,p1,{'Hits' 'FA'}); 

h.axes_invSoln_perf_ROC_std_MCC.NextPlot = 'replace';
plot(h.axes_invSoln_perf_ROC_std_MCC, h.inv_soln(h.current_inv_soln).classifier_results.std_thresh_vals, h.inv_soln(h.current_inv_soln).classifier_results.std_thresholded_MCC, 'Color', p_clr,'linewidth',2); axis(h.axes_invSoln_perf_ROC_std_MCC,[h.inv_soln(h.current_inv_soln).classifier_results.std_thresh_vals([1 end]) -.01 1.2])
h.axes_invSoln_perf_ROC_std_MCC.NextPlot = 'add';
h.line_thresh_invSoln_perf_ROC_std(2) = plot(h.axes_invSoln_perf_ROC_std_MCC,[h.slider_3D_image_thresh.Value h.slider_3D_image_thresh.Value]/h.inv_soln(h.current_inv_soln).classifier_results.stddev_img,h.axes_invSoln_perf_ROC_stepwise.YLim,'r--');
% h.line_thresh_invSoln_perf_ROC_std_txt(2) = text(h.axes_invSoln_perf_ROC_std_MCC,std_val+.5,h.axes_invSoln_perf_ROC_std_MCC.YLim(2)*.925,sprintf('%.1f',std_val),'Color','r');


title(h.axes_invSoln_perf_ROC_stepwise,sprintf('Hits & FA (Step-Wise)')); ylabel(h.axes_invSoln_perf_ROC_stepwise,'Number Sources'); xlabel(h.axes_invSoln_perf_ROC_stepwise,'Threshold Value');
title(h.axes_invSoln_perf_ROC_std,sprintf('Hits & FA (StdDev)')); ylabel(h.axes_invSoln_perf_ROC_std,'Number Sources'); xlabel(h.axes_invSoln_perf_ROC_std,'Standard Deviation');
title(h.axes_invSoln_perf_ROC_stepwise_MCC,sprintf('MCC: Step-Wise Thresholded')); ylabel(h.axes_invSoln_perf_ROC_stepwise_MCC,'MCC'); xlabel(h.axes_invSoln_perf_ROC_stepwise_MCC,'Threshold Value');
title(h.axes_invSoln_perf_ROC_std_MCC,sprintf('ROI: StdDev Thresholded')); ylabel(h.axes_invSoln_perf_ROC_std_MCC,'MCC'); xlabel(h.axes_invSoln_perf_ROC_std_MCC,'Standard Deviation');


