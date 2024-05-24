function sm_plot_ROI_FC_metrics(varargin)
global h;
h.axes_invSoln_FC_ROC_stepwise.clo; h.axes_invSoln_FC_ROC_stepwise.NextPlot = 'add';
h.axes_invSoln_FC_ROC_surg_std.clo; h.axes_invSoln_FC_ROC_surg_std.NextPlot = 'add';
h.axes_invSoln_TFR_ROC_stepwise.clo; h.axes_invSoln_TFR_ROC_stepwise.NextPlot = 'add';
p1 = []; p2 = []; % initiating plot handles
PLV_wave_mse = nan(1,3); PLV_wave_abs =  nan(1,3); PLI_wave_mse =  nan(1,3); PLI_wave_abs =  nan(1,3);


%% Initialize variable
TFR_wave_mse = []; 
TFR_wave_mse_evk = [];
TFR_wave_mse_ind = []; 
PLV_wave_mse = []; 
PLI_wave_mse = []; 

%% %%%%%% Plot step-wise thresholded TFR ROI ROC data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TFR TotPwr
if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'TotPwr_ROI_FC_stepwise_metrics')
    TFR_wave_mse = h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.ROI_FC_wave_error_mse;
    TFR_wave_abs = h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.ROI_FC_wave_error_abs;
    if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR)
%         x = h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR;
%         y = h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.fitted_ROI_TPR;
        x = [h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.ROI_perf.FPR];
        y = [h.inv_soln(h.current_inv_soln).TFR_results.TotPwr_ROI_FC_stepwise_metrics.ROI_perf.TPR];
        p1 = plot(h.axes_invSoln_TFR_ROC_stepwise,x,y,'color',h.tfr_clr(1,:),'linewidth',2);
        auc = area(h.axes_invSoln_TFR_ROC_stepwise,x,y); auc.FaceColor = p1.Color; auc.FaceAlpha = .15;
    end
end
%% TFR IndPwr
if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'IndPwr_ROI_FC_stepwise_metrics')
    TFR_wave_mse_ind = h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.ROI_FC_wave_error_mse;
    TFR_wave_abs_ind = h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.ROI_FC_wave_error_abs;
    if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR)
%         x = h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR;
%         y = h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.fitted_ROI_TPR;
        x = [h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.ROI_perf.FPR];
        y = [h.inv_soln(h.current_inv_soln).TFR_results.IndPwr_ROI_FC_stepwise_metrics.ROI_perf.TPR];
        p2 = plot(h.axes_invSoln_TFR_ROC_stepwise,x,y,'color',h.tfr_clr(2,:),'linewidth',2);
        auc = area(h.axes_invSoln_TFR_ROC_stepwise,x,y); auc.FaceColor = p2.Color; auc.FaceAlpha = .15;
    end
end
%% TFR EvkPwr
if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'EvkPwr_ROI_FC_stepwise_metrics')
    TFR_wave_mse_evk = h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.ROI_FC_wave_error_mse;
    TFR_wave_abs_evk = h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.ROI_FC_wave_error_abs;
    if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR)
%         x = h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.fitted_ROI_FPR;
%         y = h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.fitted_ROI_TPR;
        x = [h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.ROI_perf.FPR];
        y = [h.inv_soln(h.current_inv_soln).TFR_results.EvkPwr_ROI_FC_stepwise_metrics.ROI_perf.TPR];
        p3 = plot(h.axes_invSoln_TFR_ROC_stepwise,x,y,'color',h.tfr_clr(3,:),'linewidth',2);
        auc = area(h.axes_invSoln_TFR_ROC_stepwise,x,y); auc.FaceColor = p3.Color; auc.FaceAlpha = .15;
    end
end
%% extras
title(h.axes_invSoln_TFR_ROC_stepwise,sprintf('Step-Wise Thresholded'));
axis(h.axes_invSoln_TFR_ROC_stepwise,[0 1 0 1]); box on; 
plot(h.axes_invSoln_TFR_ROC_stepwise,[0 1],[0 1],'k--'); 
try
    legend(h.axes_invSoln_TFR_ROC_stepwise,[p1 p2 p3],{'TotPwr' 'IndPwr' 'EvkPwr'},'Location','southeast');
catch
end
    ylabel(h.axes_invSoln_TFR_ROC_stepwise,'True-Positive Rate'); xlabel(h.axes_invSoln_TFR_ROC_stepwise,'False-Positive Rate'); 

%% %%%%% Bar plot of TFR wave Errors
%% TFR wave error MSE
min_max = [0 20];
min_max2 = [min([TFR_wave_mse TFR_wave_mse_ind TFR_wave_mse_evk]) max([TFR_wave_mse TFR_wave_mse_ind TFR_wave_mse_evk])*1.2];
if ~isempty(min_max2)
    if min_max2(2)>min_max(2); min_max(2) = ceil(min_max2(2)); end
    b = bar(h.axes_invSoln_TFR_ROC_wave_error_mse,[TFR_wave_mse; TFR_wave_mse_ind; TFR_wave_mse_evk]);
    for v=1:length(b)
        b(v).FaceColor = h.cfg.source.src_clr(v,:);
        
        xtips2 = b(v).XEndPoints;
        ytips2 = b(v).YEndPoints;
        labels2 = string(round(b(v).YData,0));
        text(h.axes_invSoln_TFR_ROC_wave_error_mse,xtips2,ytips2,labels2,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','FontSize',8)
    end
    h.axes_invSoln_TFR_ROC_wave_error_mse.XTickLabel = {'TotPwr' 'IndPwr' 'EvkPwr'};
    h.axes_invSoln_TFR_ROC_wave_error_mse.YLim(2) = min_max(2);
    title(h.axes_invSoln_TFR_ROC_wave_error_mse,sprintf('ROI: TFR wave Error (MSE)'));
    ylabel(h.axes_invSoln_TFR_ROC_wave_error_mse,'TFR Wave Error');
    %% TFR wave error Absolute
    min_max = [0 20];
    min_max2 = [min([TFR_wave_abs TFR_wave_abs_ind TFR_wave_abs_evk]) max([TFR_wave_abs TFR_wave_abs_ind TFR_wave_abs_evk])*1.2];
    if min_max2(2)>min_max(2); min_max(2) = ceil(min_max2(2)); end
    b = bar(h.axes_invSoln_TFR_ROC_wave_error_abs,[TFR_wave_abs; TFR_wave_abs_ind; TFR_wave_abs_evk]);
    for v=1:length(b)
        b(v).FaceColor = h.cfg.source.src_clr(v,:);
        
        xtips2 = b(v).XEndPoints;
        ytips2 = b(v).YEndPoints;
        labels2 = string(round(b(v).YData,0));
        text(h.axes_invSoln_TFR_ROC_wave_error_abs,xtips2,ytips2,labels2,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','FontSize',8)
    end
    h.axes_invSoln_TFR_ROC_wave_error_abs.XTickLabel = {'TotPwr' 'IndPwr' 'EvkPwr'};
    h.axes_invSoln_TFR_ROC_wave_error_abs.YLim(2) = min_max(2);
    title(h.axes_invSoln_TFR_ROC_wave_error_abs,sprintf('ROI: TFR wave Error (Absolute)'));
    ylabel(h.axes_invSoln_TFR_ROC_wave_error_abs,'TFR Wave Error');
else
end



%% %%%%%% Plot step-wise thresholded FC ROI ROC data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLV
if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'PLV_ROI_FC_stepwise_metrics')
    PLV_wave_mse = h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.ROI_FC_wave_error_mse;
    PLV_wave_abs = h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.ROI_FC_wave_error_abs;

    if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.fitted_ROI_FPR)
        %     x = h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.fitted_ROI_FPR;
        %     y = h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.fitted_ROI_TPR;
        x = [h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.ROI_perf.FPR];
        y = [h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_stepwise_metrics.ROI_perf.TPR];
        p1 = plot(h.axes_invSoln_FC_ROC_stepwise,x,y,'color',h.tab_PLV.ForegroundColor,'linewidth',2);
        auc = area(h.axes_invSoln_FC_ROC_stepwise,x,y); auc.FaceColor = p1.Color; auc.FaceAlpha = .15;
    end
end
%% PLI
if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'PLI_ROI_FC_stepwise_metrics')
     PLI_wave_mse = h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.ROI_FC_wave_error_mse;
     PLI_wave_abs = h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.ROI_FC_wave_error_abs;
   
if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.fitted_ROI_FPR)
%     x = h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.fitted_ROI_FPR;
%     y = h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.fitted_ROI_TPR;
    x = [h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.ROI_perf.FPR];
    y = [h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_stepwise_metrics.ROI_perf.TPR];
    p2 = plot(h.axes_invSoln_FC_ROC_stepwise,x,y,'color',h.tab_PLI.ForegroundColor,'linewidth',2);
    auc = area(h.axes_invSoln_FC_ROC_stepwise,x,y); auc.FaceColor = p2.Color; auc.FaceAlpha = .15;
end
end
%% extras
title(h.axes_invSoln_FC_ROC_stepwise,sprintf('Step-Wise Thresholded'));
axis(h.axes_invSoln_FC_ROC_stepwise,[0 1 0 1]); box on; 
plot(h.axes_invSoln_FC_ROC_stepwise,[0 1],[0 1],'k--'); 
legend(h.axes_invSoln_FC_ROC_stepwise,[p1 p2],{'PLV' 'PLI'},'Location','southeast');
ylabel(h.axes_invSoln_FC_ROC_stepwise,'True-Positive Rate'); xlabel(h.axes_invSoln_FC_ROC_stepwise,'False-Positive Rate'); 
%% %%%%%% Plot surg_std thresholded FC ROI ROC data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLV
if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'PLV_ROI_FC_surg_std_metrics')
    if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics)
        PLV_wave_mse = h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.ROI_FC_wave_error_mse;
        PLV_wave_abs = h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.ROI_FC_wave_error_abs;
    if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.fitted_ROI_FPR)
%     x = h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.fitted_ROI_FPR;
%     y = h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.fitted_ROI_TPR;
    x = [h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.ROI_perf.FPR];
    y = [h.inv_soln(h.current_inv_soln).TFR_results.PLV_ROI_FC_surg_std_metrics.ROI_perf.TPR];
    p1 = plot(h.axes_invSoln_FC_ROC_surg_std,x,y,'color',h.tab_PLV.ForegroundColor,'linewidth',2);
    auc = area(h.axes_invSoln_FC_ROC_surg_std,x,y); auc.FaceColor = p1.Color; auc.FaceAlpha = .15;
    end
    else
        PLV_wave_mse = []; PLV_wave_abs =[]; x=[];y=[]; p1=[];auc=[];
    end

end
%% PLI
if isfield(h.inv_soln(h.current_inv_soln).TFR_results,'PLI_ROI_FC_surg_std_metrics')
    if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics)
        PLI_wave_mse = h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.ROI_FC_wave_error_mse;
        PLI_wave_abs = h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.ROI_FC_wave_error_abs;
        
        if ~isempty(h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.fitted_ROI_FPR)
%             x = h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.fitted_ROI_FPR;
%             y = h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.fitted_ROI_TPR;
            x = [h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.ROI_perf.FPR];
            y = [h.inv_soln(h.current_inv_soln).TFR_results.PLI_ROI_FC_surg_std_metrics.ROI_perf.TPR];
            p2 = plot(h.axes_invSoln_FC_ROC_surg_std,x,y,'color',h.tab_PLI.ForegroundColor,'linewidth',2);
            auc = area(h.axes_invSoln_FC_ROC_surg_std,x,y); auc.FaceColor = p2.Color; auc.FaceAlpha = .15;
        end
    else
        PLV_wave_mse = []; PLV_wave_abs =[]; x=[];y=[]; p1=[];auc=[];
    end
end
%% extras
title(h.axes_invSoln_FC_ROC_surg_std,sprintf('Surrogate StdDev Thresholded'));
axis(h.axes_invSoln_FC_ROC_surg_std,[0 1 0 1]); box on; 
plot(h.axes_invSoln_FC_ROC_surg_std,[0 1],[0 1],'k--'); 
legend(h.axes_invSoln_FC_ROC_surg_std,[p1 p2],{'PLV' 'PLI'}','Location','southeast');
ylabel(h.axes_invSoln_FC_ROC_surg_std,'True-Positive Rate'); xlabel(h.axes_invSoln_FC_ROC_surg_std,'False-Positive Rate');

%% %%%%% Bar plot of FC wave Errors
%% FC wave error MSE
min_max = [0 .05];
min_max2 = [min([PLV_wave_mse PLI_wave_mse]) max([PLV_wave_mse PLI_wave_mse])*1.2];
if ~isempty(min_max2)
    if min_max2(2)>min_max(2); min_max(2) = min_max2(2); end
    b = bar(h.axes_invSoln_FC_ROC_wave_error_mse,[PLV_wave_mse; PLI_wave_mse]);
    for v=1:length(b)
        b(v).FaceColor = 'none'; %h.plv_clr(v,:);
        
        xtips2 = b(v).XEndPoints;
        ytips2 = b(v).YEndPoints;
        labels2 = string(round(b(v).YData,2));
        text(h.axes_invSoln_FC_ROC_wave_error_mse,xtips2,ytips2,labels2,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','FontSize',8)
    end
    h.axes_invSoln_FC_ROC_wave_error_mse.XTickLabel = {'PLV' 'PLI'};
    h.axes_invSoln_FC_ROC_wave_error_mse.YLim(2) = min_max(2);
    title(h.axes_invSoln_FC_ROC_wave_error_mse,sprintf('ROI: FC wave Error (MSE)'));
    ylabel(h.axes_invSoln_FC_ROC_wave_error_mse,'FC Wave Error');
    %% FC wave error Absolute
    min_max = [0 .25];
    min_max2 = [min([PLV_wave_abs PLI_wave_abs]) max([PLV_wave_abs PLI_wave_abs])*1.2];
    if min_max2(2)>min_max(2); min_max(2) = min_max2(2); end
    
    b = bar(h.axes_invSoln_FC_ROC_wave_error_abs,[PLV_wave_abs; PLI_wave_abs]);
    for v=1:length(b)
        b(v).FaceColor = 'none'; %h.plv_clr(v,:);
        
        xtips2 = b(v).XEndPoints;
        ytips2 = b(v).YEndPoints;
        labels2 = string(round(b(v).YData,2));
        text(h.axes_invSoln_FC_ROC_wave_error_abs,xtips2,ytips2,labels2,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom','FontSize',8)
    end
    h.axes_invSoln_FC_ROC_wave_error_abs.XTickLabel = {'PLV' 'PLI'};
    h.axes_invSoln_FC_ROC_wave_error_abs.YLim(2) = min_max(2);
    title(h.axes_invSoln_FC_ROC_wave_error_abs,sprintf('ROI: FC wave Error (Absolute)'));
    ylabel(h.axes_invSoln_FC_ROC_wave_error_abs,'FC Wave Error');
else
end



