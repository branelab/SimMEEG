function sm_get_inv_tfr_point(varargin)
global h

rotate3d off;

if ~isempty(h.current_inv_tfr_time_point)
    
%     % setting axes and getting input from user for selecting time-frequency point on the graph
%     if h.radio_inv_plot_peak_tfr_connectivity.Value==1 && h.radio_inv_plot_true_tfr_connectivity.Value==1
%         axes(h.axes_inv_soln_tfr);
%     elseif h.radio_inv_plot_peak_tfr_connectivity.Value==0 && h.radio_inv_plot_true_tfr_connectivity.Value==1
%         axes(h.axes_inv_soln_tfr);
%     else
%         return
%     end
    
    if strcmp(varargin{1}.Type,'surface')==1  || strcmp(varargin{1}.Type,'line')==1  % user selected time-freq point on TFR plot
        %     [x,y,btn]=ginput(1);
        xyz = get (varargin{1}.Parent, 'CurrentPoint');
        x=xyz(1,1);
        y=xyz(1,2);
    else                                        % called function from another function to plot white time-freq lines
        if isfield(h,'current_inv_tfr_time_point') && isfield(h,'current_inv_tfr_freq_point')
        x = h.current_inv_tfr_time_point;
        y = h.current_inv_tfr_freq_point;
        else
            x = 1;
            y = 1;
        end
    end
    
    
    % staying within plot boundaries
    if h.radio_inv_plot_peak_tfr_connectivity.Value==1 
        if x< h.axes_inv_soln_tfr.XLim(1); x=h.axes_inv_soln_tfr.XLim(1); elseif x>h.axes_inv_soln_tfr.XLim(2); x=h.axes_inv_soln_tfr.XLim(2); end
        if y< h.axes_inv_soln_tfr.YLim(1); y=h.axes_inv_soln_tfr.YLim(1); elseif y>h.axes_inv_soln_tfr.YLim(2); y=h.axes_inv_soln_tfr.YLim(2); end
        %% current time point
        ss=find(h.sim_data.cfg.study.lat_sim<=x);
        if isempty(ss); ss=1; end
        h.current_inv_tfr_time_samp = ss(end);
        h.current_inv_tfr_time_point = h.sim_data.cfg.study.lat_sim(h.current_inv_tfr_time_samp);
        ss=find(h.inv_soln(h.current_inv_soln).TFR_results.TFR_freqs<=y);
        if isempty(ss); ss=1; end
        h.current_inv_tfr_freq_samp = ss(end);
        h.current_inv_tfr_freq_point = h.inv_soln(h.current_inv_soln).TFR_results.TFR_freqs(h.current_inv_tfr_freq_samp);
%                 display(h.current_inv_tfr_freq_point)
    elseif h.radio_inv_plot_true_tfr_connectivity.Value==1
        if x< h.axes_true_source_tfr.XLim(1); x=h.axes_true_source_tfr.XLim(1); elseif x>h.axes_true_source_tfr.XLim(2); x=h.axes_true_source_tfr.XLim(2); end
        if y< h.axes_true_source_tfr.YLim(1); y=h.axes_true_source_tfr.YLim(1); elseif y>h.axes_true_source_tfr.YLim(2); y=h.axes_true_source_tfr.YLim(2); end
        %% current time point
        ss=find(h.cfg.study.lat_sim<=x);
        h.current_inv_tfr_time_samp = ss(end);
        h.current_inv_tfr_time_point = h.cfg.study.lat_sim(h.current_inv_tfr_time_samp);
        ss=find(h.sim_data.cfg.source.TFR_results.TFR_freqs<=y);
        if isempty(ss)
            h.current_inv_tfr_freq_samp = 1;
        else 
            h.current_inv_tfr_freq_samp = ss(end);
        end
        h.current_inv_tfr_freq_point = h.sim_data.cfg.source.TFR_results.TFR_freqs(h.current_inv_tfr_freq_samp);
    end
    
    %% plot tft time-freq point as cross hairs on TFR plot
    
    % Peak Sources
    if h.radio_inv_plot_peak_tfr_connectivity.Value==1 
        if isfield(h,'current_inv_tfr_point_lines'); if any(isvalid(h.current_inv_tfr_point_lines)); delete(h.current_inv_tfr_point_lines); end; end
        h.current_inv_tfr_point_lines = plot3(h.axes_inv_soln_tfr,[x x; h.axes_inv_soln_tfr.XLim]',[h.axes_inv_soln_tfr.YLim; y y]',repmat(h.axes_inv_soln_tfr.CLim(2)*10,2,2),'--','Color',[1 1 1]*1,'linewidth',1);
        for a=1:2; h.current_inv_tfr_point_lines(a).ButtonDownFcn = @sm_get_inv_tfr_point; end
        
        try
        ss = find( h.current_inv_tfr_plot.XData <= h.current_inv_tfr_time_point); ss=ss(end);
        fs = find( h.current_inv_tfr_plot.YData <= h.current_inv_tfr_freq_point); fs=fs(end);
        pt_val = h.current_inv_tfr_plot.CData(fs,ss);
        catch
            ss = 1; fs=1; pt_val=1;
        end
        if length(h.listbox_plv_contrasts.Value)>1
            h.axes_inv_soln_tfr.Title.String = sprintf('Averaged Peak Hit Sources: Time = %.3f Freq = %.1f Val = %.2f',h.current_inv_tfr_time_point,h.current_inv_tfr_freq_point,pt_val);
        else
            h.axes_inv_soln_tfr.Title.String = sprintf('Peak Hit Sources: Time = %.3f Freq = %.1f Val = %.2f',h.current_inv_tfr_time_point,h.current_inv_tfr_freq_point,pt_val);
        end
          h.axes_inv_soln_tfr.Title.FontWeight = 'normal'; h.axes_inv_soln_tfr.Title.FontSize = 9; h.axes_inv_soln_tfr.Title.HorizontalAlignment='left';
        h.axes_inv_soln_tfr.Title.Color = 'k';% h.current_inv_tfr_point_lines.Color;
        h.axes_inv_soln_tfr.Title.Position(1:2) = [h.axes_inv_soln_tfr.XLim(1) h.axes_inv_soln_tfr.YLim(2)*1.05];
    end
    
    if h.radio_inv_plot_true_tfr_connectivity.Value==1 && isvalid(h.current_true_tfr_plot)
        % True Sources
        if isfield(h,'current_true_tfr_point_lines'); if any(isvalid(h.current_true_tfr_point_lines)); delete(h.current_true_tfr_point_lines); end; end
        h.current_true_tfr_point_lines = plot3(h.axes_true_source_tfr,[x x; h.axes_true_source_tfr.XLim]',[h.axes_true_source_tfr.YLim; y y]',repmat(h.axes_true_source_tfr.CLim(2)*10,2,2),'--','Color',[1 1 1]*1,'linewidth',1);
        for a=1:2; h.current_true_tfr_point_lines(a).ButtonDownFcn = @sm_get_inv_tfr_point; end
        
        ss = find( h.current_true_tfr_plot.XData <= h.current_inv_tfr_time_point); ss=ss(end);
        fs = find( h.current_true_tfr_plot.YData <= h.current_inv_tfr_freq_point); fs=fs(end);
        pt_val = h.current_true_tfr_plot.CData(fs,ss);
        
        if length(h.listbox_plv_contrasts.Value)>1
            h.axes_true_source_tfr.Title.String = sprintf('Averaged True Sources: Time = %.3f Freq = %.1f Val = %.2f',h.current_inv_tfr_time_point,h.current_inv_tfr_freq_point,pt_val);
        else
            h.axes_true_source_tfr.Title.String = sprintf('True Sources: Time = %.3f Freq = %.1f Val = %.2f',h.current_inv_tfr_time_point,h.current_inv_tfr_freq_point,pt_val);
        end
        h.axes_true_source_tfr.Title.FontWeight = 'normal'; h.axes_true_source_tfr.Title.FontSize = 9; h.axes_true_source_tfr.Title.HorizontalAlignment='left';
        h.axes_true_source_tfr.Title.Color = 'k';% h.current_inv_tfr_point_lines.Color;
        h.axes_true_source_tfr.Title.Position(1:2) = [h.axes_true_source_tfr.XLim(1) h.axes_true_source_tfr.YLim(2)*1.05];
    end
    
    
    %% plot connectivity graph
    sm_plot_connectivity_graph;
    disableDefaultInteractivity(h.axes_true_source_tfr); h.axes_true_source_tfr.Toolbar.Visible='off';
    disableDefaultInteractivity(h.axes_inv_soln_tfr); h.axes_inv_soln_tfr.Toolbar.Visible='off';
else
    % initially set time-freq point at 0 ms and 1st frequency
    h.current_inv_tfr_time_point=[1 1];
    
end



