function listbox_peaks_found_Callback(varargin)
global h

%% plot inv_soln TFR waves

%% plot inv_soln time-domain waves

if h.radio_3D_peak_locs.Value==1    % peaks need to be shown on map to highlight them
mrk_size2 = 50;
    try
%         if length(h.listbox_inv_solns.Value)==1
            %                 ln_clr = h.map3D_peak_locs(1).CData; %bsxfun(@mtimes,ones(length(h.current_3D_peak_idx),3),[.9 .6 .3]*.75);  %lines(length(h4.current_3D_peak_idx));
            %         ln_clr(1:3,:) = h.src_clr;
            for v=1:size(h.current_peak_swf,2)
                h.current_inv_swf_plots(v).Color = h.ln_clr(v,:);
                h.current_inv_swf_plots_fft(v).Color = h.ln_clr(v,:);
                %     if v>3
                h.current_inv_swf_plots(v).Color(4) = h.false_positive_lineAlpha;
                h.current_inv_swf_plots_fft(v).Color(4) = h.false_positive_lineAlpha;
                %     end
            end
            
            %% Updated locations for selected voxels
            sel_idx = h.current_3D_peak_idx(h.listbox_peaks_found.Value);
            %                 sel_idx = sel_idx(~isnan(sel_idx));
            if all(isnan(sel_idx))
                sel_idx = nan; sel_locs = [];
            else
                sel_locs = nan(length(sel_idx),3);
                sel_locs(~isnan(sel_idx),:) = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos(sel_idx(~isnan(sel_idx)),:);
            end
            
            if any(~isnan(sel_idx))
                
                if isfield(h,'current_peak_selected')
                    if isvalid(h.current_peak_selected)
                        h.current_peak_selected.XData = sel_locs(:,1); h.current_peak_selected.YData = sel_locs(:,2); h.current_peak_selected.ZData = sel_locs(:,3);
                        %         h.current_inv_swf_plots(h.listbox_peaks_found.Value).Color(1:3) = [1 0 1];
                        for v=h.listbox_peaks_found.Value
                            h.current_inv_swf_plots(v).Color(4) = 1;
                            h.current_inv_swf_plots_fft(v).Color(4) = 1; h.current_inv_swf_plots_fft(v).LineWidth = 3;
                            uistack(h.current_inv_swf_plots(v),'top'); uistack(h.current_inv_swf_plots_fft(v),'top');
                        end
                    else
                        h.current_peak_selected = scatter3(h.axes_3D_images,sel_locs(:,1),sel_locs(:,2),sel_locs(:,3),'o','MarkerEdgeColor',[1 0 1]*0,'SizeData',mrk_size2,'linewidth',2);
                    end
                else
                    h.current_peak_selected = scatter3(h.axes_3D_images,sel_locs(:,1),sel_locs(:,2),sel_locs(:,3),'o','MarkerEdgeColor',[1 0 1]*0,'SizeData',mrk_size2,'linewidth',2);
                end
            else
                delete(h.current_peak_selected);
            end
%         end
        
        %% Update legend
        axes(h.axes_source_waves); legend off;
        if h.radio_3D_true_locs.Value==1 && h.radio_3D_peak_locs.Value==0 && ~isempty(h.current_3D_peak_idx) % only true source waves
            xdata = h.norm_true_swf;
            lg_txt=[]; for v=1:size(xdata,2); lg_txt{v} = sprintf('True Source %.f',v); end
            legend(h.axes_source_waves,h.current_true_swf_plots(1:size(xdata,2)),lg_txt);
        elseif h.radio_3D_true_locs.Value==0 && h.radio_3D_peak_locs.Value==1  && ~isempty(h.current_3D_peak_idx) % only peak waves
            xdata = h.current_norm_peak_swf;
            sel_v = h.listbox_peaks_found.Value;
            lg_txt = num2str(h.inv_soln(h.current_inv_soln).peak_idx(sel_v)');
            legend(h.axes_source_waves,h.current_inv_swf_plots(sel_v),lg_txt);
        elseif h.radio_3D_true_locs.Value==1 && h.radio_3D_peak_locs.Value==1  && ~isempty(h.current_3D_peak_idx) % both
             xdata = h.norm_true_swf;
            lg_txt1=[]; for v=1:size(xdata,2); lg_txt1{v} = sprintf('True Source %.f',v); end
            sel_v = h.listbox_peaks_found.Value;
%             lg_txt2 = cellstr(num2str(h.inv_soln(h.current_inv_soln).peak_idx(sel_v)'))'; 
            lg_txt2 = {};
            for v=1:length(sel_v)
                lg_txt2(v) = { num2str(h.inv_soln(h.current_inv_soln).peak_idx(sel_v(v)))};
            end
            legend(h.axes_source_waves,[h.current_true_swf_plots(1:size(xdata,2)) h.current_inv_swf_plots(sel_v)],[lg_txt1(1:size(xdata,2)) lg_txt2{1:length(sel_v)}]);
        else
            lg_txt=[]; for v=1:size(xdata,2); lg_txt{v} = sprintf('True Source %.f',v); end
            legend(h.axes_source_waves,h.current_true_swf_plots(1:size(xdata,2)),lg_txt);
        end
        h.axes_source_waves.Legend.Position = [0.4875    0.2826    0.1350    0.1233];
        h.axes_source_waves.Legend.FontSize = 6.5;
        
        
    catch me
        bs_plot_peak_waves;
    end
end


sm_set_source_wave_scales('Both');
