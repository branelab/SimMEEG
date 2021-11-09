function bs_calc_fft(varargin)
global h

ln_clr =  h.map3D_peak_locs(1).CData; %lines(length(h.current_3D_peak_idx)); ln_clr(1:3,:) = h.src_clr;

h.axes_source_fft.clo;
h.axes_source_fft.NextPlot='replace';
if h.radio_3D_true_locs.Value==1
    for v=1:size(h.norm_true_swf,2)
        %% True Sources
        data = h.norm_true_swf(:,v);
        L=length(data);
        nfft = 2^nextpow2(L);
        fft_data=fft(data,nfft)/L;
        fft_freqs=h.sim_data.cfg.study.srate/2*linspace(0,1,nfft/2+1);
        
        %% plotting data
        h.current_inv_true_swf_plots_fft(v) = plot(h.axes_source_fft,fft_freqs,2*abs(fft_data(1:nfft/2+1,:)),'color',ln_clr(v,:),'linewidth',1);
        h.axes_source_fft.XLabel.String= 'Frequency (Hz)'; h.axes_source_fft.YLabel.String= 'Amplitude |y(f)|';
        h.axes_source_fft.NextPlot='add';
    end
end

if h.radio_3D_peak_locs.Value==1
    %% Peak Sources
    for v=1:size(h.current_peak_swf,2)
        data = h.current_peak_swf(:,v);
        L=length(data);
        nfft = 2^nextpow2(L);
        fft_data=fft(data,nfft)/L;
        fft_freqs=h.sim_data.cfg.study.srate/2*linspace(0,1,nfft/2+1);
        
        %% plotting data
        h.current_inv_swf_plots_fft(v) = plot(h.axes_source_fft,fft_freqs,2*abs(fft_data(1:nfft/2+1,:)),'color',ln_clr(v,:),'linewidth',2);
        h.current_inv_swf_plots_fft(v).Color(4) = h.false_positive_lineAlpha;
        h.axes_source_fft.XLabel.String= 'Frequency (Hz)'; h.axes_source_fft.YLabel.String= 'Amplitude |y(f)|';
        h.axes_source_fft.NextPlot='add';
        
        
    end
end
box on; axis on;
h.axes_source_fft.XLim = str2num(h.edit_plot_freq_int.String);
h.axes_source_fft.Title.String='FFT';

disableDefaultInteractivity(h.axes_source_fft); h.axes_source_fft.Toolbar.Visible='off'; 

