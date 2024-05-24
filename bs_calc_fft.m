function bs_calc_fft(varargin)
global h

h.axes_source_fft.clo;
h.axes_source_fft.NextPlot='replace';
if h.radio_3D_true_locs.Value==1
    for v=1:size(h.norm_true_swf,2)
        %% True Sources
        if h.radio_normalize_swf.Value == 1 % FFT on normalized waveforms
            data = h.norm_true_swf(:,v);
        else
            data = bsxfun(@minus,h.sim_data.sig_final(:,v,:),nanmean(h.sim_data.sig_final(h.sim_data.cfg.study.base_samps,v,:)))*h.sim_data.cfg.source.vx_amp(v);
            data = squeeze(nanmean(data,3));
        end
        L=length(data);
        nfft = 2^nextpow2(L);
        fft_data=fft(data,nfft)/L;
        fft_freqs=h.sim_data.cfg.study.srate/2*linspace(0,1,nfft/2+1);
        
        %% plotting data
        h.current_inv_true_swf_plots_fft(v) = plot(h.axes_source_fft,fft_freqs,2*abs(fft_data(1:nfft/2+1,:)),'color',h.ln_clr(v,:),'linewidth',1);
        h.axes_source_fft.XLabel.String= 'Frequency (Hz)'; h.axes_source_fft.YLabel.String= 'Amplitude |y(f)|';
        h.axes_source_fft.NextPlot='add';
        
        h.sim_data.cfg.source.true_evk_fft_data(:,v) = fft_data;
        h.sim_data.cfg.source.true_evk_fft_freqs = fft_freqs;

    end
end

if h.radio_3D_peak_locs.Value==1
    h.inv_soln(h.current_inv_soln).swf_evk_fft_data = [];
    %% Peak Sources
    for v=1:size(h.current_peak_swf,2)
        if h.radio_normalize_swf.Value == 1 % FFT on normalized waveforms
            data = h.current_norm_peak_swf(:,v);
        else
            data = h.current_peak_swf(:,v);
        end
        
        L=length(data);
        nfft = 2^nextpow2(L);
        fft_data=fft(data,nfft)/L;
        fft_freqs=h.sim_data.cfg.study.srate/2*linspace(0,1,nfft/2+1);
        
        %% plotting data
        h.current_inv_swf_plots_fft(v) = plot(h.axes_source_fft,fft_freqs,2*abs(fft_data(1:nfft/2+1,:)),'color',h.ln_clr(v,:),'linewidth',2);
        h.current_inv_swf_plots_fft = handle(h.current_inv_swf_plots_fft); 
        h.current_inv_swf_plots_fft(v).Color(4) = h.false_positive_lineAlpha;
        h.axes_source_fft.XLabel.String= 'Frequency (Hz)'; h.axes_source_fft.YLabel.String= 'Amplitude |y(f)|';
        h.axes_source_fft.NextPlot='add';
        
        h.inv_soln(h.current_inv_soln).swf_evk_fft_data(:,v) = fft_data; 
        h.inv_soln(h.current_inv_soln).swf_evk_fft_freqs = fft_freqs; 
        
    end
end
box on; axis on;
h.axes_source_fft.XLim = str2num(h.edit_plot_freq_int.String);
if h.radio_normalize_swf.Value == 1 % FFT on normalized waveforms
    h.axes_source_fft.Title.String='FFT: Normalized Evoked'; h.axes_source_fft.YLabel.String= 'Normalized Amplitude';
else
    h.axes_source_fft.Title.String='FFT: Evoked'; h.axes_source_fft.YLabel.String= 'Amplitude';
end

disableDefaultInteractivity(h.axes_source_fft); h.axes_source_fft.Toolbar.Visible='off'; 

