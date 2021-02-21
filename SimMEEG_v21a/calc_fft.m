function [fft_data,freqs,nfft]=calc_fft(data,srate,plot_flag,ln_color)
%function [fft_data,freqs]=calc_fft(data,srate,plot_flag,ln_color);

L=length(data);
nfft = 2^nextpow2(L);
fft_data=fft(data,nfft)/L;
freqs=srate/2*linspace(0,1,nfft/2+1);

if plot_flag==1
    if nargin<4
        subplot(2,1,1); plot(freqs,2*abs(fft_data(1:nfft/2+1)),'linewidth',2); ylabel('Amplitude |y(f)|');
        subplot(2,1,2); plot(freqs,angle(fft_data(1:nfft/2+1)),'linewidth',2); xlabel('Frequency (Hz)'); ylabel('Phase');
    else
        subplot(2,1,1); plot(freqs,2*abs(fft_data(1:nfft/2+1,:)),'color',ln_color,'linewidth',2);  ylabel('Amplitude |y(f)|');
        subplot(2,1,2); plot(freqs,unwrap(angle(fft_data(1:nfft/2+1,:))),'color',ln_color,'linewidth',2); xlabel('Frequency (Hz)'); ylabel('Phase');
    end
end

