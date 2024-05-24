function sm_plot_fft(varargin)
global h
amp_gain = 1; 

%% setting up Figure
h.fig_fft = figure(300+h.menu_synthetic_noise_cov_type.Value); clf;
ax = subplot_axes(2,2,.1,.1, 0, 0, 0);


%% Getting Data
if h.menu_sens_type.Value==1 % MEG
    amp_units = 'fT'; %amp_gain = 1e15; % converting to fT
elseif h.menu_sens_type.Value==2 % EEG
    amp_units = [char(181) 'V']; %amp_gain = 1e6; % converting to microV
end

noise_fft = calc_fft(h.sim_data.sens_noise_final(:,h.listbox_chans.Value,:),h.cfg.study.srate,0,'k');
signal_fft = calc_fft(h.sim_data.sens_sig_data(:,h.listbox_chans.Value,:),h.cfg.study.srate,0,'k');
final_fft = calc_fft(h.sim_data.sens_final(:,h.listbox_chans.Value,:),h.cfg.study.srate,0,'k');

noise_fft_evk = calc_fft(nanmean(h.sim_data.sens_noise_final(:,h.listbox_chans.Value,:),3),h.cfg.study.srate,0,'k');
signal_fft_evk = calc_fft(nanmean(h.sim_data.sens_sig_data(:,h.listbox_chans.Value,:),3),h.cfg.study.srate,0,'k');
[final_fft_evk,freqs,nfft] = calc_fft(nanmean(h.sim_data.sens_final(:,h.listbox_chans.Value,:),3),h.cfg.study.srate,0,'k');

fft_trial_amp = nanmean(2*abs(cat(4,noise_fft,signal_fft,final_fft)) ,3)*amp_gain; % for less code to plot
fft_trial_phase = nanmean(unwrap(angle(cat(4,noise_fft,signal_fft,final_fft)))  ,3); % for less code to plot
fft_evk_amp = 2*abs(cat(3,noise_fft_evk ,signal_fft_evk ,final_fft_evk))   *amp_gain; % for less code to plot
fft_evk_phase = unwrap(angle(cat(3,noise_fft_evk ,signal_fft_evk ,final_fft_evk)))  ; % for less code to plot

clear noise_fft signal_fft final_fft noise_fft_evk signal_fft_evk final_fft_evk;


%% FFT scales
fft_trial_scale = max(max(max(fft_trial_amp)));
fft_evk_scale = max(max(max(fft_evk_amp)));
amp_scale = max([fft_trial_scale fft_evk_scale])*1.1;

%% Plotting 
% setting up axes
for a=1:4
    axes(ax(a));  cla; hold on; 
    axis(ax(a),'on');
    ax(a).XLim = str2num(h.edit_plot_freq_int.String);
    box(ax(a),'on');
end


ax(1).YLim = [0 amp_scale]; ylabel(ax(1),sprintf('Amplitude (%s)',amp_units));  title(ax(1),'Trial FFT Amplitude'); 
ax(2).YLim = [0 amp_scale]; title(ax(2),'Evoked FFT Amplitude'); 
ax(3).YLim = [-180 180]; ylabel(ax(3),'Phase (degrees)');  xlabel(ax(3),'Frequency (Hz)'); title(ax(3),'Trial FFT Phases'); 
ax(4).YLim = [-180 180]; xlabel(ax(4),'Frequency (Hz)'); title(ax(4),'Trial FFT Phases'); 

clear p1;
plot_order = [3 2 1]; 
ln_width = [1 2 3];  % noise, signal, final 
ln_clr = [.5 .5 .5; 0 .6 0; .5 0 .5];  % noise, signal, final
h.fft_fig_p1 = []; h.fft_fig_p2 = []; h.fft_fig_p3 = [];  h.fft_fig_p4 = []; 
for c=plot_order
    % Trial FFT Amplitudes
    h.fft_fig_p1(c,:) = plot(ax(1),freqs,(fft_trial_amp(1:nfft/2+1,:,c)),'color',ln_clr(c,:),'linewidth',ln_width(c));
    % Evoked FFT Amplitudes
    h.fft_fig_p2(c,:) = plot(ax(2),freqs,fft_evk_amp(1:nfft/2+1,:,c),'color',ln_clr(c,:),'linewidth',ln_width(c));
    %% Trial FFT Phases
    h.fft_fig_p3(c,:) = plot(ax(3),freqs,fft_trial_phase(1:nfft/2+1,:,c),'color',ln_clr(c,:),'linewidth',ln_width(c));
    %% Evoked FFT Phases
    h.fft_fig_p4(c,:) = plot(ax(4),freqs,fft_evk_phase(1:nfft/2+1,:,c),'color',ln_clr(c,:),'linewidth',ln_width(c));
end
h.fft_fig_p1 = handle(h.fft_fig_p1); h.fft_fig_p2 = handle(h.fft_fig_p2); 
h.fft_fig_p3 = handle(h.fft_fig_p3); h.fft_fig_p4 = handle(h.fft_fig_p4); 

legend(ax(1),h.fft_fig_p1(:,1),{'Noise' 'Signal' 'Final'});
%% adding text box to list selected channels
% delete(t1_txt); delete(t1)
t1_txt = uicontrol(h.fig_fft,'Style','text', 'BackgroundColor', h.fig_fft.Color, 'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(ax(2).Position([1 3]))+.01 .9 .1 .035],...
        'FontSize',10,'HorizontalAlignment','left','String','Sensors');

h.fft_fig_chans = uicontrol(h.fig_fft,'Style','listbox', 'BackgroundColor', h.fig_fft.Color, 'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(ax(2).Position([1 3]))+.01 .05 .1 .85],...
        'FontSize',10,'HorizontalAlignment','left','String',sprintf('%s\n',h.listbox_chans.String{h.listbox_chans.Value}),...
        'enable','on','CallBack',{@highlight_chan},'max',1e4);

end

function highlight_chan(varargin)
global h
ln_clr = [.5 .5 .5; 0 .6 0; .5 0 .5];  % noise, signal, final

x_sel = h.fft_fig_chans.Value; 
%     not_sel = setdiff(1:size(p1,2),h.fft_fig_chans.Value); 
    
    for cx=1:size(ln_clr,1)
        for v=1:size(h.fft_fig_p1,2)  % resetting to original colors
        h.fft_fig_p1(cx,v).Color = ln_clr(cx,:);
        h.fft_fig_p2(cx,v).Color = ln_clr(cx,:);
        h.fft_fig_p3(cx,v).Color = ln_clr(cx,:);
        h.fft_fig_p4(cx,v).Color = ln_clr(cx,:);
        end
    end
    for cx=1:size(ln_clr,1)
        for v=1:length(x_sel)
        h.fft_fig_p1(cx,x_sel(v)).Color = [1 0 0]; uistack(h.fft_fig_p1(cx,x_sel(v)),'top');
        h.fft_fig_p2(cx,x_sel(v)).Color = [1 0 0]; uistack(h.fft_fig_p2(cx,x_sel(v)),'top');
        h.fft_fig_p3(cx,x_sel(v)).Color = [1 0 0]; uistack(h.fft_fig_p3(cx,x_sel(v)),'top');
        h.fft_fig_p4(cx,x_sel(v)).Color = [1 0 0]; uistack(h.fft_fig_p4(cx,x_sel(v)),'top');
        end
    end
end
   