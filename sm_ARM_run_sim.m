function sm_ARM_run_sim(varargin)
global h

% first gather params and set sources



%% Updating ARM params
h.cfg.ARM_params.num_sources = str2num(h.edit_ARM_num_sources.String);
h.cfg.ARM_params.num_samps = length(h.cfg.study.lat_sim);
h.cfg.ARM_params.num_trials = h.cfg.study.num_trials;
h.cfg.ARM_params.ARMorder = str2num(h.edit_ARM_order.String);
h.cfg.ARM_params.num_interactions = str2num(h.edit_ARM_num_interactions.String);
h.cfg.ARM_params.ori_constraint = h.menu_ARM_ori_constraint.String{h.menu_ARM_ori_constraint.Value};
h.cfg.ARM_params.source_locations_type = h.menu_ARM_source_locs.String{h.menu_ARM_source_locs.Value};

% [data, Arsig, ~, h.cfg.ARM_params.lambdamax, h.cfg.ARM_params.sigma] = arm_generate_signal(h.cfg.ARM_params.num_sources, h.cfg.ARM_params.num_samps*h.cfg.ARM_params.num_trials, h.cfg.ARM_params.ARMorder, h.cfg.ARM_params.num_interactions);
pad_samps = 2500; % iniital padded samples needed for ARM calculations

%% currently using narrow-band noise by filtering --> need to add other noise types broadband, notched band, pink 
data_wav = randn(h.cfg.ARM_params.num_sources, h.cfg.ARM_params.num_samps*h.cfg.ARM_params.num_trials+pad_samps);

%% filtering initial data wav
filt_type = 'bandpassfir'; filt_method = 'equiripple';
filt_dB = 30; filt_order = [];
freqs = str2num(h.edit_ARM_freq_range.String);
sfreq = freqs(1)-1; if sfreq<=0;sfreq=0.1; end % need 
efreq = freqs(2)+2; if efreq>(h.cfg.study.srate/2); efreq=(h.cfg.study.srate/2); end 
freq_band = [sfreq freqs efreq];
fprintf('Creating Filter for initial ARM waveform\n'); 
[h.ARM_filt_design, data_wav] = filter_data_new([],double(data_wav'),h.cfg.study.srate,freq_band,filt_type,filt_method,1,filt_dB,filt_order,0);
data_wav = data_wav'; 

%% calculating ARM data
[data, Arsig, ~, h.cfg.ARM_params.lambdamax, h.cfg.ARM_params.sigma] = arm_generate_signal(h.cfg.ARM_params.num_sources, h.cfg.ARM_params.num_samps*h.cfg.ARM_params.num_trials, h.cfg.ARM_params.ARMorder, h.cfg.ARM_params.num_interactions, pad_samps, data_wav);
h.cfg.ARM_params.interactions_matrix = reshape(Arsig, h.cfg.ARM_params.num_sources, h.cfg.ARM_params.num_sources, h.cfg.ARM_params.ARMorder);

%% Plotting ARM interactions
h.axes_ARM_interactions.NextPlot = 'replace';
imagesc(h.axes_ARM_interactions,mean(abs(h.cfg.ARM_params.interactions_matrix),3) ~= 0); 
h.axes_ARM_interactions.NextPlot = 'add';
Xgrid = .5:1:h.cfg.ARM_params.num_sources+.5; Yvals = repmat(h.axes_ARM_interactions.YLim,length(Xgrid),1);
plot(h.axes_ARM_interactions,[Xgrid; Xgrid],Yvals','k');
Ygrid = .5:1:h.cfg.ARM_params.num_sources+.5; Xvals = repmat(h.axes_ARM_interactions.YLim,length(Xgrid),1);
plot(h.axes_ARM_interactions,Xvals',[Ygrid; Ygrid],'k');
colormap(h.axes_ARM_interactions,[1 1 1; 0 .6 0]);
h.axes_ARM_interactions.XTick = 0:1:h.cfg.ARM_params.num_sources;
h.axes_ARM_interactions.Title.String = 'ARM Interactions';
xlabel(h.axes_ARM_interactions,'Source #'); ylabel(h.axes_ARM_interactions,'Source #'); 



%% normalizing to +/-1 scale then multiplying by ARM signal percentage
max_val = max(abs([max(data,[],'all') min(data,[],'all')]));
% max_val = max(quantile(data',.999)); % quantile helps get more consistent amplitudes across sources due to a few acute jumps in amplitudes for high freqs
data = data/max_val';  
% data = bsxfun(@mtimes,data',sig_perc); data = data'; 
% for v=1:size(data,1); data(v,:) = data(v,:).*sig_perc(v); end

%% reshaping and putting into correct order for SimMEEG
h.sim_data.ARM_source_sig_data = permute(reshape(data,[h.cfg.ARM_params.num_sources h.cfg.ARM_params.num_samps h.cfg.ARM_params.num_trials]),[2 1 3]);
tb_data = h.table_ARM_source_params.Data;

win_type = h.menu_ARM_win_type.String{h.menu_ARM_win_type.Value};

for v=1:size(h.sim_data.ARM_source_sig_data,2) % Num sources
    sig_time = tb_data(v,10:11);
    rise_time = tb_data(v,12);
    fall_time = rise_time;
    pre_perc = tb_data(v,9);
    post_perc = pre_perc;
    sig_perc = tb_data(v,8);
%     v_amp = tb_data(v,7);
    [hwin] = sm_create_amp_window(win_type, h.cfg.study.lat_sim, h.cfg.study.srate, sig_time, rise_time, fall_time, pre_perc, sig_perc, post_perc);
    data = squeeze(h.sim_data.ARM_source_sig_data(:,v,:));
    data2 = data.*repmat(hwin,1,size(data,2));
    h.sim_data.ARM_source_sig_data(:,v,:) = data2; 
%     h.sim_data.ARM_source_sig_final(:,v,:) = data2*v_amp; 
end

%% adding noise
num_sources = size(h.sim_data.ARM_source_sig_data,2); 
cfg = h.cfg; 
% changing cfg based on noise params from Panel "Study Parameters" not from Panel "Generate Sensor Noise".
sig_freqs = repmat(freqs,[num_sources 1 1]); cfg.source.sig_freqs = permute(sig_freqs,[1 3 2]); 
cfg.study.synthetic_noise_flag = h.edit_noise_flag.Value;
cfg.study.synthetic_noise_freqs = str2num(h.edit_noise_freqs.String);
cfg.study.synthetic_pink_noise_slope = h.cfg.study.pink_noise_slope;

noise_wav = sim_noise_wav(cfg, size(h.sim_data.ARM_source_sig_data,2) );
xperc = str2num(h.edit_noise_amp_perc.String)/100;

noise_wav = noise_wav * xperc;

h.sim_data.ARM_source_sig_data = h.sim_data.ARM_source_sig_data + noise_wav;




