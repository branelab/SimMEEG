function h = sm_batch_combine_sens_data(h)
snr_val = str2num(h.edit_sens_SNR.String);

dims1 = size(h.sim_data.sens_sig_data);
dims2 = size(h.sim_data.sens_noise_scaled);
if sum(dims1~=dims2)>0
    errordlg(sprintf('Dimensions of sensor signals (%.f %.f %.f) and noise (%.f %.f %.f) are not compatible',dims1,dims2));
    return;
end
% if snr_val<1 && snr_val>-1
%     errordlg(sprintf('Signal-to-Noise (dB) cannot be between -1 and 1\n\nSNR = %.f',snr_val),'WARNING!');
% end
%% simulating sensor noise WITHOUT spatiotemporal covariance
if ~isfield(h.sim_data,'sens_noise_scaled')
    sim_sens_noise();
end

%% SNR target: using sum-of-squares in active vs. control to gain sensor noise to match target SNR val (in dB)
snr_thresh=.05;
act_int = str2num(h.edit_poststim_int.String);
act_s = round( (act_int-h.cfg.study.lat_sim(1))*h.cfg.study.srate );  act_s(act_s==0)=1;
h.cfg.study.act_samps = act_s(1):act_s(2);
base_s = round( (h.cfg.study.base_int-h.cfg.study.lat_sim(1))*h.cfg.study.srate ); base_s(base_s==0)=1;
h.cfg.study.base_samps = base_s(1):base_s(2);

%% baselining signal and noise to get proper estimates for Sum-of-Squares
h.sim_data.sens_noise_scaled = bsxfun(@minus, h.sim_data.sens_noise_scaled, nanmean(h.sim_data.sens_noise_scaled(h.cfg.study.base_samps,:,:))); 
h.sim_data.sens_sig_data = bsxfun(@minus, h.sim_data.sens_sig_data, nanmean(h.sim_data.sens_sig_data(h.cfg.study.base_samps,:,:))); 

% Use mean(SNR) or max(SNR) across sensors ??
brn = h.sim_data.sens_noise_scaled;
h.sim_data.sens_final = h.sim_data.sens_sig_data + brn;

% s_ss = sum(sum(sum(h.sim_data.sens_final(h.cfg.study.act_samps,:,:).^2)));
s_ss = sum(sum(sum(h.sim_data.sens_sig_data(h.cfg.study.act_samps,:,:).^2)));
n_ss = sum(sum(sum(h.sim_data.sens_final(h.cfg.study.base_samps,:,:).^2)));
snr_ss = s_ss/n_ss; 
snr_db = 10*log10(snr_ss);
db_diff = snr_db - snr_val; 

%% while loop
t=0;    % Stops while loop after 1000 iterations
while abs(db_diff)>snr_thresh
    t=t+1;
    
    noise_gain = sqrt(10^(db_diff/10));
    brn = brn*noise_gain;
    h.sim_data.sens_final = h.sim_data.sens_sig_data + brn;
    
    s_ss = sum(sum(sum(h.sim_data.sens_sig_data(h.cfg.study.act_samps,:,:).^2)));
    n_ss = sum(sum(sum(h.sim_data.sens_final(h.cfg.study.base_samps,:,:).^2)));
    snr_ss = s_ss/n_ss;
    snr_db = 10*log10(snr_ss);
    db_diff = snr_db - snr_val;

   
    if t>200
        errordlg(sprintf('SNR (dB) = %.2f was not found to be within %.2f +/- %.2f \n\n',snr_db,snr_val,snr_thresh));
        break
    end
end

h.sim_data.sens_noise_final = brn;








