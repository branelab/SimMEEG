function [SNR, SNR_dB, SNR_hilbert_trial_act_ctrl, SNR_hilbert_evk_act_ctrl, SNR_hilbert_evk_act, SNR_hilbert_evk_ctrl] = sm_calc_SNR(data,act_samps,ctrl_samps)
% function [SNR, SNR_dB, SNR_hilbert_trial_act_ctrl, SNR_hilbert_evk_act_ctrl, SNR_hilbert_evk_act, SNR_hilbert_evk_ctrl]  = sm_calc_SNR(data,act_samps,ctrl_samps);
% This is a side program that calculates the SNR as done in SimMEEG
%   calculates sum of squares across act_samps and ctrl_samps then takes the ratio SNR = act/ctrl
%
% INPUT
%   data = [samps x chans x trials]
%   act_samps = samples in active interval for numerator of SNR
%   act_samps = samples in contrla interval for denominator of SNR
%
% OUTPUT
%   SNR = sum(sum(sum(data(act_samps,:,:).^2))) / sum(sum(sum(data(ctrl_samps,:,:).^2))); 
%   SNR_dB = 10*log10(SNR); 
%   SNR_hilbert_evk_act = sum(nanmean(envlp_data_evk(act_samps,:,:),3))./sum(nanmean(envlp_data(act_samps,:,:),3)); % envlp_data_evk = hilbert transform of evoked (Averaged-trial) data
%   SNR_hilbert_evk_ctrl = sum(nanmean(envlp_data_evk(ctrl_samps,:,:),3))./sum(nanmean(envlp_data(ctrl_samps,:,:),3));
%   SNR_hilbert_trial_act_ctrl = sum(nanmean(envlp_data(act_samps,:,:),3))./sum(nanmean(envlp_data(ctrl_samps,:,:),3));  % envlp_data = hilbert transform of data
%   SNR_hilbert_evk_act_ctrl = sum(envlp_data_evk(act_samps,:,:))./sum(nanmean(envlp_data_evk(ctrl_samps,:,:),3));
%
%   Note, use 20*log10(SNR_hilbert_evk_act) to get dB because hilbert transform results are amplitude not amplitude^2
%% act_samps vs ctrl samps
SNR = sum(sum(sum(data(act_samps,:,:).^2))) / sum(sum(sum(data(ctrl_samps,:,:).^2))); 
SNR_dB = 10*log10(SNR);

% Evoked SNR in active interval
[~,envlp_data]=calc_hilbert_amp_phase(data);
[~,envlp_data_evk]=calc_hilbert_amp_phase(nanmean(data,3));
SNR_hilbert_trial_act_ctrl  = sum(nanmean(envlp_data(act_samps,:,:),3))./sum(nanmean(envlp_data(ctrl_samps,:,:),3));
SNR_hilbert_evk_act_ctrl    = sum(envlp_data_evk(act_samps,:,:))./sum(nanmean(envlp_data_evk(ctrl_samps,:,:),3));
SNR_hilbert_evk_act         = sum(nanmean(envlp_data_evk(act_samps,:,:),3))./sum(nanmean(envlp_data(act_samps,:,:),3));
SNR_hilbert_evk_ctrl        = sum(nanmean(envlp_data_evk(ctrl_samps,:,:),3))./sum(nanmean(envlp_data(ctrl_samps,:,:),3));



