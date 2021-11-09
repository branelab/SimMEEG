function [SNR,RN,SNR_pval,SNR_dB]=bl_calc_SNR(data,lat,act_samps,ctrl_samps,BandWidth)
% function [SNR,RN,SNR_pval,SNR_dB]=bl_calc_SNR(data,lat,act_samps,ctrl_samps,BandWidth)
%
% This program calculates the Signal-to-Noie Ratio act_samps/ctrl_samps .
%       SNR = std(mean(data(:,act_samps,:),3)/std(mean(data(:,ctrl_samps,:),3)).
%
%INPUT
% data = [chans x samps x trials];
% lat = latency values in data;
% act_samps = samples for active (signal) interval
% ctrl_samps = samples for control (baseline?) interval
% BandWidth = estimated bandwidth of the signal to determine significance values degreeFredom = 2*BandWidth*Time.
%
%
%OUTPUT
% SNR = signal-to-Noise ratio 
% RN = residual noise = std(mean(data(:,ctrl_samps,:),3))
%
%   - SNR can be statistically evaluated on a F-distribution for
%   estimated degrees of freedom = 2*BandWidth*time
%
% created by A. Herdman, May 23, 2014

df = floor(2*BandWidth*(lat(act_samps(end))-lat(act_samps(1)))); % estimated degrees of freedom
SNR = std(mean(data(:,act_samps,:),3),[],2)./std(mean(data(:,ctrl_samps,:),3),[],2);
SNR_dB = 20*log10(std(mean(data(:,act_samps,:),3),[],2)./std(mean(data(:,ctrl_samps,:),3),[],2));
RN=std(mean(data(:,ctrl_samps,:),3),[],2);
SNR_pval = 1 - fcdf(SNR.^2,1,df); % probability of reponse being present



