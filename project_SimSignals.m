function sens_data = project_SimSignals(data,leadfield,vx_idx,source_amp,source_ori, avgref_flag)
%
% INPUT:
%   data = [samps x sources x trials] source waveforms to project to sensor space - should have values scaled from -1.0 to 1.0
%   leadfield.H = [sensors x sources x voxels]
%   vx_idx = indices of voxels from leadfield where sources are located
%   source_amp = source amplitudes (nA) to gain the data before projecting to sensors
%   source_ori = orientations of sources for each vx_idx
%
% OUTPUT
%   sens_data = [samps x sensors x trials]
%

global h

%% projecting source data to sensors
eeg2=zeros(size(data,1),size(leadfield.H,1),size(data,3),length(vx_idx)); 
%% conversion factor needed for EEG/MEG because Field Trips lead fields appear to be different than in Brainstorm (openmeeg) --> Not sure why this is?
% adjustment needed by user to get reasonable 150-200 fT for 35-40 nA dipoles in audiotry cortices like in Herdman et al., 2003 NeuroImage dataset
% or for EEG 10-20 microV for 100 nA dipoles in audiotry cortices like in Herdman & Takai 2013
lf_gain = str2num(h.edit_leadfield_gain.String);

%% This is not working for leadfields calculated using FieldTrip
% if max(max(max(data>1e-6)))
%     source_amp = source_amp*1e-9; % converting to Amps --> source data units must be in Amps to return sensor data in Tesla (MEG) or Volts (EEG)
% end

%% forward projecting the source data to the sensors
for v=1:size(data,2)
    data(:,v,:)=data(:,v,:)*source_amp(v);  % converting to nAmps
    lf(:,v) = leadfield.H(:,:,vx_idx(v))*source_ori(v,:)';
end
lf = lf*lf_gain;
if avgref_flag==1   % leadfields for average reference
    lf = lf - mean(lf);
end

for t=1:size(data,3)
    sens_data(:,:,t) = squeeze(data(:,:,t))*lf';
end
h.lf = lf;
h.swf = mean(data,3); 
