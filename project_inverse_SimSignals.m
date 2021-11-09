function source_data = project_inverse_SimSignals(sens_data,leadfield,vx_idx,source_ori,lf_gain)
% Inverse porjection of sensor data into sources --> reverse of project_SimSignals.m
%
% INPUT:
%   sens_data = [samps x sensors x trials] sensor waveforms to project to source space
%   leadfield.H = [sensors x sources x voxels]
%   vx_idx = indices of voxels from leadfield where sources are located
%   source_ori = orientations of sources for each vx_idx
%   lf_gain = leadfield gain adjustment that was used to get reasonable nA to uV (or fT) sensor scale 
%
% OUTPUT
%   source_data = [samps x sources x trials]
%

global h

%% projecting sens_data to sources at vx_idx with orientations of source_ori
source_data=zeros(size(sens_data,1),length(vx_idx),size(sens_data,3));

%% conversion factor needed for EEG because Field Trips lead fields appear to be different than in Brainstorm (openmeeg) --> Not sure why this is?
% adjustment needed by user to get reasonable 150-200 fT for 35-40 nA dipoles in audiotry cortices like in Herdman et al., 2003 NeuroImage dataset
% lf_gain = str2num(h.edit_leadfield_gain.String);

%% inverse projection of sensor data to source space using forward leadfields
inv_wts = nan(size(sens_data,2),length(vx_idx));
for v=1:length(vx_idx)
    inv_wts(:,v) = leadfield.H(:,:,vx_idx(v)) * source_ori(v,:)' *lf_gain;
end
% double checking that inv_wts are correct by reversing the forward projected sens_sig_data to get inverse weights (inv_wts)
% sens_sig_data = h.sim_data.sens_sig_data; sig_data = h.sim_data.sig_final;
% inv_wts = nanmean(sens_sig_data,3)' / nanmean(sig_data,3)';

for t=1:size(source_data,3)
    source_data(:,:,t) = (squeeze(sens_data(:,:,t))) / (inv_wts)'; % vector inversion process
end



