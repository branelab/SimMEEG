function sens_data = project_SimSignals(data,leadfield,vx_idx,source_amp,source_ori)
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
%% conversion factor needed for EEG because Field Trips lead fields appear to be different than in Brainstorm (openmeeg) --> Not sure why this is?
% adjustment needed by user to get reasonable 150-200 fT for 35-40 nA dipoles in audiotry cortices like in Herdman et al., 2003 NeuroImage dataset
lf_gain = str2num(h.edit_leadfield_gain.String); 

%% forward projecting the source data to the sensors
for v=1:size(data,2)
    data(:,v,:)=data(:,v,:)*source_amp(v);  % converting to nAmps
    for hv=1:size(leadfield.H,1)
        d1=(data(:,v,:)*leadfield.H(hv,1,vx_idx(v)))*source_ori(v,1)*lf_gain; % contribution of dipole for source orientation aligned to 1st axis.
        d2=(data(:,v,:)*leadfield.H(hv,2,vx_idx(v)))*source_ori(v,2)*lf_gain; % contribution of dipole for source orientation aligned to 2nd axis.
        d3=(data(:,v,:)*leadfield.H(hv,3,vx_idx(v)))*source_ori(v,3)*lf_gain; % contribution of dipole for source orientation aligned to 3rd axis.
        eeg2(:,hv,:,v)=(d1+d2+d3);
    end
end

%% adjusting Noise to yield SNR
sens_data=nansum(eeg2,4); 
sens_data(sens_data==0)=nan;

