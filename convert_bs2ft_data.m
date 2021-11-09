function  ft_data = convert_bs2ft_data(data,anatomy,cfg)
% BrainSim to FieldTrip data format

%% creating FT data format
lat=cfg.study.lat_sim;
% chantype={repmat({'eeg'},[length(anatomy.leadfield.label) 1])};
% chanunit={repmat({'uV'},[length(anatomy.leadfield.label) 1])};
chantype = anatomy.sens.chantype;
chanunit= anatomy.sens.chanunit;
chan_labels={anatomy.leadfield.label'};
s_time={repmat({lat},[1 size(data,3)])};
clear xdata; for t=1:size(data,3); xdata(t)={squeeze(data(:,:,t)')}; end
xdata={xdata};
samp_info={[ones(size(data,3),1) ones(size(data,3),1)+1]};
trial_info={ones(size(data,3),1)};
hdr=struct('Fs',cfg.study.srate,'nChans',length(anatomy.leadfield.label'),'label', chan_labels,'nTrials', size(data,3),'nSamples',length(lat),'nSamplesPre',round(abs(lat(1))*cfg.study.srate),'chantype',chantype,'chanunit',chanunit);
ft_data=struct('hdr',hdr,'label',chan_labels,'time',s_time,'trial',xdata,'fsample',cfg.study.srate,...
    'sampleinfo',samp_info,'trialinfo',trial_info);
