function h = sm_batch_sm_load_real_sensor_file(varargin)
% this program loads in real data from a .mat file with the following variables:
%
%   sensor_file.mat
%       data = [sensor x samples x trials]
%       srate = sample rate
%       lat = latency values of samples
h = varargin{1};

% msgbox('Not yet implemented for commnad-line batch scripting. Please use "Synthetic Data" simulations options','sm_batch_sm_load_real_sensor_file.m'); 
% 
% return

study_trials = str2num(h.edit_num_trials.String);
study_srate = str2num(h.edit_srate.String);

if h.monte_carlo_flag == 1 && h.menu_monte_synthetic_real_data.Value==4  % load Real Sensor Data in Monte-Carlo sims
    [fpath,fname] = fileparts(h.monte_real_dsname);
    h.waitfor_txt.String = sprintf('Loading Real Sensor data from file: %s\n',fullfile(fpath,fname)); drawnow;
else                                                                     % load Real Sensor Data in manual mode
    [fname,fpath]=uigetfile('*.mat','Load Real Sensor Dataset');
    h.waitfor_panel.Visible='on';
    h.waitfor_txt.String = sprintf('Loading Real Sensor data from file: %s\n',fullfile(fpath,fname)); drawnow;
end

x = load(fullfile(fpath,fname),'data','srate','lat');
% h.anatomy.sens.bad_sensors = find(isnan(squeeze(x.data(:,1,1)))==1);    % finding bad sensors for data with nans
% update_anatomy_fields;
% if ~isempty(bad_chans) % checking for bad channels being set to nans
%     x.data = sm_interpolate_data(x.data,bad_chans);
% end

%% data structure needs to be [chans x samples x trials];
xidx = find(size(x.data) == length(x.lat));
if xidx==1 % samples are in first index --> nee to shift dims
    x.data = permute(x.data,[2 1 3]);
elseif xidx==3
    x.data = permute(x.data,[1 3 2]);
end

%%
if x.srate > study_srate  % downsampling the data
    x.data = resample(double(x.data), study_srate , x.srate, 'Dimension',2);
    x.lat = resample(double(x.lat), study_srate , x.srate, 'Dimension',2);
    x.srate = study_srate;
    
    % trying to align timings - if not selecting first samples
    % see if all study.lat_sim are within x.lat values
    lat_sim = h.cfg.study.lat_sim;
    if x.lat(1)<lat_sim(1) && x.lat(end)>lat_sim(end)   % study.lat_sim lies within x.lat so data can be extracted
        xs = find(x.lat<lat_sim(1)); xs1 = xs(end);
        xs = find(x.lat<lat_sim(end)); xs2 = xs(end);
        x.data = x.data(:,xs1:xs2,:);
        x.lat = lat_sim;
    elseif length(x.lat) >= length(lat_sim)
        x.data = x.data(:,1:length(lat_sim),:);
        x.lat = lat_sim;
    end
else
end

if size(x.data,2)>length(h.cfg.study.lat_sim)
    x.data = x.data(:,1:length(h.cfg.study.lat_sim),:);
end

%% randomly selecting number of trials
if x.srate==study_srate && size(x.data,2)<=length(h.cfg.study.lat_sim) && size(x.data,3)>=study_trials && size(x.data,1)>=length(h.anatomy.leadfield.label) % >= num of sensors used for leadfields
    
    nt = randperm(size(x.data,3));  % randomly selecting trials
    chan_idx = 1:length(h.anatomy.leadfield.label); % selecting first sensors
    xsamps = 1:length(h.cfg.study.lat_sim);     % selecting first samples
    fprintf('Sensor indices for sensor %s data: ',varargin{3}); fprintf('%.f ',chan_idx); fprintf('\n');
    fprintf('Randomly selecting %.f of %.f trials\n',study_trials,size(x.data,3))
    
    switch varargin{3}
        case 'noise'
            h.sim_data.sens_noise = permute(x.data(chan_idx,xsamps,nt(1:study_trials)),[2 1 3]);
            h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
            h.sim_data.sens_noise_final = h.sim_data.sens_noise_scaled;
            h.sim_data.sens_noise_final_org = h.sim_data.sens_noise_final;
            h.sim_data.sens_noise_type = h.menu_noise_projection.String{h.menu_noise_projection.Value};
        case 'final'
            h.sim_data.sens_final = permute(x.data(chan_idx,xsamps,nt(1:study_trials)),[2 1 3]);
            h.sim_data.sens_sig_data = zeros(size(h.sim_data.sens_final));
            h.sim_data.sens_noise_final = zeros(size(h.sim_data.sens_final));
            
            h.sim_data.sens_final_org = h.sim_data.sens_final;  % setting orignal data to newly load M/EEG dataset
            h.sim_data.sens_noise_final_org = h.sim_data.sens_noise_final;  % removing noise
            h.sim_data.sens_sig_data_org = h.sim_data.sens_sig_data;   % removing sensor data from projected sources
            h.sim_data.sens_noise_type      = 'None - Loaded in Real Full Data';
            h.sim_data.source_waveform_type = 'None - Loaded in Real Full Data';
    end
    
    
else
    s1 = sprintf('Some Data File parameters are NOT compatible with Study Parameters:\n\n                         Data File          Study Parameters\n');
    if x.srate~=study_srate; s2 = sprintf('Sample Rate    %.f       ~=      %s   (No)\n',x.srate,h.edit_srate.String);
    else; s2 = sprintf('Sample Rate     %.f       ==      %s   (Yes)\n',x.srate,h.edit_srate.String); end
    
    if size(x.data,2)~=length(h.cfg.study.lat_sim); s3 = sprintf('Num Samples   %.f       <       %.f   (No)\n',size(x.data,2),length(h.cfg.study.lat_sim));
    else;  s3 = sprintf('Num Samples   %.f       >=      %.f   (Yes)\n',size(x.data,2),length(h.cfg.study.lat_sim));  end
    
    if size(x.data,1)<length(h.anatomy.leadfield.label); s4 = sprintf('Num Sensors   %.f      <       %.f   (No)\n',size(x.data,1),length(h.anatomy.leadfield.label));
    else;  s4 = sprintf('Num Sensors    %.f      >=      %.f   (Yes)\n',size(x.data,1),length(h.anatomy.leadfield.label));   end
    
    if size(x.data,3)<study_trials; s5 = sprintf('Num Trials    %.f       <       %.f   (No)\n',size(x.data,3),study_trials);
    else;  s5 = sprintf('Num Trials        %.f       >=      %.f   (Yes)\n',size(x.data,3),study_trials);   end
    sx = [s1, s2, s3, s4, s5];
    w = warndlg(sprintf('%s\n',sx));
    w.Position(3:4)=[400 200]; htext = findobj(w, 'Type', 'Text'); htext.FontSize = 11; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
end

