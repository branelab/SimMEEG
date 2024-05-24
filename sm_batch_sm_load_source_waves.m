function h = sm_batch_sm_load_source_waves(h)
% this program loads in source waveform data from a .mat file with the following variables:
%
%   source_wave.mat
%       source_data = [sources x samples x trials]
%       srate = sample rate
%       lat = latency values of samples


study_trials = str2num(h.edit_num_trials.String);
study_srate = str2num(h.edit_srate.String);

[fpath,fname,fext] = fileparts(h.monte_params.monte_real_dsname);
x = load(fullfile(fpath,fname),'source_data','srate','lat');
if ~isfield(x,'source_data')
    fprintf('Variable "source_data" not found in file: %s\n',fullfile(fpath,[fname fext]));
else
    
    if x.srate > study_srate  % downsampling the data
        
        x.source_data = resample(double(x.source_data), study_srate , x.srate, 'Dimension',2);
        x.lat = resample(double(x.lat), study_srate , x.srate, 'Dimension',2);
        x.srate = study_srate;
        
        % trying to align timings - if not selecting first samples
        % see if all study.lat_sim are within x.lat values
        lat_sim = h.cfg.study.lat_sim;
        if x.lat(1)<lat_sim(1) && x.lat(end)>lat_sim(end)   % study.lat_sim lies within x.lat so data can be extracted
            xs = find(x.lat<lat_sim(1)); xs1 = xs(end);
            xs = find(x.lat<lat_sim(end)); xs2 = xs(end);
            x.source_data = x.source_data(:,xs1:xs2,:);
            x.lat = lat_sim;
        elseif length(x.lat) >= length(lat_sim)
            x.source_data = x.source_data(:,1:length(lat_sim),:);
            x.lat = lat_sim;
        end
    else
    end
    
    %% randomly selecting number of trials
    if x.srate==study_srate && size(x.source_data,2)<=length(h.cfg.study.lat_sim) && size(x.source_data,3)>=study_trials && size(x.source_data,1)>=size(h.cfg.source.vx_locs,1)
        nt = randperm(size(x.source_data,3));  % randomly selecting trials
        s_idx = 1:size(h.cfg.source.vx_locs,1) ; % selecting first sensors
        xsamps = 1:length(h.cfg.study.lat_sim);     % selecting first samples
        fprintf('Source indices for source waveforms: '); fprintf('%.f ',s_idx); fprintf('\n');
        fprintf('Randomly selecting %.f of %.f trials\n',study_trials,size(x.source_data,3))
        h.sim_data.sig_final = permute(x.source_data(s_idx,xsamps,nt(1:study_trials)),[2 1 3]);
        %     h.sim_data.sig_final = (h.sim_data.sig_final./max(max(max(abs(h.sim_data.sig_final))))); % scaled between -1 to 1
        h.sim_data.sig_wav = zeros(size(h.sim_data.sig_final));
        h.sim_data.sig_win = h.sim_data.sig_wav;
        h.sim_data.prepost_wav = h.sim_data.sig_wav;
        h.sim_data.prepost_win = h.sim_data.sig_wav;
        h.sim_data.source_waveform_type = 'Real Sources';
        
        %% overwriting sig_final_org waves
        h.sim_data.sig_final_org = h.sim_data.sig_final;
        %         h.sim_data,'sig_final_org'
        
        
    else
        s1 = sprintf('Some Data File parameters are NOT compatible with Study Parameters:\n\n                         Data File          Study Parameters\n');
        if x.srate~=study_srate; s2 = sprintf('Sample Rate    %.f       ~=      %s   (No)\n',x.srate,h.edit_srate.String);
        else; s2 = sprintf('Sample Rate     %.f       ==      %s   (Yes)\n',x.srate,h.edit_srate.String); end
        if size(x.source_data,2)~=length(h.cfg.study.lat_sim); s3 = sprintf('Num Samples   %.f       <       %.f   (No)\n',size(x.source_data,2),length(h.cfg.study.lat_sim));
        else;  s3 = sprintf('Num Samples   %.f       >=      %.f   (Yes)\n',size(x.source_data,2),length(h.cfg.study.lat_sim));  end
        if size(x.source_data,1)~=size(h.cfg.source.vx_locs,1); s4 = sprintf('Num Sources   %.f      <      %.f   (No)\n',size(x.source_data,1),size(h.cfg.source.vx_locs,1));
        else;  s4 = sprintf('Num Sensors    %.f      <=      %.f   (Yes)\n',size(x.source_data,1),size(h.cfg.source.vx_locs,1));   end
        if size(x.source_data,3)<study_trials; s5 = sprintf('Num Trials    %.f       <       %.f   (No)\n',size(x.source_data,3),study_trials);
        else;  s5 = sprintf('Num Trials        %.f       >=      %.f   (Yes)\n',size(x.source_data,3),study_trials);   end
        sx = [s1, s2, s3, s4, s5];
        w = warndlg(sprintf('%s\n',sx));
        w.Position(3:4)=[400 200]; htext = findobj(w, 'Type', 'Text'); htext.FontSize = 11; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
    end
end
