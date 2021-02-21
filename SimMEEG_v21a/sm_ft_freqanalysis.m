function sm_ft_freqanalysis(varargin)
global h


foilim = str2num(h.edit_preproc_ft_tfr_foilim.String);


%% ft_freqanalysis %
cfg=[];
cfg.method      = h.menu_preproc_ft_tfr_method.String{h.menu_preproc_ft_tfr_method.Value};
cfg.output      = h.menu_preproc_ft_tfr_output.String{h.menu_preproc_ft_tfr_output.Value};  % gives the complex Fourier spectra
% if length(foilim)==2    % calculating based on frequency band
%   cfg.foilim      = foilim;
% elseif length(foilim)~=2
cfg.foi = foilim;
cfg.taper       = h.menu_preproc_ft_tfr_taper.String{h.menu_preproc_ft_tfr_taper.Value};
cfg.tapsmofrq   = str2num(h.edit_preproc_ft_tfr_tapsmofrq.String);  % frequency smoothing of tapers, e.g., 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
% end
switch h.menu_preproc_ft_tfr_method.String{h.menu_preproc_ft_tfr_method.Value}
    case {'mtmfft'}
    case {'mtmconvol' }
        cfg.t_ftimwin   = repmat(str2num(h.edit_preproc_ft_tfr_t_ftimwin.String),1,length(foilim));  % frequency smoothing of tapers, e.g., 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
        cfg.toi   = sprintf('%s%%',h.edit_preproc_ft_tfr_toi.String);  
    case {'wavelet'}
        cfg.toi   = sprintf('%s%%',h.edit_preproc_ft_tfr_toi.String);  
    case {'tfr'}
        cfg.toi   = sprintf('%s%%',h.edit_preproc_ft_tfr_toi.String);  
    case {'mvar'}
end

if strcmp(cfg.output,'fourier')
    cfg.keeptrials  = 'yes';      % in order to separate the conditions again afterwards, we need to keep the trials. This is not otherwise necessary to compute the common filter
    cfg.keeptapers  = 'yes';
end

if isempty(h.current_fieldtrip_tfr_data)
    fprintf('Performing Time-Frequency Analysis of sensor data using "ft_freqanalysis.m"\n');
    ft_data = convert_bs2ft_data(h.sim_data.sens_final,h.anatomy,h.cfg);
    h.ft_data = ft_data;
    h.fieldtrip_tfr_data = ft_freqanalysis(cfg, ft_data);
else
    answ = questdlg('Would you like to overwrite current FieldTrip TFR data?','Data Exists!','yes','no','no');
    switch answ
        case 'yes'
            fprintf('Performing Time-Frequency Analysis of sensor data using "ft_freqanalysis.m"\n');
            ft_data = convert_bs2ft_data(h.sim_data.sens_final,h.anatomy,h.cfg);
            h.ft_data = ft_data;
            h.fieldtrip_tfr_data = ft_freqanalysis(cfg, ft_data);
        case 'no'
    end
    
end