function data_int = sm_interpolate_data(data,bad_chans)
% UNDER CONSTRUCTION - currently not implemented because "Warning: inconsistent sampleinfo" when running --> [neighbours] = ft_prepare_neighbours(cfg,ft_data);
% This program uses Fiel Trips programs for interpolated bad channels that have been set as "nans"
global h


%% Check to see if channel layout exists for current anatomy

if ~isfield(h.anatomy.sens,'layout') && strcmpi(h.anatomy.sens.type,'meg')
cfg = [];
cfg.grad = h.anatomy.sens;
cfg.projection = 'orthographic';
cfg.channel = 'M';
h.anatomy.sens.layout = ft_prepare_layout(cfg); 

elseif ~isfield(h.anatomy.sens,'layout') && strcmpi(h.anatomy.sens.type,'eeg')
cfg = [];
cfg.elec = h.anatomy.sens;
cfg.projection = 'orthographic';
cfg.channel = 'all';
h.anatomy.sens.layout = ft_prepare_layout(cfg); 

end

%% Check to see if neighbours calculated for sensor positions
% convert to FT data 
ft_data = convert_bs2ft_data(data,h.anatomy,h.cfg);
cfg = [];
cfg.method = 'distance'; %, 'triangulation' or 'template' 
cfg.neighbourdist = 4; % number, maximum distance between neighbouring 
cfg.layout = h.anatomy.sens.layout;
[neighbours] = ft_prepare_neighbours(cfg,ft_data);

