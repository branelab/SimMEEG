function sens_out = sm_make_sens_compatible(sens,sens_type)


sens_out = sens; 


switch sens_type
    case {'meg','MEG'}
       chan_unit = 'T'; senstype = 'meg';
        
    case {'eeg','EEG'}
        chan_unit = 'V'; senstype = 'eeg';
end


if ~isfield(sens_out,'label'); sens_out.label = []; end
if ~isfield(sens_out,'type'); sens_out.type = senstype; end
if ~isfield(sens_out,'chantype'); sens_out.chantype = repmat({sens_type},1,length(sens_out.label)); end
if ~isfield(sens_out,'chanunit'); sens_out.chanunit = repmat({chan_unit},1,length(sens_out.label)); end
if ~isfield(sens_out,'unit'); sens_out.type = 'm'; end
if ~isfield(sens_out,'chanpos'); sens_out.chanpos = []; end
% if ~isfield(sens_out,'chanori'); sens_out.chanori = []; end     % only for meg sensors
% if ~isfield(sens_out,'coilori'); sens_out.coilori = []; end     % only for meg sensors
% if ~isfield(sens_out,'coilpos'); sens_out.coilpos = []; end     % only for meg sensors
if ~isfield(sens_out,'elecpos'); sens_out.elecpos = sens_out.chanpos; end     % only for meg sensors
if ~isfield(sens_out,'tra'); sens_out.tra = []; end
meg_idx = strcmp('megaxial',sens_out.chantype(:));  % need to convert to 'meg' to work with ft_plot_sens.m
sens_out.chantype(meg_idx) = {'meg'}; 
if ~isfield(sens_out,'good_sensors'); sens_out.good_sensors = 1:length(sens_out.label); end
if ~isfield(sens_out,'bad_sensors'); sens_out.bad_sensors = []; end

% sens_out.tra = sens_out.tra(length(sens_out.label),length(sens_out.label));

      
        