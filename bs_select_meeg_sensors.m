function sens_out = bs_select_meeg_sensors(cfg,sens)
% This function crates a new FieldTrip sensor structure with only those sensors selected by:
%  cfg.sens_idx = 1:66 
% cfg.sens_type = 'meg' or 'eeg'

sens_out = sens; 

switch cfg.sens_type
    case {'meg','MEG'}
        sens_out.chanori  = sens.chanori(cfg.sens_idx,:);
        sens_out.chanpos  = sens.chanpos(cfg.sens_idx,:);
        sens_out.chantype = sens.chantype(cfg.sens_idx);
        sens_out.chanunit = sens.chanunit(cfg.sens_idx);
        sens_out.coilori  = sens.coilori(cfg.sens_idx,:);
        sens_out.coilpos  = sens.coilpos(cfg.sens_idx,:);
        sens_out.elecpos  = sens.elecpos(cfg.sens_idx,:);
        sens_out.label = sens.label(cfg.sens_idx);
        sens_out.tra = sens.tra(cfg.sens_idx,cfg.sens_idx);
        sens_out.bad_sensors = setdiff(1:size(sens.chanpos,1),cfg.sens_idx);
        sens_out.good_sensors = cfg.sens_idx;
        sens_out.cfg=sens_out;
    case {'eeg','EEG'}
        sens_out.chanpos  = sens.chanpos(cfg.sens_idx,:);
        sens_out.chantype = sens.chantype(cfg.sens_idx);
        sens_out.chanunit = sens.chanunit(cfg.sens_idx);
        sens_out.elecpos  = sens.elecpos(cfg.sens_idx,:);
        sens_out.label = sens.label(cfg.sens_idx);
        sens_out.tra = sens.tra(cfg.sens_idx,cfg.sens_idx);
        sens_out.bad_sensors = setdiff(1:size(sens.chanpos,1),cfg.sens_idx);
        sens_out.good_sensors = cfg.sens_idx;
        sens_out.cfg=sens_out;
end


