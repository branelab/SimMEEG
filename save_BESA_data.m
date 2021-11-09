function save_BESA_data(varargin)
global h

[fname,fpath]=uiputfile( '*.avr','Save BrainSim Study',sprintf('%s_BESA_data_%s',h.edit_study_name.String,h.menu_sens_type.String{h.menu_sens_type.Value}) );
% hm = msgbox(sprintf('Writing to BESA (.avr) format %s\n Please wait...',fname));
h.waitfor_panel.Visible='on'; h.waitfor_txt.String =sprintf('Writing to BESA (.avr) format %s\n Please wait...',fname); drawnow;

data_matrix = permute(squeeze(nanmean(h.sim_data.sens_final,3)),[2 1]); % only saving averaged data right now
time_samples = h.cfg.study.lat_sim;
channel_labels = h.anatomy.sens.label;
data_scale_factor = 1; % saving units are in uV
time_scale_factor = 1000; % saving units are in ms so need to convert sec-->ms

status = besa_save2Avr(fpath, fname, data_matrix, ...
    time_samples, channel_labels, data_scale_factor, time_scale_factor);
% close(hm);
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
