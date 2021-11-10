


sens_data = h.sim_data.sens_noise_final;
% sens_data = h.sim_data.sens_sig_data;
leadfield = h.inv_soln(h.current_inv_soln).leadfield;

vx_idx = h.sim_data.cfg.source.vx_idx; 
source_ori = h.sim_data.cfg.source.vx_ori; 
%% conversion factor needed for EEG because Field Trips lead fields appear to be different than in Brainstorm (openmeeg) --> Not sure why this is?
% adjustment needed by user to get reasonable 150-200 fT for 35-40 nA dipoles in audiotry cortices like in Herdman et al., 2003 NeuroImage dataset
lf_gain = str2num(h.edit_leadfield_gain.String);
source_data = project_inverse_SimSignals(sens_data,leadfield,vx_idx,source_ori,lf_gain);

figure(2); clf; plot(squeeze(nanmean(h.sim_data.sig_final,3)),'k','linewidth',2); hold on; plot(squeeze(nanmean(source_data(:,:,:),3)));

