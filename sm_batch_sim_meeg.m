function h = sm_batch_sim_meeg(h)
% project simulated data

%% removing original sens data
if isfield(h.sim_data,'sens_final_org')
    h.sim_data = rmfield(h.sim_data,'sens_final_org');
    h.sim_data = rmfield(h.sim_data,'sens_noise_final_org');
    h.sim_data = rmfield(h.sim_data,'sens_sig_data_org');
end

%% checking if source data exist
h.sim_data.cfg = h.monte_params.cfg;
if ~isfield(h.sim_data,'sig_final') % source data have not yet been simulated
    src.Tag=''; hobj=''; update_cfg(src,hobj); run_sim(src,hobj);    % simulating source data
end

%% projecting source sig_final to sensors
leadfield = h.anatomy.leadfield;
h.sim_data.sens_sig_data = sm_batch_project_SimSignals(h.sim_data.sig_final,leadfield,h.monte_params.cfg.source.vx_idx,h.monte_params.cfg.source.vx_amp,h.monte_params.cfg.source.vx_ori,h);
h = sm_batch_combine_sens_data(h);
fprintf('Finished Simulating M/EEG sensor data\n');
