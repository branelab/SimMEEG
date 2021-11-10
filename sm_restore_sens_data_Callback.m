function sm_restore_sens_data_Callback(varargin)
% restores sensor data back to orginal data before preprocessing was applied
global h

if ~isfield(h.sim_data,'sens_final_org')
    h.sim_data.sens_final_org = h.sim_data.sens_final;
    h.sim_data.sens_noise_final_org = h.sim_data.sens_noise_final;
    h.sim_data.sens_sig_data_org = h.sim_data.sens_sig_data;
elseif isfield(h.sim_data,'sens_final_org')
    h.sim_data.sens_final = h.sim_data.sens_final_org;
    h.sim_data.sens_noise_final = h.sim_data.sens_noise_final_org;
    h.sim_data.sens_sig_data = h.sim_data.sens_sig_data_org;
end