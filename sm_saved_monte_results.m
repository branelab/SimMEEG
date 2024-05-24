function sm_saved_monte_results(varargin)
global h

sname = varargin{end};
save(sname,'-struct', 'h', 'cfg','monte_params'); % default variables to save
study_name = h.cfg.study.study_name; save(sname,'study_name','-append'); % default

%% save "sim_data.sens_final" also as FieldTrip's Data format
if h.radio_monte_save_FT_data.Value == 1
    ft_data = convert_bs2ft_data(h.sim_data.sens_final,h.anatomy,h.cfg);
    save(sname,'ft_data','-append');
end

%% saving selected results
clear sim_data;
%% sim_data
fnames = h.listbox_monte_saved_sim_data.String(h.listbox_monte_saved_sim_data.Value);
for f = 1:length(fnames)
try    fn = fnames(f); sim_data.(fn{1}) = h.sim_data.(fn{1}); end
end

%% true source
fnames = h.listbox_monte_saved_true_source.String(h.listbox_monte_saved_true_source.Value);
if isfield(sim_data,'cfg')
if isfield(sim_data.cfg,'source'); sim_data.cfg = rmfield(sim_data.cfg,'source'); end % resetting cfg below with selected fieldnames
end

for f = 1:length(fnames)
    try fn = fnames(f); sim_data.cfg.source.(fn{1}) = h.sim_data.cfg.source.(fn{1}); end
end

if isfield(sim_data.cfg.source,'TFR_results')
    fnames = h.listbox_monte_saved_true_source_TFR_results.String(h.listbox_monte_saved_true_source_TFR_results.Value);
        sim_data.cfg.source = rmfield(sim_data.cfg.source,'TFR_results'); % resetting cfg below with selected fieldnames
    for f = 1:length(fnames)
        try fn = fnames(f); sim_data.cfg.source.TFR_results.(fn{1}) = h.sim_data.cfg.source.TFR_results.(fn{1}); end
    end
end
save(sname,'sim_data','-append'); 

%% inv_soln
fnames = h.listbox_monte_saved_inv_soln.String(h.listbox_monte_saved_inv_soln.Value)';
fnames{2,1} = {};
inv_soln = struct(fnames{:});
for a=1:length(h.inv_soln)
    for f = 1:length(fnames)
       try fn = fnames(1,f); inv_soln(a).(fn{1}) = h.inv_soln(a).(fn{1}); end
    end
end
fnames_soln = h.listbox_monte_saved_inv_soln_soln.String(h.listbox_monte_saved_inv_soln_soln.Value)';
if isfield(inv_soln,'soln')
    inv_soln = rmfield(inv_soln,'soln'); %inv_soln(1).soln = struct(fnames_soln{:});
    fnames_soln{2,1} = {};
    for a=1:length(h.inv_soln)
        for f = 1:length(fnames_soln)
            try fn = fnames_soln(1,f); inv_soln(a).soln.(fn{1}) = h.inv_soln(a).soln.(fn{1}); end
        end
    end
end


fnames_tfr = h.listbox_monte_saved_inv_soln_TFR_results.String(h.listbox_monte_saved_inv_soln_TFR_results.Value)';
if isfield(inv_soln,'TFR_results')
    inv_soln = rmfield(inv_soln,'TFR_results'); %inv_soln(1).soln = struct(fnames_soln{:});
    fnames_tfr{2,1} = {};
    for a=1:length(h.inv_soln)
        for f = 1:length(fnames_tfr)
            try fn = fnames_tfr(1,f); inv_soln(a).TFR_results.(fn{1}) = h.inv_soln(a).TFR_results.(fn{1}); end
        end
    end
end


save(sname,'inv_soln','-append');
