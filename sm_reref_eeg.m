function sm_reref_eeg(varargin)
global h

if ~isempty(h.sim_data)
    if ~isempty(h.sim_data.sens_final) && h.btn_reref_EEG.Value ==1 % REREF
        dims = size(h.sim_data.sens_final);
%         ref_chans = 1:size(h.anatomy.sens.chanpos,1);
        ref_chans = h.listbox_chans.Value; 
        ex_chans = setxor(1:dims(2), ref_chans);
        if ~isfield(h.sim_data, 'sens_final_org')
            h.sim_data.sens_final_noref = h.sim_data.sens_final;
        end
        [data]=bl_reref(h.sim_data.sens_final,ref_chans,ex_chans,1);
        h.sim_data.sens_final = data; 
        
        %% update btn
        if length(ref_chans)<=2
            ref_names = deblank(sprintf('%s ', h.anatomy.sens.label{ref_chans})); 
        elseif length(ref_chans)>2 && length(ref_chans) < length(h.listbox_chans.String)
             ref_names = sprintf('%s ...', deblank(sprintf('%s ', h.anatomy.sens.label{ref_chans(1:2)}))); 
       elseif length(ref_chans) == length(h.listbox_chans.String)
            ref_names = 'Avg Ref';
        end
            h.btn_reref_EEG.String = sprintf('Ref: %s',ref_names);
        h.btn_reref_EEG.BackgroundColor = [.9 1 .9];

    else    % revert back to non-ref
        h.sim_data.sens_final = h.sim_data.sens_final_noref;
        h.btn_reref_EEG.String = 'Ref: None';
        h.btn_reref_EEG.BackgroundColor = [1 .9 .9];
    end
    h.fcn_plot_sens_data();
end

