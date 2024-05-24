function h = sm_batch_menu_sens_type_CallBack(h)

if h.menu_sens_type.Value == 1  % MEG
if h.menu_sens_montage.Value > length(h.anatomy.leadfield_meg_vol) ; h.menu_sens_montage.Value = 1; end
    
    switch h.menu_head_model.String{h.menu_head_model.Value}
        case 'Volume'
            h.anatomy.leadfield = h.anatomy.leadfield_meg_vol(h.menu_sens_montage.Value);
            h.anatomy.headmodel = h.anatomy.headmodel_meg_vol(h.menu_sens_montage.Value);
        case 'Cortical Surface'
            h.anatomy.leadfield = h.anatomy.leadfield_meg_cortex(h.menu_sens_montage.Value);
            h.anatomy.headmodel = h.anatomy.headmodel_meg_cortex(h.menu_sens_montage.Value);
    end
    
    if ~isempty(h.anatomy.leadfield)
        h.sensors_txt.String = sprintf('%.f sensors (%s) with %.f used for leadfield',size(h.anatomy.sens_meg(h.menu_sens_montage.Value).chanpos,1),h.anatomy.sens_meg(h.menu_sens_montage.Value).type,size(h.anatomy.leadfield.H,1));
        h.edit_leadfield_gain.String = '1e6';
        % find sensor labels of Leadfields that were used from sens_meg
        cfg.sens_idx=[];
        for v=1:length(h.anatomy.sens_meg(h.menu_sens_montage.Value).label)
            if sum(strcmp(h.anatomy.sens_meg(h.menu_sens_montage.Value).label(v),h.anatomy.leadfield.label))==1
                cfg.sens_idx = [cfg.sens_idx v];
            end
        end
    else
        cfg.sens_idx = find ( startsWith(h.anatomy.sens_meg(h.menu_sens_montage.Value).label,'M')==1) ;
        h.sensors_txt.String = sprintf('%.f sensors (%s) -- No Lead Field',size(h.anatomy.sens_meg(h.menu_sens_montage.Value).chanpos,1),h.anatomy.sens_meg(h.menu_sens_montage.Value).type);
    end
    h.edit_yscale.String = '-500 500'; h.edit_yscale_txt.String = 'Y Scale (fT)';
    cfg.sens_type = 'meg'; h.anatomy.sens = bs_select_meeg_sensors(cfg,h.anatomy.sens_meg(h.menu_sens_montage.Value));    % selecting only channels used in Leadfields for forward projecting data
    h.listbox_chans.String = h.anatomy.sens.label;
    
elseif h.menu_sens_type.Value == 2     % EEG
 if h.menu_sens_montage.Value > length(h.anatomy.leadfield_eeg_vol) ; h.menu_sens_montage.Value = 1; end
   
    switch h.menu_head_model.String{h.menu_head_model.Value}
        case 'Volume'
            h.anatomy.leadfield = h.anatomy.leadfield_eeg_vol(h.menu_sens_montage.Value);
            h.anatomy.headmodel = h.anatomy.headmodel_eeg_vol(h.menu_sens_montage.Value);
        case 'Cortical Surface'
            h.anatomy.leadfield = h.anatomy.leadfield_eeg_cortex(h.menu_sens_montage.Value);
            h.anatomy.headmodel = h.anatomy.headmodel_eeg_cortex(h.menu_sens_montage.Value);
    end
    
    if ~isempty(h.anatomy.leadfield)
        h.sensors_txt.String = sprintf('%.f sensors (%s) with %.f used for leadfield',size(h.anatomy.sens_eeg(h.menu_sens_montage.Value).chanpos,1),h.anatomy.sens_eeg(h.menu_sens_montage.Value).type,size(h.anatomy.leadfield.H,1));
        h.edit_yscale.String = '-100 100'; h.edit_yscale_txt.String = ['Y Scale (' char(181) 'V):'];
        
        % find sensor labels of Leadfields that were used from sens_eeg
        cfg.sens_idx=[];
        for v=1:length(h.anatomy.sens_eeg(h.menu_sens_montage.Value).label)
            if sum(strcmp(h.anatomy.sens_eeg(h.menu_sens_montage.Value).label(v),h.anatomy.leadfield.label))==1
                cfg.sens_idx = [cfg.sens_idx v];
            end
        end
    else
        cfg.sens_idx = find( strcmpi(h.anatomy.sens_eeg(h.menu_sens_montage.Value).chantype,'eeg')==1 );
        h.sensors_txt.String = sprintf('%.f sensors (%s) -- No Lead Field',size(h.anatomy.sens_eeg(h.menu_sens_montage.Value).chanpos,1),h.anatomy.sens_eeg(h.menu_sens_montage.Value).type);
    end
    
    h.edit_leadfield_gain.String = '1e-2';
    cfg.sens_type = 'eeg'; h.anatomy.sens = bs_select_meeg_sensors(cfg,h.anatomy.sens_eeg(h.menu_sens_montage.Value));    % selecting only channels used in Leadfields for forward projecting data
    h.listbox_chans.String = h.anatomy.sens.label;
    
elseif h.menu_sens_type.Value == 3      % MEEG  = combined MEG and EEG
    h.anatomy.sens = [];
end
