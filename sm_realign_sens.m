function sm_realign_sens(varargin)
global h;

if h.menu_sens_type.Value == 1      %MEG
    msgbox('Changing MEG sensors is not available yet.');
elseif h.menu_sens_type.Value == 2      %EEG
    
    elec_aligned = h.anatomy.sens_eeg;
    
    %% re-aligning electrodes
    hm=msgbox('Align EEG channels to scalp');
    cfg           = [];
    cfg.method    = 'interactive';
    cfg.elec      = elec_aligned;
    cfg.headshape = h.anatomy.mesh_volumes(1);    % scalp mesh
    elec_aligned  = ft_electroderealign(cfg,elec_aligned);
    if isvalid(hm); close(hm); end
    %% project to surface
    cfg           = [];
    cfg.method    = 'project';
    cfg.elec      = elec_aligned;
    cfg.headshape = h.anatomy.mesh_volumes(1);    % scalp mesh
    elec_aligned2  = ft_electroderealign(cfg,elec_aligned);
    %% double check alignment
    hm=msgbox('Check Projected Channels');
    cfg           = [];
    cfg.method    = 'interactive';
    cfg.elec      = elec_aligned2;
    cfg.headshape = h.anatomy.mesh_volumes(1);    % scalp mesh
    elec_aligned2  = ft_electroderealign(cfg,elec_aligned2);
    if isvalid(hm); close(hm); end
    
    %         h.anatomy.elec_aligned = elec_aligned2;
    %% Selecting only EEG good channels for creating leadfields
    chan_size = {sprintf('1:%.f',length(elec_aligned2.label))};
    chan_idx = inputdlg({'Enter Good EEG Scalp Channels index:'},'Select Good EEG Scalp Channels',1,chan_size);
    cfg.sens_type = 'eeg'; cfg.sens_idx = str2num(chan_idx{:});
    elec_good = bs_select_meeg_sensors(cfg,elec_aligned2);
    %% double check alignment
    hm=msgbox('Double check alignment of selected channels');
    cfg           = [];
    cfg.method    = 'interactive';
    cfg.elec      = elec_good;
    cfg.headshape = h.anatomy.mesh_volumes(1);    % scalp mesh
    elec_good  = ft_electroderealign(cfg,elec_good);
    h.anatomy.elec_good = elec_good;
    sens_eeg = elec_good;
    h.anatomy.sens_eeg = elec_good;
    if isvalid(hm); close(hm); end
    %% reset channel location table
    h.sens_table.Data(:,1) = h.anatomy.sens_eeg.label;
    h.sens_table.Data(:,2:4) = num2cell(h.anatomy.sens_eeg.chanpos,[size(h.anatomy.sens_eeg.chanpos,1) size(h.anatomy.sens_eeg.chanpos,2)]);
end


