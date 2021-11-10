function sm_load_sensors(varargin)

global h


btn = questdlg(sprintf('Load Channels or Sensors\n\nThis will overwrite these variables.\n\nDo you want to continue?',h.anat_file),'Load Sensors?','Yes','No','No');

switch btn
    case 'Yes'
        sens_type = questdlg(sprintf('Which Type?'),'Load Sensors?','MEG','EEG','MEG');
%         sens_type = h.menu_sens_type.String{h.menu_sens_type.Value};
        switch sens_type
            case 'EEG'
                % user selects EEG electrode file
                [fname,fpath]=uigetfile('*.*','Load EEG Electrodes from File');
                h.anatomy.elec_file = fname;
                %% load EEG sensors
                elec = ft_read_sens(fullfile(fpath,fname));
                %% convert to mm
                [elec] = ft_convert_units(elec, 'mm');
                
                h.anatomy.elec = elec;
                h.anatomy.elec_org = elec;
                elec_aligned = elec;
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
%                 %% save sensors in anatomy .mat file
%                 %         if ~strcmp(h.anat_file,'ANATOMY_DEFAULT_MEEG_adult.mat')
%                 [~,fname2,~]=fileparts(h.anat_file);
%                 sname = fullfile(h.anat_path, sprintf('%s.mat',fname2));
%                 if exist(sname,'file')
%                     save(sname,'elec','elec_good','sens_eeg','-append');
%                 else
%                     save(sname,'elec','elec_good','sens_eeg');
%                 end
%                 %         end
                
            case 'MEG'
                [fname]=uigetdir('Load MEG Electrodes from dataset File');
                %% Load MEG sensors
                grad = ft_read_sens(fname, 'senstype', 'meg');
                %% convert to mm
                [grad] = ft_convert_units(grad, 'mm');
                grad_aligned = grad;
                %% making sure num coilpos = num labels
                grad_aligned.coilori = grad_aligned.coilori(1:size(grad_aligned.label,1),:);
                grad_aligned.coilpos = grad_aligned.coilpos(1:size(grad_aligned.label,1),:);
                grad_aligned.tra     = grad_aligned.tra(1:size(grad_aligned.label,1),1:size(grad_aligned.label,1));
                grad_aligned.elecpos = grad_aligned.chanpos;
                %% re-aligning electrodes
                hm=msgbox('Align MEG sensors to scalp');
                cfg           = [];
                cfg.method    = 'interactive';
                cfg.keepchannel  = 'yes';
                cfg.elec      = grad_aligned;
                cfg.headshape = h.anatomy.mesh_volumes(1);    % scalp mesh
                [grad_aligned2]  = ft_electroderealign(cfg,grad_aligned);
                close(hm);
                %% Selecting only MEG good channels for creating leadfields
                chan_size = {sprintf('1:%.f',length(grad_aligned2.label))};
                chan_idx = inputdlg({'Enter Good MEG Channels index:'},'Select Good MEG Channels',1,chan_size);
                cfg.sens_type = 'meg'; cfg.sens_idx = str2num(chan_idx{:});
                grad_good = bs_select_meeg_sensors(cfg,grad_aligned2);
                %% double check alignment
                hm=msgbox('Double check alignment of selected sensors');
                cfg           = [];
                cfg.method    = 'interactive';
                cfg.elec      = grad_good;
                cfg.headshape = h.anatomy.mesh_volumes(1);    % scalp mesh
                grad_good  = ft_electroderealign(cfg,grad_good);
                h.anatomy.grad_good = grad_good;
                sens_meg = grad_good;
                h.anatomy.sens_meg = grad_good;
                close(hm);
%                 %% save sensors in anatomy .mat file
%                 %         if ~strcmp(h.anat_file,'ANATOMY_DEFAULT_MEEG_adult.mat')
%                 [~,fname2,~]=fileparts(h.anat_file);
%                 sname = fullfile(h.anat_path, sprintf('%s.mat',fname2));
%                 if exist(sname,'file')
%                     save(sname,'grad','grad_good','sens_meg','-append');
%                 else
%                     save(sname,'grad','grad_good','sens_meg');
%                 end
%                 %         end
        end

update_anatomy_fields();
    case 'No'
end
