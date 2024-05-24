function sm_create_leadfield(varargin)
global h

h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Calculating Lead Fields\n\n%s\n\nThis may take several minutes',h.anat_file); drawnow;

btn = questdlg(sprintf('Creating Lead Fields for all existing Head Models\n\nThis will overwrite these variables.\n\nDo you want to continue?'),'Create Lead Fields?','Yes','No','No');

switch btn
    case 'Yes'
        
        %% EEG Volume Whole-brain Forward Solution (leadfield)
        if ~isempty(h.anatomy.sens_eeg) && ~isempty(h.anatomy.headmodel_eeg_vol) && h.menu_sens_type.Value == 2 
            cfg            = [];
            %     cfg.elec       = h.anatomy.sens_eeg;
            cfg.elec       = h.anatomy.sens;     % This takes current sensors so it considers the reduced # channels
            cfg.channel   = {'eeg'}; %{'eeg' '-M1' '-M2'};    % this iw where you cna exlcude bad channels
            cfg.headmodel  = h.anatomy.headmodel_eeg_vol(h.menu_sens_montage.Value);
            cfg.headmodel.tissue = {'scalp' 'skull' 'brain'};
            cfg.resolution = 5;
            cfg.unit       = 'mm'; % same unit as above, i.e. in cm
            % cfg.normalize   = 'yes'; % to remove depth bias --> This improves reconstruction of signal amplitudes coming from deeper sources
            cfg.normalize   = 'no'; % to remove depth bias --> This improves reconstruction of signal amplitudes coming from deeper sources
            leadfield = ft_prepare_leadfield(cfg);
            % converting lead field into a more meaningful matrix representation of [chans x source locations].
            % only selecting inside voxels;
            x=cell2mat(leadfield.leadfield);
            leadfield.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
            leadfield.voxel_pos=leadfield.pos(leadfield.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
            leadfield_eeg_vol = leadfield;
            h.anatomy.leadfield_eeg_vol = leadfield;
            h.anatomy.leadfield = leadfield; h.menu_sens_type.Value = 2;
            
            % %     if ~strcmp(h.anat_file,'ANATOMY_DEFAULT_MEEG_adult_128chans.mat')
            %         [~,fname2,~]=fileparts(h.anat_file);
            %         sname = fullfile(h.anat_path,sprintf('%s.mat',fname2));
            %         save(sname,'leadfield_eeg_vol','-append');
            % %     end
        else
            warning('Lead Field for EEG not calculated. Please load EEG sensors and/or create Head Model');
        end
        %% EEG Cortex(Surface) Forward Solution (leadfield)
        if ~isempty(h.anatomy.sens_eeg) && ~isempty(h.anatomy.headmodel_eeg_cortex) && length(h.anatomy.mesh_volumes)>3 && h.menu_sens_type.Value == 2 
            cfg            = [];
            %     cfg.elec       = h.anatomy.sens_eeg;
            cfg.elec       = h.anatomy.sens;     % This takes current sensors so it considers the reduced # channels
            cfg.channel   = {'eeg'}; %{'eeg' '-M1' '-M2'};    % this iw where you cna exlcude bad channels
            cfg.headmodel  = h.anatomy.headmodel_eeg_cortex(h.menu_sens_montage.Value);
            cfg.sourcemodel.pos = h.anatomy.mesh_volumes(4).pos;
            % cfg.normalize   = 'yes'; % to remove depth bias --> This improves reconstruction of signal amplitudes coming from deeper sources
            cfg.normalize   = 'no'; % to remove depth bias --> This improves reconstruction of signal amplitudes coming from deeper sources
            leadfield = ft_prepare_leadfield(cfg);
            % converting lead field into a more meaningful matrix representation of [chans x source locations].
            % only selecting inside voxels;
            x=cell2mat(leadfield.leadfield');
            leadfield.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
            leadfield.voxel_pos=leadfield.pos(leadfield.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
            leadfield.ori = normals(leadfield.voxel_pos, h.anatomy.mesh_volumes(4).tri);
            leadfield_eeg_cortex = leadfield;
            h.anatomy.leadfield_eeg_cortex = leadfield;
            h.anatomy.leadfield = leadfield; h.menu_sens_type.Value = 2;
            
            % %     if ~strcmp(h.anat_file,'ANATOMY_DEFAULT_MEEG_adult_128chans.mat')
            %         [~,fname2,~]=fileparts(h.anat_file);
            %         sname = fullfile(h.anat_path,sprintf('%s.mat',fname2));
            %         save(sname,'leadfield_eeg_cortex','-append');
            % %     end
        else
            h.anatomy.leadfield_eeg_cortex = [];
            warning('Lead Field for EEG Cortex not calculated. Please load EEG sensors, create Head Model, or load in cortical surface as mesh_volume(4)');
        end
        
        %% MEG Volume Whole-brain Forward Solution (leadfield)
        if ~isempty(h.anatomy.sens_meg) && ~isempty(h.anatomy.headmodel_meg_vol) && h.menu_sens_type.Value == 1 
            % selecting only MEG sensors, excluding ref sensors
            cfg.sens_type = 'meg'; cfg.sens_idx = find(startsWith(h.anatomy.sens_meg.label,'M')==1);
            sens_meg = bs_select_meeg_sensors(cfg,h.anatomy.sens_meg);
            sens_meg = ft_convert_units(sens_meg,'mm'); 
            cfg            = [];
            cfg.grad       = sens_meg;
            cfg.channel   = {'meg'}; %{'eeg' '-M1' '-M2'};    % this is where you can exlcude bad channels
            cfg.headmodel  = h.anatomy.headmodel_meg_vol(h.menu_sens_montage.Value);
            cfg.headmodel = ft_convert_units(cfg.headmodel,'mm');
            cfg.resolution = 5;
            cfg.unit       = 'mm'; % same unit as above, i.e. in cm
            % cfg.normalize   = 'yes'; % to remove depth bias --> This improves reconstruction of signal amplitudes coming from deeper sources
            cfg.normalize   = 'no'; % to remove depth bias --> This improves reconstruction of signal amplitudes coming from deeper sources
            leadfield = ft_prepare_leadfield(cfg);
            % converting lead field into a more meaningful matrix representation of [chans x source locations].
            % only selecting inside voxels;
            x=cell2mat(leadfield.leadfield);
            leadfield.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
            leadfield.voxel_pos=leadfield.pos(leadfield.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
            leadfield_meg_vol = leadfield;
            h.anatomy.leadfield_meg_vol = leadfield;
            h.anatomy.leadfield = leadfield;  h.menu_sens_type.Value = 1;
            % %     if ~strcmp(h.anat_file,'ANATOMY_DEFAULT_MEEG_adult_128chans.mat')
            %         [~,fname2,~]=fileparts(h.anat_file);
            %         sname = fullfile(h.anat_path,sprintf('%s.mat',fname2));
            %         save(sname,'leadfield_meg_vol','-append');
            % %     end
        else
            warning('Lead Field for MEG not calculated. Please load MEG sensors and/or create Head Model');
        end
        %% MEG Cortex Whole-brain Forward Solution (leadfield)
        if ~isempty(h.anatomy.sens_meg) && ~isempty(h.anatomy.headmodel_meg_cortex) && length(h.anatomy.mesh_volumes)>3 && h.menu_sens_type.Value == 1 
            cfg.sens_type = 'meg'; cfg.sens_idx = find(startsWith(h.anatomy.sens_meg.label,'M')==1);
            sens_meg = bs_select_meeg_sensors(cfg,h.anatomy.sens_meg);
            sens_meg = ft_convert_units(sens_meg,'mm'); 
            cfg            = [];
            cfg.grad       = sens_meg;
            cfg.channel   = {'meg'}; %{'eeg' '-M1' '-M2'};    % this is where you can exlcude bad channels
            cfg.headmodel  = h.anatomy.headmodel_meg_cortex(h.menu_sens_montage.Value);
            cfg.headmodel = ft_convert_units(cfg.headmodel,'mm');
            cfg.sourcemodel.pos = h.anatomy.mesh_volumes(4).pos;
            % cfg.normalize   = 'yes'; % to remove depth bias --> This improves reconstruction of signal amplitudes coming from deeper sources
            cfg.normalize   = 'no'; % to remove depth bias --> This improves reconstruction of signal amplitudes coming from deeper sources
            leadfield = ft_prepare_leadfield(cfg);
            % converting lead field into a more meaningful matrix representation of [chans x source locations].
            % only selecting inside voxels;
            x=cell2mat(leadfield.leadfield');
            leadfield.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
%             leadfield.voxel_pos=leadfield.pos(leadfield.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
            leadfield.voxel_pos=leadfield.pos; % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
            leadfield.ori = normals(leadfield.voxel_pos, h.anatomy.mesh_volumes(4).tri);
%             leadfield_meg_cortex = leadfield;
            h.anatomy.leadfield_meg_cortex = leadfield;
            h.anatomy.leadfield = leadfield;  h.menu_sens_type.Value = 1;
            % %     if ~strcmp(h.anat_file,'ANATOMY_DEFAULT_MEEG_adult_128chans.mat')
            %         [~,fname2,~]=fileparts(h.anat_file);
            %         sname = fullfile(h.anat_path,sprintf('%s.mat',fname2));
            %         save(sname,'leadfield_meg_cortex','-append');
            % %     end
        else
            h.anatomy.leadfield_meg_cortex = [];
            warning('Lead Field for MEG Cortex not calculated. Please load MEG sensors, create Head Model, or load in cortical surface as mesh_volume(4)');
        end
        
    case 'No'
end

h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');





