function sm_create_headmodel(varargin)
global h

h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Calculating LeadFields\n\n%s\n\nThis may take several minutes',h.anat_file); drawnow;

btn = questdlg(sprintf('Head Models for EEG and MEG are only created for the following:\n\nMEG = localspheres (Volume)\nEEG = concentricspheres (Volume)\n\nThis will overwrite these variables.\n\nDo you want to continue?'),'Create Head Models?','Yes','No','No');

switch btn
    case 'Yes'
        if ~isempty(h.anatomy.sens_eeg)
            %% Head model
            cfg = [];
            cfg.tissue = {'scalp' 'skull' 'brain'};
            cfg.conductivity = [0.3300, 0.04, 0.3300];   % based on Lew et al. NeuroImage 76 (2013) 282-293
            cfg.elec = h.anatomy.sens_eeg;
            cfg.method = 'concentricspheres'; % 'bemcp'; % 'openmeeg'; % 'concentricspheres'; %'singlesphere'; % 'bemcp'; %
            % headmodel_eeg_skull = ft_prepare_headmodel(cfg, mesh_volumes([1 2 3]));
            %      cfg.fitind = 3; % fit center of sphere to brain hull
            %    headmodel_eeg_vol = ft_prepare_headmodel(cfg, h.anatomy.mesh_volumes(1:3));
            headmodel_eeg_vol = ft_prepare_headmodel(cfg, h.anatomy.mesh_volumes(2));
            origin = headmodel_eeg_vol.o;
            radius = headmodel_eeg_vol.r;
            %% re-adjusting sphere to fit brain hull range in the anterior-posterior direction to capture more brain area
            maxradius = range(h.anatomy.mesh_volumes(2).pos(:,1))/2;
            expansion_gain = maxradius/max(radius);
            headmodel_eeg_vol.r = headmodel_eeg_vol.r * expansion_gain;
            rx = max(radius - headmodel_eeg_vol.r);
            headmodel_eeg_vol.o(3) = headmodel_eeg_vol.o(3)+rx;
            figure(999); clf; hold on;
            ft_plot_mesh(h.anatomy.mesh_volumes(1),'Facecolor','none','EdgeAlpha',.1,'edgecolor','k');
            ft_plot_mesh(h.anatomy.mesh_volumes(2),'Facecolor','b','EdgeColor','none','FaceAlpha',.1);
            ft_plot_headmodel(headmodel_eeg_vol,'Facecolor','r','FaceAlpha',.1); view(180,0);
            title('EEG Head Model Concentric Spheres');
            legend({'scalp' 'brain hull' 'headmodel'});
            
            h.anatomy.headmodel_eeg_vol = headmodel_eeg_vol;
            h.anatomy.headmodel = h.anatomy.headmodel_eeg_vol;
            h.menu_sens_type.Value = 2;
            
%             %% save anatomy .mat file
%             % if ~strcmp(h.anat_file,'ANATOMY_DEFAULT_MEEG_adult_128chans.mat')
%             [~,fname2,~]=fileparts(h.anat_file);
%             sname = fullfile(h.anat_path,sprintf('%s.mat',fname2));
%             try
%                 save(sname,'headmodel_eeg*','-append');
%             catch
%                 save(sname,'headmodel_eeg*');
%             end
%             % end
            
        end
        
        if ~isempty(h.anatomy.sens_meg)
            %% Head model
            % adult dewar
            h99 = figure(99);
            cfg = [];
            cfg.tissue = {'scalp' 'skull' 'brain'}; cfg.conductivity = [0.3300, 0.04, 0.3300];   % based on Lew et al. NeuroImage 76 (2013) 282-293
            cfg.grad = h.anatomy.sens_meg;
            cfg.method = 'localspheres'; % 'openmeeg'; %'concentricspheres'; %'singlesphere'; % 'bemcp'; %
            headmodel_meg = ft_prepare_headmodel(cfg,h.anatomy.mesh_volumes(3));
            headmodel_meg_vol = headmodel_meg;
            
            close(h99);
            
            % save headmodel_eeg headmodel_eeg;
            h.anatomy.headmodel_meg_vol = headmodel_meg_vol;
            h.anatomy.headmodel = h.anatomy.headmodel_meg_vol;
            h.menu_sens_type.Value = 1;
            
%             %% save anatomy .mat file
%             [~,fname2,~]=fileparts(h.anat_file);
%             sname = fullfile(h.anat_path,sprintf('%s.mat',fname2));
%             try
%                 save(sname,'headmodel_meg*','-append');
%             catch
%                 save(sname,'headmodel_meg*');
%             end
            
        end
        
        update_anatomy_fields();
    case 'No'
end

h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');

