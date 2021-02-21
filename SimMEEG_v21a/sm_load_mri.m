function sm_load_mri(varargin)
global h

h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Loading MRI\n\n%s\n\nThis may take several minutes to segment the mri',h.anat_file); drawnow;
btn = questdlg(sprintf('Using Field Trip to load MRI and segmentation\n\nThis will overwrite these variables.'),'Load MRI?','Yes','No','No');

switch btn
    case 'Yes'
        
        [fname,fpath]=uigetfile({'*.m*','MRI File (*.mri, or *.mat)'; '*.*','All files'});
        h.anat_path = fpath;
        h.anat_file = fname;
        
        [~,fname2,fext]=fileparts(h.anat_file);
        
        switch fext
            
            case '.mri' % if .mri then load it, segment it, and create mesh_volumes
                % clearing anatomy
                h.anatomy = struct('mri','','sens_eeg','','sens_meg','','mesh_volumes',[],'headmodel_eeg_vol','','headmodel_meg_vol','','headmodel_eeg_cortex','','headmodel_meg_cortex','','leadfield_eeg_vol','','leadfield_eeg_cortex','','leadfield_meg_vol','','leadfield_meg_cortex','','sens','','headmodel','','leadfield','','mesh_cortex','');
                
                %% following FieldTrip's Tutorial on generating Headmodel @ https://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_bem/
                %% Read MRI
                mri = ft_read_mri(fullfile(fpath,fname));
                h.anatomy.mri = mri;
                
                seg_file = fullfile(h.anat_path,'segmentedmri.mat');
                if exist(seg_file,'file')   % load previously segmented mri
                    fprintf('Loading "segmentedmri" from file: %s\n',seg_file);
                    load(seg_file);
                else % segment MRI
                    %% Segmentation of MRI
                    cfg           = [];
                    cfg.output    = {'brain','skull','scalp'};
                    segmentedmri = ft_volumesegment(cfg, h.anatomy.mri);
                    h.anatomy.segmentedmri = segmentedmri;
                end
                
                %% Mesh
                cfg=[];
                cfg.tissue={'scalp','skull','brain'};
                cfg.numvertices = [3000 2000 1000];
                mesh_volumes = ft_prepare_mesh(cfg,segmentedmri);
                % reorder so that it is
                h.anatomy.mesh_volumes = mesh_volumes;
                h.anatomy.mesh_volumes_org = mesh_volumes;
                
%                 %% save anatomy .mat file
%                 %         if ~strcmp(h.anat_file,'ANATOMY_DEFAULT_MEEG_adult.mat')
%                 sname = fullfile(h.anat_path,sprintf('%s.mat',fname2));
%                 if exist(sname,'file')
%                     save(sname,'mri','segmentedmri','mesh_volumes','-append');
%                 else
%                     save(sname,'mri','segmentedmri','mesh_volumes');
%                 end
%                 %         end
                
                %% plot mesh
                figure; clf; hold on;
                ft_plot_mesh(h.anatomy.mesh_volumes(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
                ft_plot_mesh(h.anatomy.mesh_volumes(2),'edgecolor','none','facealpha',0.4);
                ft_plot_mesh(h.anatomy.mesh_volumes(3),'edgecolor','none','facecolor',[0.4 0.6 0.4]); view(180,0);
        end
        update_anatomy_fields();
    case 'No'
end

h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');


