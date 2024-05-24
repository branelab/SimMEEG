function sm_set_ARM_params(varargin)
global h

%% showing all ARM params
for a=1:length(h.panel_ARM_params.Children); h.panel_ARM_params.Children(a).Visible = 'on'; end
h.axes_ARM_interactions.clo;


%% Setting Visible based on "Combine Source Data" Menu input
if h.menu_ARM_add.Value == 1    % hiding all ARM parameters
    for a=1:length(h.panel_ARM_params.Children); h.panel_ARM_params.Children(a).Visible = 'off'; end
    h.menu_ARM_add_txt.Visible = 'on'; h.menu_ARM_add.Visible = 'on';

elseif h.menu_ARM_add.Value == 3    % adding waveforms together then only 3 sources are allowed
    % set num_sources = 3 and disable edit box
    h.edit_ARM_num_sources.String = '3';
    h.edit_ARM_num_sources.Enable = 'off';
    h.menu_ARM_ori_constraint_txt.Visible = 'off';
    h.menu_ARM_ori_constraint.Visible = 'off';
    h.menu_ARM_source_locs_txt.Visible = 'off';
    h.menu_ARM_source_locs.Visible = 'off';
    % limiting number of interaction to those possible
    %% Maximum number of interactions can't excced num_sources^2
    n_sources = str2num(h.edit_ARM_num_sources.String);
    max_int = n_sources^2 - n_sources;
    if str2num(h.edit_ARM_num_interactions.String)>max_int
        h.edit_ARM_num_interactions.String = num2str(max_int); % possible interactions are all possible combinations because ARM model simulates symmetrical and assymetrical interactions meaning that asymetrical will have interactions for source1-->source 2 but not source2-->source1
    end
    
    varargin(end) = {'Add Waveforms'};
else
    h.edit_ARM_num_sources.Enable = 'on';
end

%% setting specific called functions
switch varargin{end}
    case 'Num Interactions'
        %% Maximum number of interactions can't excced num_sources^2
        n_sources = str2num(h.edit_ARM_num_sources.String);
        max_int = n_sources^2 - n_sources;
        if str2num(h.edit_ARM_num_interactions.String)>max_int
            h.edit_ARM_num_interactions.String = num2str(max_int); % possible interactions are all possible combinations because ARM model simulates symmetrical and assymetrical interactions meaning that asymetrical will have interactions for source1-->source 2 but not source2-->source1
        end
    case 'Ori Constraint'
        sm_ARM_set_source_ori(h.menu_ARM_ori_constraint.String{h.menu_ARM_ori_constraint.Value});    % setting orientations
        sm_ARM_change_source_params('Update');
    case 'Locations'
        if h.menu_ARM_source_locs.Value == 1 % Randomly selecting sources locations within leadfield grid with a minimum spatial distance set by h.edit_ARM_min_spatial_dist
            num_sources = str2num(h.edit_ARM_num_sources.String);
            dist_thresh = str2num(h.edit_ARM_min_spatial_dist.String);
            rn = randperm(size(h.anatomy.leadfield.H,3));
            seed_idx = rn(1);   % first voxel as seed then find other rnadom voxels that are > dist_thresh
            vx_pos = h.anatomy.leadfield.voxel_pos;
            k=0;
            hm = waitbar(0,sprintf('Number of Distant Sources Found (Total = %.f)',num_sources));
           while k<1000
                k=k+1;
                rn = randperm(size(h.anatomy.leadfield.H,3)); 
%                 rn_diff = setdiff(rn,vx_idx,'stable'); 
%                 rn = [vx_idx rn_diff];
                [vdist, ~]= pdist2(vx_pos(rn(1:num_sources),:),vx_pos(rn(1:num_sources),:),'euclidean','Smallest',2); % finds smallest Euclidean distane among sources
                waitbar(sum(vdist(2,:)>dist_thresh)/num_sources, hm)
                    num_found = sum(vdist(2,:)>dist_thresh);
               fprintf('Distant Sources Found = %.f\n',num_found);
                    %% break loop
                 if  num_found == num_sources
                     vx_idx = rn(1:num_sources);
                      break; 
                 end
                 if k>=1000 
                     vx_idx = rn(1:num_sources);
                     warndlg(sprintf('Only %.f/%.f sources were found to be %.f mm distant from each other\n\n Try reducing "Number Sources" or "Min Distance"',num_found,num_sources,dist_thresh)); 
                 end
           end
           
            close(hm); 
            h.cfg.ARM_params.vx_idx = vx_idx;
            h.cfg.ARM_params.vx_locs = h.anatomy.leadfield.voxel_pos(h.cfg.ARM_params.vx_idx,:);
            
            sm_ARM_set_source_ori(h.menu_ARM_ori_constraint.String{h.menu_ARM_ori_constraint.Value});    % setting orientations
            h.cfg.ARM_params.vx_amp = randi(str2num(h.edit_ARM_amp_range.String),1,num_sources);
            h.cfg.ARM_params.sig_amp_perc = randi(str2num(h.edit_ARM_sig_perc_range.String),1,num_sources); % place holder for now
            h.cfg.ARM_params.prepost_amp_perc = randi(str2num(h.edit_ARM_prepost_perc_range.String),1,num_sources); % place holder for now
            h.cfg.ARM_params.sig_latency = repmat(str2num(h.edit_ARM_sig_latency.String),num_sources,1);
            h.cfg.ARM_params.sig_risetime = repmat(str2num(h.edit_ARM_sig_risetime.String),num_sources,1);
            h.cfg.ARM_params.source_locs_file = '';
            h.edit_ARM_num_sources.Enable = 'on';
            
            
        elseif h.menu_ARM_source_locs.Value == 2 % load source locations from file]
            [fname,fpath] = uigetfile({'*.csv;*.xlsx','CSV XLS'},'Load Source Locations File');
            h.cfg.ARM_params.source_locs_file = fullfile(fpath,fname);
            if any(fname~=0)
                [~,~,fext] = fileparts(fname);
                switch fext
                    case '.csv'
                        [vx_locs] = textread(fullfile(fpath,fname),'','delimiter',',');
                    case '.xlsx'
                        [vx_locs] =  xlsread(fullfile(fpath,fname));
                end
                % in Source file: clmns 1:3 = source locations, clmn 4 = source amplitudes (nAm)
                h.cfg.ARM_params.vx_idx = find_nearest_voxel(vx_locs(:,1:3),h.anatomy.leadfield.voxel_pos);
                h.cfg.ARM_params.vx_locs = h.anatomy.leadfield.voxel_pos(h.cfg.ARM_params.vx_idx,:);
                if size(h.cfg.ARM_params.vx_locs,2) == 3  % only locations --> randomly assign loaded orientations and amplitudes
                    h.cfg.ARM_params.vx_amp = randi(str2num(h.edit_ARM_amp_range.String),1,num_sources);
                    h.cfg.ARM_params.sig_amp_perc = randi(str2num(h.edit_ARM_sig_perc_range.String),1,num_sources); % place holder for now
                    h.cfg.ARM_params.prepost_amp_perc = randi(str2num(h.edit_ARM_prepost_perc_range.String),1,num_sources); % place holder for now
                    h.cfg.ARM_params.sig_latency = repmat(str2num(h.edit_ARM_sig_latency.String),num_sources,1);
                    h.cfg.ARM_params.sig_risetime = repmat(str2num(h.edit_ARM_sig_risetime.String),num_sources,1);
                    sm_ARM_set_source_ori('Random');
                elseif  size(h.cfg.ARM_params.vx_locs,2) == 4   % only locations & amplitudes --> randomly assign loaded orientations
                    h.cfg.ARM_params.vx_amp = vx_locs(:,4)';
                    h.cfg.ARM_params.sig_amp_perc = randi(str2num(h.edit_ARM_sig_perc_range.String),1,num_sources); % place holder for now
                    h.cfg.ARM_params.prepost_amp_perc = randi(str2num(h.edit_ARM_prepost_perc_range.String),1,num_sources); % place holder for now
                    h.cfg.ARM_params.sig_latency = repmat(str2num(h.edit_ARM_sig_latency.String),num_sources,1);
                    h.cfg.ARM_params.sig_risetime = repmat(str2num(h.edit_ARM_sig_risetime.String),num_sources,1);
                    sm_ARM_set_source_ori('Random');
                elseif  size(h.cfg.ARM_params.vx_locs,2) == 6    % only locations & loaded orientations  --> randomly assign amplitudes
                    h.cfg.ARM_params.vx_amp = randi([30 100],1,size(vx_locs,1));
                    h.cfg.ARM_params.vx_ori = vx_locs(:,4:6);
                    h.cfg.ARM_params.sig_amp_perc = randi(str2num(h.edit_ARM_sig_perc_range.String),1,num_sources); % place holder for now
                    h.cfg.ARM_params.prepost_amp_perc = randi(str2num(h.edit_ARM_prepost_perc_range.String),1,num_sources); % place holder for now
                    h.cfg.ARM_params.sig_latency = repmat(str2num(h.edit_ARM_sig_latency.String),num_sources,1);
                    h.cfg.ARM_params.sig_risetime = repmat(str2num(h.edit_ARM_sig_risetime.String),num_sources,1);
               elseif size(h.cfg.ARM_params.vx_locs,2) == 7 % Locations & loaded orientations  % amplitude in File
                    h.cfg.ARM_params.vx_amp = vx_locs(:,size(h.cfg.ARM_params.vx_locs,2))';
                    h.cfg.ARM_params.vx_ori = vx_locs(:,4:6);
                    h.cfg.ARM_params.sig_amp_perc = randi(str2num(h.edit_ARM_sig_perc_range.String),1,num_sources); % place holder for now
                    h.cfg.ARM_params.prepost_amp_perc = randi(str2num(h.edit_ARM_prepost_perc_range.String),1,num_sources); % place holder for now
                    h.cfg.ARM_params.sig_latency = repmat(str2num(h.edit_ARM_sig_latency.String),num_sources,1);
                    h.cfg.ARM_params.sig_risetime = repmat(str2num(h.edit_ARM_sig_risetime.String),num_sources,1);
                elseif size(h.cfg.ARM_params.vx_locs,2) == 11   % file was saved by SimMEEG
                    [x,y,z] = sph2cart(deg2rad(vx_locs(:,4)),deg2rad(vx_locs(:,5)),1);
                    h.cfg.ARM_params.vx_ori = [x y z];
                    h.cfg.ARM_params.vx_amp = vx_locs(:,6)';
                    h.cfg.ARM_params.sig_amp_perc = vx_locs(:,7)'; % place holder for now
                    h.cfg.ARM_params.prepost_amp_perc = vx_locs(:,8)'; % place holder for now
                    h.cfg.ARM_params.sig_latency = vx_locs(:,9:10)';
                    h.cfg.ARM_params.sig_risetime = vx_locs(:,11)';
                else
                    warndlg(sprintf('Source Locations, Orientations, & Amplitudes were NOT loaded correctly\nPlease check that the file has data columns in one of these formats:\n[x, y, z]\n[x, y, z, amp]\n[x, y, z, xori, yori, zori, amp]\n'))
                    h.cfg.ARM_params.vx_amp = [];
                end
            end
            h.edit_ARM_num_sources.String = num2str(length(h.cfg.ARM_params.vx_idx));
            h.cfg.ARM_params.sig_amp_perc = randi(str2num(h.edit_ARM_sig_perc_range.String),1,str2num(h.edit_ARM_num_sources.String)); % place holder for now
            h.edit_ARM_num_sources.Enable = 'off';
        end
        sm_ARM_change_source_params('Update');
    case 'Add Waveforms'
        num_sources = length(h.cfg.source.vx_idx); 
        h.cfg.ARM_params.vx_idx = h.cfg.source.vx_idx;  % setting locs to original Synthetic Source locations
        h.cfg.ARM_params.vx_locs = h.cfg.source.vx_locs;
        h.cfg.ARM_params.vx_ori = h.cfg.source.vx_ori;    
        %% but still randomizing amplitudes, sig_perc, and prepost_perc
        h.cfg.ARM_params.vx_amp = randi(str2num(h.edit_ARM_amp_range.String),1,num_sources);
        h.cfg.ARM_params.sig_amp_perc = randi(str2num(h.edit_ARM_sig_perc_range.String),1,num_sources); % place holder for now
        h.cfg.ARM_params.prepost_amp_perc = randi(str2num(h.edit_ARM_prepost_perc_range.String),1,num_sources); % place holder for now
        h.cfg.ARM_params.sig_latency = repmat(str2num(h.edit_ARM_sig_latency.String),num_sources,1);
        h.cfg.ARM_params.sig_risetime = repmat(str2num(h.edit_ARM_sig_risetime.String),num_sources,1);
        h.cfg.ARM_params.source_locs_file = '';
        h.edit_ARM_num_sources.Enable = 'on';
        sm_ARM_change_source_params('Update');
  
end


%% Updating ARM params
h.cfg.ARM_params.combine_type = h.menu_ARM_add.String{h.menu_ARM_add.Value};
h.cfg.ARM_params.num_sources = str2num(h.edit_ARM_num_sources.String);
h.cfg.ARM_params.num_samps = length(h.cfg.study.lat_sim);
h.cfg.ARM_params.num_trials = h.cfg.study.num_trials;
h.cfg.ARM_params.ARMorder = str2num(h.edit_ARM_order.String);
h.cfg.ARM_params.num_interactions = str2num(h.edit_ARM_num_interactions.String);
h.cfg.ARM_params.ori_constraint = h.menu_ARM_ori_constraint.String{h.menu_ARM_ori_constraint.Value};
h.cfg.ARM_params.source_locations_type = h.menu_ARM_source_locs.String{h.menu_ARM_source_locs.Value};



% plot source locations on 3D MRI
h.fcn_handle.plot_3D_mri(); 

