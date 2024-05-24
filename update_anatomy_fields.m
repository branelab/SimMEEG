function update_anatomy_fields(varargin)
global h

h.anatomy_file_txt.String = h.anat_file;
h.btn_select_source_locs.Enable = 'on';

try
    h.mri_txt.String = sprintf('MRI Size: %.f x %.f x %.f Unit(%s) Coord(%s)',size(h.anatomy.mri.anatomy),h.anatomy.mri.unit,h.anatomy.mri.coordsys);
catch
    try
    h.mri_txt.String = sprintf('MRI Size: %.f x %.f x %.f Unit(%s) Coord(%s)',h.anatomy.mri.dim,h.anatomy.mri.unit,h.anatomy.mri.coordsys);
    catch
        h.mri_txt.String  = 'No MRI ';
        warndlg(sprintf('\nNo MRI Exist!\n\nPlease load an MRI'));
    end
end

try
    h.headmodel_txt.String = sprintf('Head model Type: %s',h.anatomy.headmodel.type); % starting with EEG
catch
    h.headmodel_txt.String = 'No Head Models';
end
try
    h.leadfield_txt.String = sprintf('%.f sensors x  %.f orientations x %.f dipoles',size(h.anatomy.leadfield.H));
    h.menu_head_model.Enable = 'on';
catch
    h.leadfield_txt.String = 'No Lead Fields'; h.menu_head_model.Enable = 'inactive';
end

try
    h.mri_path.String = h.anat_path;
catch
     h.mri_path.String = ' ';
end

%% Sensor String
if ~isempty(h.anatomy.leadfield) && ~isempty(h.anatomy.sens)
    h.sensors_txt.String = sprintf('%.f sensors (%s) with %.f used for leadfield',size(h.anatomy.sens.chanpos,1),h.anatomy.sens.type,size(h.anatomy.leadfield.H,1));
    h.btn_run_sim_meeg.Enable = 'on'; h.btn_run_sim_meeg.ForegroundColor = 'k';
    h.btn_sim_sens_noise.Enable = 'on'; h.btn_sim_sens_noise.ForegroundColor = 'k';
elseif isempty(h.anatomy.leadfield) && ~isempty(h.anatomy.sens)
    h.sensors_txt.String = sprintf('%.f sensors (%s) -- No Lead Field',size(h.anatomy.sens.chanpos,1),h.anatomy.sens.type);
    h.menu_sens_type.Enable = 'on';
    h.btn_run_sim_meeg.Enable = 'inactive'; h.btn_run_sim_meeg.ForegroundColor = h.btn_run_sim_meeg.BackgroundColor*.7;
    h.btn_sim_sens_noise.Enable = 'inactive'; h.btn_sim_sens_noise.ForegroundColor = h.btn_sim_sens_noise.BackgroundColor*.7;
    
else
    h.sensors_txt.String = 'No Sensors'; h.menu_sens_type.Enable = 'inactive';
end

%% orientation constraint
% try
if ~isfield(h.anatomy.leadfield,'ori')    % cortical surface leadfields exist
    h.menu_head_model.Enable = 'on'; %h.menu_inv_headmodel.Enable = 'on';
    h.menu_head_model.BackgroundColor = [1 1 1];  h.menu_head_model.ForegroundColor = [1 1 1]*0;
    h.pwd = pwd;
    cd(fullfile(h.FieldTrip_dir,'private\'));
%     for n=1:length(h.anatomy.leadfield_eeg_cortex); h.anatomy.leadfield_eeg_cortex(n).ori = normals(h.anatomy.mesh_cortex.pos,h.anatomy.mesh_cortex.tri,'vertex'); end% finding orientations normalized to cortical surface "cortically contrained"
       h.anatomy.leadfield.ori = normals(h.anatomy.mesh_cortex.pos,h.anatomy.mesh_cortex.tri,'vertex'); % finding orientations normalized to cortical surface "cortically contrained"
    cd(h.pwd);
% else
%     h.menu_head_model.Enable = 'inactive'; %h.menu_inv_headmodel.Enable = 'inactive';
%     h.menu_head_model.Value = 1;
%     h.menu_head_model.BackgroundColor = [1 1 1]*1;  h.menu_head_model.ForegroundColor = [1 1 1]*0;
end
% catch
% end

%% Sensor Menu
if ~isempty(h.anatomy.sens_eeg) && ~isempty(h.anatomy.sens_meg)
    h.menu_sens_type.Enable = 'on';
else
    if ~isempty(h.anatomy.sens_meg); h.menu_sens_type.Value = 1; 
    elseif ~isempty(h.anatomy.sens_eeg); h.menu_sens_type.Value = 2;
    else
        warndlg(sprintf('\nNo Sensors Exist!\n\nPlease load MEG or EEG sensors'));
    end
    h.menu_sens_type.Enable = 'inactive';
end

try

    %% HeadModel Menu
if isempty(h.anatomy.headmodel_eeg_vol) && isempty(h.anatomy.headmodel_eeg_cortex) &&  ...
        isempty(h.anatomy.headmodel_meg_vol) && isempty(h.anatomy.headmodel_meg_cortex)
    warndlg( sprintf('\nNo Head Models Exist\n\nPlease load or create Head Models\n') );
    h.menu_head_model.Enable = 'inactive';  % No Headmodels found
    h.menu_inv_headmodel.Enable = 'inactive';  % No Headmodels found
elseif  ( ~isempty(h.anatomy.headmodel_eeg_vol) || ~isempty(h.anatomy.headmodel_meg_vol) ) && ... % Volume only
        ( isempty(h.anatomy.headmodel_eeg_cortex) && isempty(h.anatomy.headmodel_meg_cortex) )
    h.menu_head_model.Enable = 'inactive'; h.menu_head_model.Value = 1;
    h.menu_inv_headmodel.Enable = 'inactive'; h.menu_inv_headmodel.Value = 1; % Volume only
    h.menu_ori_normal.Enable = 'inactive'; h.menu_ori_normal.Value=1;  h.menu_ori_normal.ForegroundColor=[1 1 1]*.5;
elseif  ( isempty(h.anatomy.headmodel_eeg_vol) && isempty(h.anatomy.headmodel_meg_vol) ) && ... % Cortex only
        ( ~isempty(h.anatomy.headmodel_eeg_cortex) || ~isempty(h.anatomy.headmodel_meg_cortex) )
    h.menu_head_model.Enable = 'inactive'; h.menu_head_model.Value = 2; % Cortex only
    h.menu_inv_headmodel.Enable = 'inactive'; h.menu_inv_headmodel.Value = 2; % Cortex only
elseif  ( ~isempty(h.anatomy.headmodel_eeg_vol) || ~isempty(h.anatomy.headmodel_meg_vol) ) && ... % Volume & Cortex exist
        ( ~isempty(h.anatomy.headmodel_eeg_cortex) || ~isempty(h.anatomy.headmodel_meg_cortex) )
    h.menu_head_model.Enable = 'on'; %h.menu_head_model.Value = 1; % Volume & Cortex
    h.menu_inv_headmodel.Enable = 'on'; %h.menu_inv_headmodel.Value = 1; % Volume & Cortex
    h.menu_ori_normal.Enable = 'on';  h.menu_ori_normal.ForegroundColor=[1 1 1]*0;
end
catch
    warning('Headmodel not present in h.anatomy. Please Create HeadModels.');
end

% update channel list box
try
    if h.menu_sens_type.Value==1 % MEG
%         h.listbox_chans.String = h.anatomy.sens_meg.label; h.listbox_chans.ForegroundColor = h.sens_clr;
        h.listbox_chans.String = h.anatomy.sens.label; h.listbox_chans.ForegroundColor = h.sens_clr;
    elseif h.menu_sens_type.Value==2 % EEG
%         h.listbox_chans.String = h.anatomy.sens_eeg.label; h.listbox_chans.ForegroundColor = h.chan_clr;
        h.listbox_chans.String = h.anatomy.sens.label; h.listbox_chans.ForegroundColor = h.chan_clr;
    end
catch
end

%% set GOOD and BAD channels
try
    if isfield(h.anatomy.sens,'bad_sensors')
        if ~isempty(h.anatomy.sens.bad_sensors)
            h.anatomy.sens.good_sensors = setdiff(1:length(h.anatomy.sens.label),h.anatomy.sens.bad_sensors);
        else
            h.anatomy.sens.good_sensors = 1:size(h.anatomy.leadfield.H,1);
        end
    end
    if length(h.anatomy.leadfield.label)~=length(h.anatomy.sens.label)
        
        h.anatomy.sens.bad_sensors = [];
        h.anatomy.sens.good_sensors = 1:size(h.anatomy.sens.chanpos,1); %size(h.anatomy.leadfield.H,1);
        h.anatomy.sens.leadfield_idx = (1:length(h.anatomy.leadfield.label))';
        % find same labels
        sens_idx = zeros(length(h.anatomy.leadfield.label),1); % for MEG there are reference sensors thus needed to identify those over scalp and those used as reference sensors
        for v=1:length(h.anatomy.leadfield.label)
            if ~isempty(find(strcmpi(h.anatomy.sens.label,h.anatomy.leadfield.label(v)))==1)
            sens_idx(v) = 1;
            end
        end
        h.anatomy.sens.sens_idx = find(sens_idx==1); 
    else
        h.anatomy.sens.leadfield_idx = (1:length(h.anatomy.leadfield.label))';
        h.anatomy.sens.sens_idx = (1:length(h.anatomy.leadfield.label))';
    end
catch
end

%% showing #sensor for loaded MEG and EEG
if h.menu_sens_type.Value == 1
        sens_num = [];
        h.menu_sens_montage.Visible = 'on';
        for v=1:length(h.anatomy.sens_meg)
            sens_num(v) = length(h.anatomy.sens_meg(v).label);
        end
        h.menu_sens_montage.String = num2str(sens_num'); 
elseif h.menu_sens_type.Value == 2
        sens_num = [];
        h.menu_sens_montage.Visible = 'on';
        for v=1:length(h.anatomy.sens_eeg)
            sens_num(v) = length(h.anatomy.sens_eeg(v).label);
        end
        h.menu_sens_montage.String = num2str(sens_num'); 
end
        
%% update monte sensor montage #
if ~isempty(h.anatomy.sens_meg)
    sens_num = [];
    for v=1:length(h.anatomy.sens_meg)
        sens_num(v) = length(h.anatomy.sens_meg(v).label);
    end
    h.listbox_monte_MEG_sens_montage.String = num2str(sens_num');
end
if ~isempty(h.anatomy.sens_eeg)
    sens_num = [];
    for v=1:length(h.anatomy.sens_eeg)
        sens_num(v) = length(h.anatomy.sens_eeg(v).label);
    end
    h.listbox_monte_EEG_sens_montage.String = num2str(sens_num');
end


%% Update Headmodel Menu
hdm_name = {'vol' 'cortex' 'spheres_cortex' 'spheres_vol'};
hdm_flag = false(length(hdm_name),1);
for a=1:length(hdm_name)
    hdm = sprintf('headmodel_%s_%s', h.anatomy.sens.type, hdm_name{a});
    fn = fieldnames(h.anatomy); fn = fn(contains(fn,hdm));
    for f=1:length(fn)
        if ~isempty(h.anatomy.(fn{f}))
            hdm_flag(a) = true;
        end
    end
end
hdm_string = {"Volume" "Cortical Surface" "Spheres(Volume)" "Spheres(Cortical)"}; 
h.menu_head_model.String = hdm_string(hdm_flag);




