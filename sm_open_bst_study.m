function [h] = sm_open_bst_study(h)
% This function will open a BST study protocol with ProtocolName, if the ProtocolName doesn't exist it will be created.
% If the "SimMEEG_default_zip" filename is given then SimMEEG's default Anatomy will be uploaded and used given that the use_SimMEEG_default_flag=1.

%% Prepare test folders
% Start Brainstorm with GUI
if ~isappdata(0, 'BrainstormRunning')
%     bst_startup(h.Brainstorm_dir, h.bst.gui_flag, h.bst.BrainstormDbDir); % running as nogui because it was called from SimMEEG directly
    bst_startup(h.Brainstorm_dir, h.bst.gui_flag); % 
    h.bst.BrainstormDbDir = bst_get('BrainstormDbDir');
else
    h.bst.BrainstormDbDir = bst_get('BrainstormDbDir');
end


%% ===== START =====
%% Check inputs
if isempty(h.bst.BrainstormDbDir)
    error('You must specify an empty directory in input.');
end
% Create test folder
if ~isfolder(h.bst.BrainstormDbDir)
    if ~mkdir(h.bst.BrainstormDbDir)
        error(['Could not create folder: ' h.bst.BrainstormDbDir]);
    end
end

%% default SimMEEG Antomy/Data for Brainstorms Database
h.bst.use_SimMEEG_default_flag =1; % default for now is to only use SimMEEG Anatomy Protocol in BST
h.bst.SimMEEG_default_zip = 'C:\BRANELab\Software\SimMEEG\anatomy\SimMEEG_Anatomy.zip';
[fpath,fname,fect] = fileparts(h.bst.SimMEEG_default_zip);
if exist(h.bst.SimMEEG_default_zip,'file')
else
    [fname,fpath] = uigetfile('*.zip','Select Brainstorm Anatomy.zip file');
    fext = '.zip';
    h.bst.SimMEEG_default_zip = fullfile(fpath,fname);
end
h.bst.ProtocolName = fname;

if isappdata(0, 'BrainstormRunning')
    h.bst.BrainstormDbDir = bst_get('BrainstormDbDir');
else
    h.bst.BrainstormDbDir = fullfile(h.simmeeg_dir,'SimMEEG_temp_db');
end
bst_set('BrainstormDbDir', h.bst.BrainstormDbDir);

%% Unzip default Anatomy if doesn't exist
[fpath,h.bst.ProtocolName,fext]= fileparts(h.bst.SimMEEG_default_zip);
h.bst.iProtocol = bst_get('iProtocol',h.bst.ProtocolName);
if h.bst.use_SimMEEG_default_flag==1 && isempty(bst_get('iProtocol',h.bst.ProtocolName)) % flag to unzip SimMEEG default Anatomy/Study Protocol in Brainstorm
    %     SimMEEG_default_zip = 'C:\BRANELab\Software\SimMEEG\anatomy\SimMEEG_BST_Anatomy.zip';
    % study_dir = fullfile(BrainstormDbDir,ProtocolName);
    % fprintf('Unzipping Default "anat" and "data" from %s \ninto Brainstorm Temp database: %s\n',SimMEEG_default_zip,study_dir);
    % unzip(SimMEEG_default_zip,study_dir);
    import_protocol(h.bst.SimMEEG_default_zip)
    
elseif h.bst.use_SimMEEG_default_flag==2 % Convert current SimMEEG anat into Brainstorm format and save into brainstorm temp dir
    % Not implemented yet still need to write code for this
elseif h.bst.use_SimMEEG_default_flag==0  % Add temp Protocol study and
    % Delete existing protocol
    gui_brainstorm('DeleteProtocol', h.bst.ProtocolName);
    % Create new protocol
    gui_brainstorm('CreateProtocol', h.bst.ProtocolName, 0, 0, h.bst.BrainstormDbDir);
else
    fprintf('Brainstorm Protocol Study already Exists: %s\nIf you want to start a new study, please delete this Protocol/folder.\n\n',fullfile(h.bst.BrainstormDbDir,h.bst.ProtocolName)); 
end

h.bst.ProtocolInfo = bst_get('ProtocolInfo');
h.bst.iProtocol = bst_get('iProtocol',h.bst.ProtocolName);

%% find study with same sensors and headmodels - set same headmodel as current inverse solution type
% h.bst.iStudy = h.bst.ProtocolInfo.iStudy;
% h.bst.Study = bst_get('Study',h.bst.iStudy);
[h] = sm_bst_update_study(h);
