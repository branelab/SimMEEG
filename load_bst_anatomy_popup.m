function load_bst_anatomy_popup(varargin)
%% Francois to do -->
% 1) replace following command with brainstorm_db set by user in brainstorm GUI:  " bst.brainstorm_db =  uigetdir(h.bst_subj_anat_dir,'Set Brainstorm''s database directory'); "
% 2) set "bst.study_menu.Value" to default to the study directory set by user in brainstorm GUI

global h
global bst
h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Loading Brainstorm Anatomy'); drawnow;

if ~isappdata(0, 'BrainstormRunning')   % Brainstorm not running
    w=warndlg(sprintf('\n\nPlease run Brainstorm first and\nthen call SimMEEG within Brainstorm.\n'),'Brainstorm is NOT Running!');
    w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = 11; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
    return
else
    
    %% Initializaing variables
    bst.converted_flag = 0; % 0= anatomy not converted yet; % anatomy converted to FieldTrip format and stored in h.anatomy
    bst.subj_anat_names = ' ';
    bst.subj_data_names = ' ';
    
    bst.mri_files = ' ';
    bst.scalp_files = ' ';
    bst.skull_files = ' ';
    bst.brain_hull_files = ' ';
    bst.brain_cortex_files = ' ';
    bst.brain_wm_files = ' ';
    bst.brain_pial_files = ' ';
    
    bst.sens_meg_files = ' ';
    bst.sens_eeg_files = ' ';
    bst.hdm_eeg_vol_files = ' ';
    bst.hdm_eeg_cortex_files = ' ';
    bst.hdm_meg_vol_files = ' ';
    bst.hdm_meg_cortex_files = ' ';
    
    %% Get brainstorm_db directory
    bst.brainstorm_db =  uigetdir(h.brainstorm_db,'Set Brainstorm''s database directory'); h.brainstorm_db = bst.brainstorm_db;
    fnames = cellstr(ls(bst.brainstorm_db)); bst.study_names = fnames(3:end);
    
    figure(bst.popup_fig);
    
    %% uicontrol default position/sizes
    bst.txt_pos = [.01 .95 .3 .035];
    bst.menu_pos = [.01 .95 .5 .035];
    
    bst.txt_pos2 = [.01 .91 .3 .08];
    bst.menu_pos2 = [.01 .91 .6 .08]; txt_spacing = 1.75;
    
    %% Menu: Study
    bst.study_menu_txt = uicontrol(bst.popup_fig,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('Select Study:'));
    bst.study_menu = uicontrol(bst.popup_fig,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(bst.study_menu_txt.Position([1 3]))+.01 bst.study_menu_txt.Position(2) bst.menu_pos(3) bst.menu_pos(4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.study_names,'Callback',@sm_select_bst_study);
    %% Menu: Subject Anatomy Directory
    bst.subj_anat_menu_txt = uicontrol(bst.popup_fig,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos(1) bst.study_menu_txt.Position(2)-.04 bst.txt_pos(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('Select Subject Anatomy:'));
    bst.subj_anat_menu = uicontrol(bst.popup_fig,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.subj_anat_menu_txt.Position([1 3]))+.01 bst.subj_anat_menu_txt.Position(2) bst.menu_pos(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.subj_anat_names,'Callback',@sm_select_bst_subj_anat);
    %% Menu: Subject Data Directory
    bst.subj_data_menu_txt = uicontrol(bst.popup_fig,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos(1) bst.subj_anat_menu_txt.Position(2)-.04 bst.txt_pos(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('Select Subject Data:'));
    bst.subj_data_menu = uicontrol(bst.popup_fig,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.subj_data_menu_txt.Position([1 3]))+.01 bst.subj_data_menu_txt.Position(2) bst.menu_pos(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.subj_data_names,'Callback',@sm_select_bst_subj_data);
    
    %% %%%%% Panel: Anatomy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delete(bst.anatomy_panel)
    bst.anatomy_panel = uipanel(bst.popup_fig,'Title','Anatomy Files','FontSize',12,'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.01 bst.subj_data_menu.Position(2)-.365 .98 .35],'FontWeight','bold');
    %% Menu: MRI files
    bst.anat_mri_file_txt = uicontrol(bst.anatomy_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('MRI file:'));
    bst.anat_mri_file = uicontrol(bst.anatomy_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_mri_file_txt.Position([1 3]))+.01 bst.anat_mri_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.mri_files);
    %% Menu: Scalp files
    bst.anat_scalp_file_txt = uicontrol(bst.anatomy_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2(1) bst.anat_mri_file_txt.Position(2)-(bst.txt_pos2(4)*txt_spacing) bst.txt_pos2(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('Scalp file:'));
    bst.anat_scalp_file = uicontrol(bst.anatomy_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_scalp_file_txt.Position([1 3]))+.01 bst.anat_scalp_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.scalp_files);
    %% Menu: Skull files
    bst.anat_skull_file_txt = uicontrol(bst.anatomy_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2(1) bst.anat_scalp_file_txt.Position(2)-(bst.txt_pos2(4)*txt_spacing) bst.txt_pos2(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('Skull file:'));
    bst.anat_skull_file = uicontrol(bst.anatomy_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_skull_file_txt.Position([1 3]))+.01 bst.anat_skull_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.skull_files);
    %% Menu: Brain Hull files
    bst.anat_brain_hull_file_txt = uicontrol(bst.anatomy_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2(1) bst.anat_skull_file_txt.Position(2)-(bst.txt_pos2(4)*txt_spacing) bst.txt_pos2(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('Brain Hull file:'));
    bst.anat_brain_hull_file = uicontrol(bst.anatomy_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_brain_hull_file_txt.Position([1 3]))+.01 bst.anat_brain_hull_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.brain_hull_files);
    %% Menu: Brain Cortex files
    bst.anat_brain_cortex_file_txt = uicontrol(bst.anatomy_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2(1) bst.anat_brain_hull_file_txt.Position(2)-(bst.txt_pos2(4)*txt_spacing) bst.txt_pos2(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('Brain Cortex file:'));
    bst.anat_brain_cortex_file = uicontrol(bst.anatomy_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_brain_cortex_file_txt.Position([1 3]))+.01 bst.anat_brain_cortex_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.brain_cortex_files);
    %% Menu: Brain Cortex files
    bst.anat_brain_wm_file_txt = uicontrol(bst.anatomy_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2(1) bst.anat_brain_cortex_file_txt.Position(2)-(bst.txt_pos2(4)*txt_spacing) bst.txt_pos2(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('White Matter file:'));
    bst.anat_brain_wm_file = uicontrol(bst.anatomy_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_brain_wm_file_txt.Position([1 3]))+.01 bst.anat_brain_wm_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.brain_wm_files);
    %% Menu: Pial Cortex files
    bst.anat_brain_pial_file_txt = uicontrol(bst.anatomy_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2(1) bst.anat_brain_wm_file_txt.Position(2)-(bst.txt_pos2(4)*txt_spacing) bst.txt_pos2(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('Brain Pial file:'));
    bst.anat_brain_pial_file = uicontrol(bst.anatomy_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_brain_pial_file_txt.Position([1 3]))+.01 bst.anat_brain_pial_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.brain_pial_files);
    
    
    
    
    sm_select_bst_study();
    
    %% %%%%% Panel: HeadModels & LeadFields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delete(bst.headmodel_panel)
    bst.headmodel_panel = uipanel(bst.popup_fig,'Title','Head Models & Lead Fields','FontSize',12,'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.01 bst.anatomy_panel.Position(2)-.365 .98 .35],'FontWeight','bold');
    %% Menu: MEG Sensor files
    bst.anat_sens_meg_file_txt = uicontrol(bst.headmodel_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('MEG Sensor file:'));
    bst.anat_sens_meg_file = uicontrol(bst.headmodel_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_sens_meg_file_txt.Position([1 3]))+.01 bst.anat_sens_meg_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.sens_meg_files);
    %% Menu: MEG HeadModel Volume files
    bst.anat_hdm_meg_vol_file_txt = uicontrol(bst.headmodel_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2(1) bst.anat_sens_meg_file_txt.Position(2)-(bst.txt_pos2(4)*txt_spacing) bst.txt_pos2(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('HeadModel MEG Volume file:'));
    bst.anat_hdm_meg_vol_file = uicontrol(bst.headmodel_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_hdm_meg_vol_file_txt.Position([1 3]))+.01 bst.anat_hdm_meg_vol_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.hdm_meg_vol_files);
    %% Menu: MEG HeadModel Cortex files
    bst.anat_hdm_meg_cortex_file_txt = uicontrol(bst.headmodel_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2(1) bst.anat_hdm_meg_vol_file_txt.Position(2)-(bst.txt_pos2(4)*txt_spacing) bst.txt_pos2(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('HeadModel MEG Cortex file:'));
    bst.anat_hdm_meg_cortex_file = uicontrol(bst.headmodel_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_hdm_meg_cortex_file_txt.Position([1 3]))+.01 bst.anat_hdm_meg_cortex_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.hdm_meg_cortex_files);
    %% Menu: EEG Sensor files
    bst.anat_sens_eeg_file_txt = uicontrol(bst.headmodel_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2(1) bst.anat_hdm_meg_cortex_file_txt.Position(2)-(bst.txt_pos2(4)*txt_spacing*2) bst.txt_pos2(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('EEG Sensor file:'));
    bst.anat_sens_eeg_file = uicontrol(bst.headmodel_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_sens_eeg_file_txt.Position([1 3]))+.01 bst.anat_sens_eeg_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.sens_eeg_files);
    %% Menu: EEG HeadModel Volume files
    bst.anat_hdm_eeg_vol_file_txt = uicontrol(bst.headmodel_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2(1) bst.anat_sens_eeg_file_txt.Position(2)-(bst.txt_pos2(4)*txt_spacing) bst.txt_pos2(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('HeadModel EEG Volume file:'));
    bst.anat_hdm_eeg_vol_file = uicontrol(bst.headmodel_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_hdm_eeg_vol_file_txt.Position([1 3]))+.01 bst.anat_hdm_eeg_vol_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.hdm_eeg_vol_files);
    %% Menu: EEG HeadModel Cortex files
    bst.anat_hdm_eeg_cortex_file_txt = uicontrol(bst.headmodel_panel,'Style','text','Foregroundcolor','k','Units','normalize',...
        'Position',[bst.txt_pos2(1) bst.anat_hdm_eeg_vol_file_txt.Position(2)-(bst.txt_pos2(4)*txt_spacing) bst.txt_pos2(3:4)],'FontSize',10,'HorizontalAlignment','right',...
        'String',sprintf('HeadModel EEG Cortex file:'));
    bst.anat_hdm_eeg_cortex_file = uicontrol(bst.headmodel_panel,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize','Visible','on',...
        'Position',[sum(bst.anat_hdm_eeg_cortex_file_txt.Position([1 3]))+.01 bst.anat_hdm_eeg_cortex_file_txt.Position(2) bst.menu_pos2(3:4)],...
        'FontSize',8,'HorizontalAlignment','center','String',bst.hdm_eeg_cortex_files);
    
    %% %%%%% Convert Anatomy from BST to FieldTrip %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Button: Convert Anatomy
    bst.btn_convert_bst2ft = uicontrol(bst.popup_fig,'BackgroundColor',[.8 .8 1],'ForegroundColor',[0 0 0],'Style','pushbutton','Units','normalize',...
        'Position',[.3 .075 .3 .04],'Visible','on','UserData',1,...
        'FontSize',10,'HorizontalAlignment','center','String','Convert & Save','Callback',@sm_convert_bst2ft);
    %% Button: Close
    bst.btn_close_bst_anat = uicontrol(bst.popup_fig,'BackgroundColor',[1 .8 .8],'ForegroundColor',[0 0 0],'Style','pushbutton','Units','normalize',...
        'Position',[1-.1 bst.study_menu_txt.Position(2) .075 .03],'Visible','on','UserData',1,...
        'FontSize',10,'HorizontalAlignment','center','String','Close','Callback',@close_bst_anatomy);
    
    sm_select_bst_subj_anat();
    sm_select_bst_subj_data();
    
    waitfor(bst.popup_fig);
end
%%
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');


function sm_select_bst_study(varargin)
global bst

bst.subj_data_menu.Value = 1;
bst.subj_anat_menu.Value = 1;

bst.study_dir = fullfile(bst.brainstorm_db,bst.study_names{bst.study_menu.Value});
bst.study_anat_dir = fullfile(bst.study_dir,'anat');
bst.study_data_dir = fullfile(bst.study_dir,'data');
fnames = cellstr(ls(bst.study_anat_dir)); bst.subj_anat_names = fnames(3:end);
fnames = cellstr(ls(bst.study_data_dir)); bst.subj_data_names = fnames(3:end);
bst.subj_anat_menu.String = bst.subj_anat_names;
bst.subj_data_menu.String = bst.subj_data_names;
sm_select_bst_subj_anat(); sm_select_bst_subj_data();

function sm_select_bst_subj_anat(varargin)
global bst

%% Anatomy Files: subject files within anatomy directory
bst.subj_anat_dir = fullfile(bst.study_anat_dir,bst.subj_anat_menu.String{bst.subj_anat_menu.Value});
fnames = cellstr(ls(bst.subj_anat_dir)); fnames= fnames(3:end); fnames = [{'none'}; fnames];
% populate all antomy file menus with all names in directory
bst.anat_mri_file.String = fnames;
bst.anat_scalp_file.String = fnames;
bst.anat_skull_file.String = fnames;
bst.anat_brain_hull_file.String = fnames;
bst.anat_brain_cortex_file.String = fnames;
bst.anat_brain_wm_file.String = fnames;
bst.anat_brain_pial_file.String = fnames;

%% MRI file --> set mri file to likely default of 'subjectimage'
xidx = find(startsWith(fnames,'subjectimage')==1);
if isempty(xidx); bst.anat_mri_file.Value = 1; else; bst.anat_mri_file.Value = xidx; end
%% Scalp file
xidx = find(startsWith(fnames,'tess_head')==1);
if isempty(xidx); bst.anat_scalp_file.Value = 1; else; bst.anat_scalp_file.Value = xidx; end
%% Skull file
xidx = find(startsWith(fnames,'tess_outerskull')==1);
if isempty(xidx); bst.anat_skull_file.Value = 1; else; bst.anat_skull_file.Value = xidx; end
%% Brain Hull file
xidx = find(startsWith(fnames,'tess_innerskull')==1);
if isempty(xidx); bst.anat_brain_hull_file.Value = 1; else; bst.anat_brain_hull_file.Value = xidx; end
%% Brain Cortex file
xidx = find(startsWith(fnames,'tess_cortex_pial_low')==1);
if isempty(xidx); bst.anat_brain_cortex_file.Value = 1; else; bst.anat_brain_cortex_file.Value = xidx; end
%% Brain White Matter file
xidx = find(startsWith(fnames,'tess_cortex_white_low')==1);
if isempty(xidx); bst.anat_brain_wm_file.Value = 1; else; bst.anat_brain_wm_file.Value = xidx; end
%% Brain Pial file
xidx = find(startsWith(fnames,'tess_cortex_pial_cereb_low')==1);
if isempty(xidx); bst.anat_brain_pial_file.Value = 1; else; bst.anat_brain_pial_file.Value = xidx; end
function sm_select_bst_subj_data(varargin)
global bst

%% Head Models & Lead Fields: subject files within subject data directory
bst.subj_data_dir = fullfile(bst.study_data_dir,bst.subj_data_menu.String{bst.subj_data_menu.Value});
% fnames = cellstr(ls(bst.subj_data_dir)); fnames= fnames(3:end); fnames = [{'none'}; fnames];
% if sum(ismember(fnames,'@default_study'))>0
%     bst.subj_data_dir = fullfile(bst.subj_data_dir,'@default_study');
%     fnames = cellstr(ls(bst.subj_data_dir)); fnames= fnames(3:end); fnames = [{'none'}; fnames];
% else
% end
[fpath,fname,~] = fileparts(bst.subj_data_dir);
if ~strcmpi(fname,'@default_study')
    bst.subj_data_dir = fullfile(fpath,fname);
    bst.subj_data_dir = uigetdir(bst.subj_data_dir,'Set Subject''s Data Directory to get Channels & HeadModels');
    fnames = cellstr(ls(bst.subj_data_dir)); fnames= fnames(3:end); fnames = [{'none'}; fnames];
else
    fnames = cellstr(ls(bst.subj_data_dir)); fnames= fnames(3:end); fnames = [{'none'}; fnames];
end


% populate all menus with all names in directory
bst.anat_sens_meg_file.String = fnames;
bst.anat_hdm_meg_vol_file.String = fnames;
bst.anat_hdm_meg_cortex_file.String = fnames;
bst.anat_sens_eeg_file.String = fnames;
bst.anat_hdm_eeg_vol_file.String = fnames;
bst.anat_hdm_eeg_cortex_file.String = fnames;


%% Sens MEG file
xidx = find(startsWith(fnames,'channel')==1);
if isempty(xidx); bst.anat_sens_meg_file.Value = 1; else; bst.anat_sens_meg_file.Value = xidx; end
%% HeadModel & LeadField Volume MEG file
yidx = regexpi(fnames,regexptranslate('wildcard','headmodel_vol*meg*.mat'));
xidx = false(size(yidx)); for t=1:length(yidx); if yidx{t}==1; xidx(t)=true; end; end
xidx = find(xidx==1);
if isempty(xidx); bst.anat_hdm_meg_vol_file.Value = 1; else; bst.anat_hdm_meg_vol_file.Value = xidx(1); end
%% HeadModel & LeadField Cortex MEG file
yidx = regexpi(fnames,regexptranslate('wildcard','headmodel_surf*meg*.mat'));
xidx = false(size(yidx)); for t=1:length(yidx); if yidx{t}==1; xidx(t)=true; end; end
xidx = find(xidx==1);
if isempty(xidx); bst.anat_hdm_meg_cortex_file.Value = 1; else; bst.anat_hdm_meg_cortex_file.Value = xidx(1); end

%% Sens EEG file
xidx = find(startsWith(fnames,'channel')==1);
if isempty(xidx); bst.anat_sens_eeg_file.Value = 1; else; bst.anat_sens_eeg_file.Value = xidx; end
%% HeadModel & LeadField Volume MEG file
yidx = regexpi(fnames,regexptranslate('wildcard','headmodel_vol*eeg*.mat'));
xidx = false(size(yidx)); for t=1:length(yidx); if yidx{t}==1; xidx(t)=true; end; end
xidx = find(xidx==1);
if isempty(xidx); bst.anat_hdm_eeg_vol_file.Value = 1; else; bst.anat_hdm_eeg_vol_file.Value = xidx(1); end
%% HeadModel & LeadField Cortex MEG file
yidx = regexpi(fnames,regexptranslate('wildcard','headmodel_surf*eeg*.mat'));
xidx = false(size(yidx)); for t=1:length(yidx); if yidx{t}==1; xidx(t)=true; end; end
xidx = find(xidx==1);
if isempty(xidx); bst.anat_hdm_eeg_cortex_file.Value = 1; else; bst.anat_hdm_eeg_cortex_file.Value = xidx(1); end

function close_bst_anatomy(varargin)
global h; global bst;

answ = questdlg('Are you sure you want to close Brainstorm Anatomy popup?','Close?','Yes','No','No');

switch answ
    case 'Yes'
        if isempty(h.anatomy) || bst.converted_flag==0  % not converted
            answ2 = questdlg(sprintf('Anatomy data has not been converted?\n\nDo you want to close without converting?\n'),'Close?','Yes','No','No');
            switch answ2
                case 'Yes'
                    close(bst.popup_fig)
                case 'No'
            end
        else
            close(bst.popup_fig)
        end
    case 'No'
end





