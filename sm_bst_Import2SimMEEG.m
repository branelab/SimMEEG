function sm_bst_Import2SimMEEG(varargin)

global h

%% get Studies in Protocol
h.bst.Protocol = bst_get('ProtocolInfo');
h.bst.Subjects = bst_get('ProtocolSubjects');

%% delete old Figure if present
if isfield(h.bst,'fig_bst2SimMEEG')
    if isvalid(h.bst.fig_bst2SimMEEG.main_fig)
        close(h.bst.fig_bst2SimMEEG.main_fig); h.bst = rmfield(h.bst,'fig_bst2SimMEEG');
    end
end
% if isfield(h.bst.fig_bst2SimMEEG,'Study')
%     h.bst = rmfield(h.bst.fig_bst2SimMEEG,'Study');
% end

%% %%%%% Create new Figure  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h.bst.fig_bst2SimMEEG.main_fig = uifigure; h.bst.fig_bst2SimMEEG.main_fig.Position = [500 150 700 800];
h.bst.fig_bst2SimMEEG.main_fig.MenuBar = 'none';
h.bst.fig_bst2SimMEEG.main_fig.Name = 'Import from Brainstorm to SimMEEG';
fig_pos = h.bst.fig_bst2SimMEEG.main_fig.Position;

main_pos = [fig_pos(3:4) fig_pos(3:4)];    % position for normalizing pixels using
set(h.bst.fig_bst2SimMEEG.main_fig,'WindowKeyPressFcn',@KeyPress);
set(h.bst.fig_bst2SimMEEG.main_fig,'WindowKeyReleaseFcn',@KeyPress);

%% Text: Select Datasets
h.bst.fig_bst2SimMEEG.txt_select_ds = uilabel(h.bst.fig_bst2SimMEEG.main_fig,'Position',[.01 .95 .49 .05].*main_pos, ...
    'Text','Loading Datasets ... This can take time.','FontSize',14, 'FontColor','r');

drawnow;
%% ListBox: Study Tree
% delete(h.bst.fig_bst2SimMEEG.study_tree);
h.bst.fig_bst2SimMEEG.study_tree = uitree(h.bst.fig_bst2SimMEEG.main_fig,'Position',[.01 .01 .49 .94].*main_pos);
h.bst.fig_bst2SimMEEG.study_tree.SelectionChangedFcn = @(src, event)select_all(src, event);
h.bst.fig_bst2SimMEEG.Study.study_node = uitreenode(h.bst.fig_bst2SimMEEG.study_tree,'Text',h.bst.Protocol.Comment,'NodeData',[]);
%% Btn: Import
% delete(h.bst.fig_bst2SimMEEG.btn_import_bst2simmeeg)
h.bst.fig_bst2SimMEEG.btn_import_bst2simmeeg = uibutton(h.bst.fig_bst2SimMEEG.main_fig, ...
    "Text","Import Data", "Position",[.51 .91 .15 .04].*main_pos, "FontColor", [0 .5 0],...
    "BackgroundColor", [.9 1 .9], 'Visible', 'off',...
    "ButtonPushedFcn",@(src,event)sm_bst_import_bst2simmeeg);

%% drawnow;
drawnow;

%% Create tree

if isempty(h.bst.Subjects.Subject)
    h.bst.fig_bst2SimMEEG.Study.subj_node(1) = uitreenode(h.bst.fig_bst2SimMEEG.Study.study_node,'Text', 'No Data in Brainstorm Protocol', 'NodeData', [1]);
else
    %% Subj
    for n=1:length(h.bst.Subjects.Subject)
        SubjectFile = h.bst.Subjects.Subject(n).FileName;
        CondFiles = bst_get('StudyWithSubject',   SubjectFile);
        h.bst.fig_bst2SimMEEG.Study.Subj(n).node = uitreenode(h.bst.fig_bst2SimMEEG.Study.study_node,'Text', h.bst.Subjects.Subject(n).Name, 'NodeData', [n]);
        data.Type = 'Subject'; data.Name = h.bst.Subjects.Subject(n).Name;
        h.bst.fig_bst2SimMEEG.Study.subj(n).node.UserData = data;
        %% Condition
        for c=1:length(CondFiles)
            if ~isempty(CondFiles(c).Condition)
                h.bst.fig_bst2SimMEEG.Study.Subj(n).Cond(c).node = uitreenode(h.bst.fig_bst2SimMEEG.Study.Subj(n).node,'Text', CondFiles(c).Condition{:}, 'NodeData', [n c]);
                data.Type = 'Condition'; data.Name = CondFiles(c).Condition{:};
                h.bst.fig_bst2SimMEEG.Study.Subj(n).Cond(c).node.UserData = data;
                DataFiles = CondFiles(c).Data;
                if ~isempty(DataFiles) && length(DataFiles)>1
                    %% Find uniue datasets - filenames separated by "data_cond_*.mat"
                    dname = {DataFiles.FileName};
                    [~,dsname] = fileparts({DataFiles.FileName});

                    ds = DataFiles;
                    for d=1:length(dsname)  % parsing out condition from filename
                        dc = strsplit(dsname{d},'_');
                        dcidx = find(contains(dc,'trial'), 1);
                        if isempty(dcidx) % not trial data
                            ds(d).condition = dsname{d}(6:end);
                        elseif ~isempty(dcidx) && length(dc)>3 % trial data but pre-processed info in name - essentially another condition
                            ds(d).condition = [dc{[1:dcidx-1 dcidx+1:end]}];
                        else
                            ds(d).condition = dc{2};
                        end
                        sidx = findstr( ds(d).Comment, '(');
                        if ~isempty(sidx)
                            ds(d).dataset_name = ds(d).Comment(1:sidx(1)-2);
                        else
                            ds(d).dataset_name = ds(d).Comment;
                        end
                        ds(d).file_idx = d;

                    end
                    % DataSets under "Condtiions"
                    ds_code = unique({ds.condition});
                    for d=1:length(ds_code)
                        tidx = find(strcmp({ds.condition}, ds_code{d}));
                        sidx1 = findstr( ds(tidx(1)).Comment, '('); sidx2 = findstr( ds(tidx(1)).Comment, ')');
                        if isempty(sidx1)
                            ds_name = sprintf('%s (%.f)',ds(tidx(1)).dataset_name, length(tidx));
                        else
                            dt =  ds(tidx(1)).Comment([1:sidx1-2 sidx2+1:end]);
                            ds_name = sprintf('%s (%.f)', dt, length(tidx));
                        end

                        h.bst.fig_bst2SimMEEG.Study.Subj(n).Cond(c).DataSet(d).node = uitreenode(h.bst.fig_bst2SimMEEG.Study.Subj(n).Cond(c).node,'Text', ds_name, 'NodeData', [n c d]);
                        data.Type = 'DataSet'; data.Name = ds_name;
                        h.bst.fig_bst2SimMEEG.Study.Subj(n).Cond(c).DataSet(d).node.UserData = data;
                        %% Trials
                        for t=1:length(tidx)  % adding
                            ds_txt = ds(tidx(t)).Comment;
                            if  ds(tidx(t)).BadTrial ==1; ds_txt = [ds_txt ' [BAD]']; end % bad trial
                            h.bst.fig_bst2SimMEEG.Study.Subj(n).Cond(c).DataSet(d).Trial(t).node = uitreenode(h.bst.fig_bst2SimMEEG.Study.Subj(n).Cond(c).DataSet(d).node, 'Text', ds_txt, 'NodeData', [n c d t]);
                            data.Type = 'Trial'; data.Name = file_fullpath(ds(tidx(t)).FileName);
                            h.bst.fig_bst2SimMEEG.Study.Subj(n).Cond(c).DataSet(d).Trial(t).node.UserData = data;
                        end
                    end
                elseif ~isempty(DataFiles) && length(DataFiles)==1
                    h.bst.fig_bst2SimMEEG.Study.Subj(n).Cond(c).DataSet(1).Trial(1).node  = uitreenode(h.bst.fig_bst2SimMEEG.Study.Subj(n).Cond(c).node, 'Text', DataFiles.Comment, 'NodeData', [n c 1 1]);
                    data.Type = 'Trial'; data.Name = file_fullpath(DataFiles.FileName);
                    h.bst.fig_bst2SimMEEG.Study.Subj(n).Cond(c).DataSet(1).Trial(1).node.UserData = data;
                else
                end
            end
        end
    end
end

%% Expand nodes
for n=1:length(h.bst.fig_bst2SimMEEG.Study.Subj)
    expand(h.bst.fig_bst2SimMEEG.Study.Subj(n).node);
end
expand(h.bst.fig_bst2SimMEEG.study_tree); drawnow;

%% update figure
h.bst.fig_bst2SimMEEG.study_tree.Multiselect = "on";

h.bst.fig_bst2SimMEEG.txt_select_ds.Text = 'Select Brainstorm Datasets:';
h.bst.fig_bst2SimMEEG.txt_select_ds.FontColor = 'k';
h.bst.fig_bst2SimMEEG.btn_import_bst2simmeeg.Visible = 'on';

end

function KeyPress(Source, EventData)

switch EventData.EventName
    case 'WindowKeyPress'
        fprintf('Pressed: %s\n',EventData.Key)
        Source.UserData = EventData.Key;
    case 'WindowKeyRelease'
        fprintf('Released: %s\n',EventData.Key)
        Source.UserData = '';
end

end

function select_all(varargin)
% select all nodes underneath varargin source
global h

sel_nodes = varargin{2}.SelectedNodes;
prev_nodes = varargin{2}.PreviousSelectedNodes;
sel_node = setxor(sel_nodes, prev_nodes);
if length(sel_node)>1
    sel_node = sel_nodes;
end

% select all trial nodes below
if ~isempty(sel_node.UserData) && strcmp(sel_node.UserData.Type,'DataSet')   % ordinate dataset node
    x = sel_node.NodeData;
    new_nodes = [h.bst.fig_bst2SimMEEG.Study.Subj(x(1)).Cond(x(2)).DataSet(x(3)).Trial(:).node]';
    if strcmp(h.bst.fig_bst2SimMEEG.main_fig.UserData,'control')
        if any(ismember(new_nodes, prev_nodes))     % when ctrl key pressed on already selected nodes - remove nodes from list
            collapse(sel_node);
            prev_nodes = setxor(new_nodes, prev_nodes);
            prev_nodes = setxor(sel_node, prev_nodes);
            h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes = [prev_nodes];
        else
            expand(sel_node);
            h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes = [sel_node; prev_nodes; new_nodes];
        end
    else
        expand(sel_node);
        h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes = [sel_node; new_nodes];
    end
end

end


function sm_bst_import_bst2simmeeg(varargin)
% imports selected files from BST into SimMEEG as "sens_final" and converts to SimMEEG format
global h
if isempty(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes)
    warndlg('Please select Brainstorm datasets', 'No Datasets Selected');
else
    %% Updating Import Figure
    %     h.bst.fig_bst2SimMEEG.study_tree.Enable = 'off';
    %     h.bst.fig_bst2SimMEEG.btn_import_bst2simmeeg.Enable = 'off';
    %     h.bst.fig_bst2SimMEEG.txt_select_ds.Text = 'Importing Brainstorm Datasets into SimMEEG ...';
    %     h.bst.fig_bst2SimMEEG.txt_select_ds.FontColor = 'r';
    %     drawnow;

    %% overwriting h.sim_data
    h.sim_data = struct('sig_final', [], 'sig_wav', [],'prepost_wav', [], 'noise_wav', [], 'cfg', [], 'prepost_win', [], ...
        'sig_win', [], 'source_waveform_type', [], 'intersource_correlations', [], 'sens_noise_type', 'Brainstorm Imported', 'sens_noise', [], 'sens_noise_scaled', [], ...
        'sens_sig_data', [], 'sens_final', [], 'noise_gain', [], 'sens_noise_final', []);

    %% resetting SimMEEG study parameters
    % load data
    for f=1:length(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes)
        ftype{f} = h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes(f).UserData.Type;
        fname{f} = h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes(f).UserData.Name;
    end
    tidx = find(contains(ftype,'Trial'));
    fn = load(file_fullpath(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes(tidx(1)).UserData.Name));

    %% get channel info
    ChannelFile = file_fullpath(bst_get('ChannelFileForStudy', file_fullpath(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes(tidx(1)).UserData.Name)));
    chans = bst_get('ChannelFile',ChannelFile);
    ch = load(ChannelFile);
    if any(contains( chans.Channel.Modalities, 'MEG') )
        chan_idx = find(contains({ch.Channel.Type},'MEG'));
        unit_scale = 1e15;  % convert to micoV or femtoT
    elseif any(contains( chans.Channel.Modalities, 'EEG') )
        chan_idx = find(contains({ch.Channel.Type},'EEG'));
        unit_scale = 1e6;  % convert to micoV or femtoT
    else
        warndlg('Sensors file does not contain "MEG" or "EEG" channels','No MEG/EEG sensors');
    end
    % rescaling units
    % check sensor units V, microV, ... SimMEEG needs to be in femtoT
    sigma_val = mean(std(fn.F(chan_idx,:),[],2));
    if sigma_val>1
        unit_scale = 1;   % convert to micoV or femtoT
    elseif sigma_val>1e-3 && sigma_val<1  % milli
        unit_scale = 1e3;   % convert to micoV or femtoT
    elseif sigma_val>1e-6 && sigma_val<1e-3  % micro
        unit_scale = 1e6;  % convert to micoV or femtoT
    elseif sigma_val>1e-9 && sigma_val<1e-6  % nan
        unit_scale = 1e9;  % convert to micoV or femtoT
    elseif sigma_val>1e-12 && sigma_val<1e-9  % pico
        unit_scale = 1e12;  % convert to micoV or femtoT
    elseif sigma_val>1e-15 && sigma_val<1e-12  % femto
        unit_scale = 1e15;  % convert to micoV or femtoT
    end

    %% update SimMEEG study info
    h.edit_srate.String = num2str(mean(round(1./diff(fn.Time))));
    h.edit_dur.String = sprintf('%.3f %.3f', fn.Time([1 end]));
    h.edit_num_trials.String = num2str(size(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes,1));

    src.Tag = 'bst2sm';
    h.fcn_update_cfg(src,[]);
    h.cfg.study.lat_sim = fn.Time;
    h.btn_reref_EEG.Value = 0;

    %% loading Data: sens_final = [samps x sens x trials]
    dims = size(fn.F);
    chan_idx = 1:size(h.anatomy.sens.chanpos,1);
    h.sim_data.sens_final = nan(dims(2), length(chan_idx), size(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes,1));
    for t=1:length(tidx)
        x = load(file_fullpath(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes(tidx(t)).UserData.Name));
        h.sim_data.sens_final(:,:,t) = x.F(chan_idx,:)'*unit_scale;    % convert from V to microV for SimMEEG
    end
    nan_data = zeros(size(h.sim_data.sens_final));
    h.sim_data.sens_noise_final = nan_data;
    h.sim_data.sens_noise_scaled = nan_data;
    h.sim_data.sens_noise = nan_data;
    h.sim_data.sens_sig_data = nan_data;
    h.sim_data.sens_final_noref = h.sim_data.sens_final;

    %% blanking signal data
    h.sim_data.sig_final = nan_data;
    h.sim_data.sig_final_org = nan_data;
    h.sim_data.sig_wav = nan_data;
    h.sim_data.sig_win = nan_data;
    h.sim_data.prepost_wav = nan_data;
    h.sim_data.prepost_win = nan_data;
    h.sim_data.source_waveform_type = 'nan';

    h.sim_data.cfg = h.cfg;

    %% Updating Import Figure
    h.bst.fig_bst2SimMEEG.study_tree.Enable = 'on';
    h.bst.fig_bst2SimMEEG.btn_import_bst2simmeeg.Enable = 'on';
    h.bst.fig_bst2SimMEEG.txt_select_ds.Text = 'Select Brainstorm Datasets:';
    h.bst.fig_bst2SimMEEG.txt_select_ds.FontColor = 'k';
    % assignin('base','h',h);
    drawnow;
end


end
