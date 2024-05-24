function sm_bst_Import2SimMEEG(varargin)

global h

%% get Studies in Protocol
h.bst.Protocol = bst_get('ProtocolInfo');

%% delete old Figure if present
if isfield(h.bst,'fig_bst2SimMEEG')
    if isvalid(h.bst.fig_bst2SimMEEG.main_fig)
        close(h.bst.fig_bst2SimMEEG.main_fig); h.bst = rmfield(h.bst,'fig_bst2SimMEEG');
    end
end

%% Create new Figure
h.bst.fig_bst2SimMEEG.main_fig = uifigure;
h.bst.fig_bst2SimMEEG.main_fig.MenuBar = 'none';
h.bst.fig_bst2SimMEEG.main_fig.Name = 'Import from Brainstorm to SimMEEG';
fig_pos = h.bst.fig_bst2SimMEEG.main_fig.Position;
main_pos = [fig_pos(3:4) fig_pos(3:4)];    % position for normalizing pixels using
%% Text: Select Datasets
h.bst.fig_bst2SimMEEG.txt_select_ds = uilabel(h.bst.fig_bst2SimMEEG.main_fig,'Position',[.01 .91 .49 .05].*main_pos, ...
    'Text','Select Brainstorm Datasets:','FontSize',14);
%% ListBox: Study Tree
% delete(h.bst.fig_bst2SimMEEG.study_tree);
h.bst.fig_bst2SimMEEG.study_tree = uitree(h.bst.fig_bst2SimMEEG.main_fig,'Position',[.01 .01 .49 .9].*main_pos);
h.bst.fig_bst2SimMEEG.study_tree.SelectionChangedFcn = @(src, event)select_all(src, event);
h.bst.fig_bst2SimMEEG.Study.study_node = uitreenode(h.bst.fig_bst2SimMEEG.study_tree,'Text',h.bst.Protocol.Comment,'NodeData',[]);
%% Btn:
% delete(h.bst.fig_bst2SimMEEG.btn_import_bst2simmeeg)
h.bst.fig_bst2SimMEEG.btn_import_bst2simmeeg = uibutton(h.bst.fig_bst2SimMEEG.main_fig, ...
    "Text","Import Data", "Position",[.51 .91 .15 .05].*main_pos, "FontColor", [0 .5 0],...
    "BackgroundColor", [.9 1 .9], ...
    "ButtonPushedFcn",@(src,event)sm_bst_import_bst2simmeeg);



%% 
datadir = h.bst.Protocol.STUDIES;
x = dir(datadir); x = x([x.isdir]); 
subj = x(~contains({x.name},'.'));

for n=1:length(subj)
    %% Subj
    h.bst.Study.subj(n).name = subj(n).name;
    h.bst.Study.subj(n).datadir = fullfile(h.bst.Protocol.STUDIES, subj(n).name);
    %% Conditions within subjdir
    x = dir(h.bst.Study.subj(n).datadir); x = x([x.isdir]);
    condlist = x(~contains({x.name},'.'));
    if isempty(condlist)    % search for datafile within subjdir
        h.bst.Study.subj(n).cond(1).name = 'none';
        h.bst.Study.subj(n).cond(1).datadir = h.bst.Study.subj(n).datadir;
        x = dir(sprintf('%s\\data_*.mat',h.bst.Study.subj(n).cond(1).datadir)); 
        %% Need to parse out datasets with different names within this diretory
        if isempty(x)
            h.bst.Study.subj(n).cond(1).datasets(1) = 'none';
        else
            for t=1:length(x)
                h.bst.Study.subj(n).cond(t).datasets(1) = {x(t).name};
            end
        end
    else
        
        for c=1:length(condlist)
            h.bst.Study.subj(n).cond(c).name = condlist(c).name;
            h.bst.Study.subj(n).cond(c).datadir = fullfile(h.bst.Protocol.STUDIES, subj(n).name, h.bst.Study.subj(n).cond(c).name);
            %% DataSets within condition folder
            datalist = dir(sprintf('%s\\data_*.mat',h.bst.Study.subj(n).cond(c).datadir));
            if isempty(datalist)
                h.bst.Study.subj(n).cond(c).dataset(1).name = 'none';
                h.bst.Study.subj(n).cond(c).dataset(1).trials(1) = {'none'};
            else
                %%
                
                
                
                %% get dataset filenames
                for d=1:length(datalist)
                    h.bst.Study.subj(n).cond(c).dataset(d).name = [];
                    h.bst.Study.subj(n).cond(c).dataset(d).trials(t) = [];
                end
                
            end
            
        end
        
        
        
    end
    
    h.bst.Study.subj(n).cond
    
end


%% Create tree

if isempty(subj)
    h.bst.fig_bst2SimMEEG.Study.subj_node(1) = uitreenode(h.bst.fig_bst2SimMEEG.Study.study_node,'Text', 'No Data in Brainstorm Protocol', 'NodeData', [1]);
else
    
    [subj, xidx] = sort(subj(sidx));
    pidx = sidx(xidx);
    
    pStudy = h.bst.ProtocolStudy.Study(pidx);
    
    n=0;
    %% Condition Level
    for nx=1:size(pStudy,2)
        h.bst.fig_bst2SimMEEG.Study.subj_node(nx) = uitreenode(h.bst.fig_bst2SimMEEG.Study.study_node,'Text', subj{nx}, 'NodeData', [nx]);
        
        if isempty(pStudy(nx).Data)
            %    h.bst.fig_bst2SimMEEG.Study.cond_node(n) = uitreenode(h.bst.fig_bst2SimMEEG.Study.study_node,'Text', '' , 'NodeData', [n]);
        else
            n = n+1;
            str = pStudy(nx).FileName; 
            bidx = findstr(str,'/'); 
            c_cond = str(bidx(1)+1:bidx(2)-1);
            h.bst.fig_bst2SimMEEG.Study.cond_node(n) = uitreenode(h.bst.fig_bst2SimMEEG.Study.subj_node(nx),'Text', c_cond, 'NodeData', [n]);
            %% find unique dataset names
            dsc = {pStudy(nx).Data.Comment};
            ds = struct('Name','','FileName','','idx',[]);
            for d=1:size(dsc,2)  % getting start of ds names
                str = pStudy(nx).Data(d).Comment;
                bidx = findstr(str,' '); 
                ds(d).Name = str(1:bidx(1)-1);
                ds(d).FileName = pStudy(nx).Data(d).FileName;
                ds(d).idx = d;
            end
            %% unique dataset names
            ds_txt = unique({ds.Name});
            for m=1:length(ds_txt)  % number of datasets within a study file "cond_node"
                h.bst.fig_bst2SimMEEG.Study.Datasets(n,m) = uitreenode(h.bst.fig_bst2SimMEEG.Study.cond_node(n),'Text', ds_txt{m}, 'NodeData', [n m]);
                ds_idx(m).idx = find(strcmp( {ds.Name}, ds_txt{m}));
                for t=1:length(ds_idx(m).idx)
                    h.bst.fig_bst2SimMEEG.Study.Dataset_Trials(n,m).trials(t) = uitreenode(h.bst.fig_bst2SimMEEG.Study.Datasets(n,m),...
                        'Text', pStudy(nx).Data(ds_idx(m).idx(t)).Comment, 'NodeData', [n m t], 'UserData', ds(t).FileName);
                end
            end
        end
    end
    for n=1:length(h.bst.fig_bst2SimMEEG.Study.cond_node)
        expand(h.bst.fig_bst2SimMEEG.Study.cond_node(n));
    end
end
expand(h.bst.fig_bst2SimMEEG.study_tree); drawnow;


h.bst.fig_bst2SimMEEG.study_tree.Multiselect = "on";



end

function select_all(varargin)
% select all nodes underneath varargin source
global h

sel_nodes = varargin{2}.SelectedNodes;
if length(sel_nodes)==1     % only one node selected
    % select all trial nodes below
    if length(sel_nodes.NodeData)==2    % ordinate dataset node
        expand(sel_nodes);
        h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes = h.bst.fig_bst2SimMEEG.Study.Dataset_Trials(sel_nodes.NodeData(1), sel_nodes.NodeData(2)).trials;
    end
else
end


end


function sm_bst_import_bst2simmeeg(varargin)
% imports selected files from BST into SimMEEG as "sens_final" and converts to SimMEEG format
global h
if isempty(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes)
    warndlg('Please select Brainstorm datasets', 'No Datasets Selected');
else
    %% overwriting h.sim_data
    h.sim_data = struct('sig_final', [], 'sig_wav', [],'prepost_wav', [], 'noise_wav', [], 'cfg', [], 'prepost_win', [], ...
        'sig_win', [], 'source_waveform_type', [], 'intersource_correlations', [], 'sens_noise_type', 'Brainstorm Imported', 'sens_noise', [], 'sens_noise_scaled', [], ...
        'sens_sig_data', [], 'sens_final', [], 'noise_gain', [], 'sens_noise_final', []);
    
    %% resetting SimMEEG study parameters
    % load data
    fn = load(file_fullpath(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes(1).UserData));
    
    % get channel info
    cfile = file_fullpath(bst_get('ChannelFileForStudy', file_fullpath(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes(1).UserData)));
    chans = bst_get('ChannelFile',cfile);
    if any(contains( chans.Channel.Modalities, 'MEG')) || any(contains( chans.Channel.Modalities, 'EEG'))
    else
        warndlg('Sensors file does not contain "MEG" or "EEG" channels','No MEG/EEG sensors');
    end
    %% rescaling units
    % check sensor units V, microV, ... SimMEEG needs to be in femtoT
    if mean(std(fn.F,[],2))>1   % femto
        unit_scale = 1;   % convert to femtoT
    elseif mean(std(fn.F,[],2))<1e-3 && mean(std(fn.F,[],2))<1  % milli
        unit_scale = 1e3;   % convert to femtoT
    elseif mean(std(fn.F,[],2))<1e-6 && mean(std(fn.F,[],2))<1e-3  % micro
        unit_scale = 1e6;  % convert to femtoT
    elseif mean(std(fn.F,[],2))<1e-9 && mean(std(fn.F,[],2))<1e-6  % nan
        unit_scale = 1e9;  % convert to femtoT
    elseif mean(std(fn.F,[],2))<1e-12 && mean(std(fn.F,[],2))<1e-9  % pico
        unit_scale = 1e12;  % convert to femtoT
    elseif mean(std(fn.F,[],2))<1e-15 && mean(std(fn.F,[],2))<1e-12  % femto
        unit_scale = 1e15;  % convert to femtoT
    end
    
    % update SimMEEG study info
    h.edit_srate.String = num2str(mean(round(1./diff(fn.Time))));
    h.edit_dur.String = num2str(fn.Time([1 end]));
    h.edit_num_trials.String = num2str(size(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes,1));
    
    h.fcn_update_cfg();
    h.cfg.study.lat_sim = fn.Time;
    h.btn_reref_EEG.Value = 0;
    
    %% loading Data: sens_final = [samps x sens x trials]
    dims = size(fn.F);
    chan_idx = 1:size(h.anatomy.sens.chanpos,1);
    h.sim_data.sens_final = nan(dims(2), length(chan_idx), size(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes,1));
    for t=1:size(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes,1)
        x = load( file_fullpath(h.bst.fig_bst2SimMEEG.study_tree.SelectedNodes(t).UserData) );
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
    
end


end

%
%
%
% % Brainstorm's Input files for running analyses
% h.bst.import_files= [];
%
% SubjectNames = fileparts(h.bst.Study.BrainStormSubject);   % returns Subject name
% [~, h.bst.condition] = fileparts(h.bst.subj_data_dir);   % returns Subject name
%
% %% Select trial data to compute covariances
% h.bst.import_files = bst_process('CallProcess', 'process_select_files_data', h.bst.import_files, [], ...
%     'subjectname',   SubjectNames, ...
%     'condition',     h.bst.condition, ...
%     'tag',           'sens_final', ...
%     'includebad',    0, ...
%     'includeintra',  0, ...
%     'includecommon', 0);
