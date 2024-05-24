function sm_export2brainstorm(varargin)
global h

[h] = sm_open_bst_study(h);  % opens study in brainstorm, if ProtocolName doesn't exist then it creates a new Protocol and Study.

%%
iStudy = h.bst.iStudy;
% Get study
sStudy = bst_get('Study', iStudy);

h.bst.Study = bst_get('Study',h.bst.iStudy);
%% get averaged filename if exists, if not then export data
SubjectNames = fileparts(h.bst.Study.BrainStormSubject);   % returns Subject name

%% deleting previous data the re-exporting them
if h.radio_overwrite_bst_data.Value == 1
    for n=1:length(h.bst.datafiles)
        fname = fullfile(h.bst.ProtocolInfo.STUDIES, h.bst.datafiles(n).FileName);
        delete(fname);
    end
else
end


% Update database
bst_set('Study', iStudy, sStudy);
panel_protocols('UpdateNode', 'Study', iStudy);

%% Export to Brainstorm
% Returned files
NewFiles = {};
% No simulated data
if isempty(h) || ~isfield(h, 'sim_data') || isempty(h.sim_data)
    return;
end
% Get available fields
allFields = {'sens_final', 'sig_final', 'sig_wav', 'prepost_wav', 'noise_wav', 'prepost_win', 'sig_win', ...
    'sens_noise', 'sens_noise_scaled', 'sens_sig_data', 'sens_noise_final', 'sens_final_org', 'sens_noise_final_org', 'sens_sig_data_org'};
availFields = {};
for i = 1:length(allFields)
    if isfield(h.sim_data, allFields{i}) && ~isempty(h.sim_data.(allFields{i}))
        availFields{end+1} = allFields{i};
    end
end
if isempty(availFields)
    return;
end
% Ask fields to import
% isSelect = ismember(availFields, {'sens_final', 'sig_final'});
isSelect = ismember(availFields, {'sens_final'});
if ~any(isSelect)
    isSelect = strcmpi(availFields, 'sig_final');
end

if h.bst.gui_flag == 1  % if running using gui
    isSelect = java_dialog('checkbox', 'Select variables to import:', 'Import SimMEEG output', [], availFields, isSelect);
    if isempty(isSelect)
        return;
    end
else
    
end

importGroups = availFields(isSelect == 1);
% Timestamp
c = clock;
strTime = sprintf('_%02.0f%02.0f%02.0f_%02.0f%02.0f', c(1)-2000, c(2:5));
% Progress bar
nTrials = length(importGroups) * size(h.sim_data.(importGroups{1}),3);
bst_progress('start', 'SimMEEG', 'Importing simulated files...', 0, nTrials);
% Detect previous executions of simmeeg in this study
iRun = 1;
if ~isempty(sStudy.Matrix)
    while any(cellfun(@(c)and(length(c)>3, strcmpi(c(1:3),sprintf('%02d-',iRun))), {sStudy.Matrix.Comment}))
        iRun = iRun + 1;
    end
end
% Load channel file
if isempty(sStudy.Channel)  % If using @default_anat
    sStudy.Channel = bst_get('ChannelForStudy',h.bst.iStudy);
    sStudy.iChannel = 1;
end
if ~isempty(sStudy.Channel)
    ChannelMat = in_bst_channel(sStudy.Channel.FileName);
    nChannels = length(ChannelMat.Channel);
    try
    iChannels = cellfun(@(c)find(strcmpi(c,{ChannelMat.Channel.Name}),1), h.anatomy.sens.label);
    catch
        for v=1:nChannels
            iChannels(v) = contains(h.anatomy.sens.label{v},ChannelMat.Channel(v).Name);
        end
    end

end
% Loop to import groups of trials
for iGroup = 1:length(importGroups)
    % List of trials to import
    trials = h.sim_data.(importGroups{iGroup});
    nSources = size(trials,2);
    % Save matrix
    if ismember(importGroups{iGroup}, {'sig_final', 'sig_wav', 'prepost_wav', 'noise_wav', 'prepost_win', 'sig_win'})
        % Create a "matrix" structure
        sMat = db_template('matrixmat');
        sMat.Time        = h.sim_data.cfg.study.lat_sim;
        sMat.Description = cell(nSources,1);
        for iSource = 1:nSources
            sMat.Description{iSource} = sprintf('Source %d', iSource);
        end
        % Add history entry
        sMat = bst_history('add', sMat, 'process', 'Generated with SimMEEG');
        % Add extra cfg structure
        sMat.cfg = h.sim_data.cfg;
        % Loop on each trial
        for iTrial = 1:size(trials,3)
            % Trial values
            sMat.Comment = sprintf('%02d-%s (#%d)', iRun, importGroups{iGroup}, iTrial);
            sMat.Value   = trials(:,:,iTrial)';
            % Output filename
            %                 FileName = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ['matrix_simmeeg', strTime, '_', strrep(importGroups{iGroup}, '_', '-'), sprintf('_trial%03d', iTrial)], 0);
            FileName = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ['matrix_simmeeg', '_', strrep(importGroups{iGroup}, '_', '-'), sprintf('_trial%03d', iTrial)], 0);
            fpath = fileparts(FileName);
            FileName = fullfile(fpath, ['data_simmeeg', '_', strrep(importGroups{iGroup}, '_', '-'), sprintf('_trial%03d.mat', iTrial)]);
            % Save on disk
            bst_save(FileName, sMat, 'v6');
            NewFiles{end+1} = FileName;
            % Add structure to database
            sNew = db_template('Matrix');
            sNew.FileName = file_short(FileName);
            sNew.Comment  = sMat.Comment;
            iItem = length(sStudy.Matrix) + 1;
            sStudy.Matrix(iItem) = sNew;
            % Increment progress bar
            bst_progress('inc',1);
        end
        % Save MEG/EEG data
    else
        % Create a "matrix" structure
        sMat = db_template('datamat');
        sMat.Time = h.sim_data.cfg.study.lat_sim;
        sMat.ChannelFlag = ones(nChannels,1);
        sMat.Device = 'SimMEEG';
        % Add history entry
        sMat = bst_history('add', sMat, 'process', 'Generated with SimMEEG');
        % Add extra cfg structure
        sMat.cfg = h.sim_data.cfg;
        % Loop on each trial
        for iTrial = 1:size(trials,3)
            % Trial values
            sMat.Comment = sprintf('%02d-%s (#%d)', iRun, importGroups{iGroup}, iTrial);
            sMat.F = zeros(nChannels, size(trials,1));
            sMat.F(iChannels,:) = trials(:,:,iTrial)';
            % Output filename
            %                 FileName = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ['data_simmeeg', strTime, '_', strrep(importGroups{iGroup}, '_', '-'), sprintf('_trial%03d', iTrial)], 0);
            FileName = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ['data_simmeeg', '_', strrep(importGroups{iGroup}, '_', '-'), sprintf('_trial%03d', iTrial)], 0);
            fpath = fileparts(FileName);
            FileName = fullfile(fpath, ['data_simmeeg', '_', strrep(importGroups{iGroup}, '_', '-'), sprintf('_trial%03d.mat', iTrial)]);
            
            % Save on disk
            bst_save(FileName, sMat, 'v6');
            NewFiles{end+1} = FileName;
            % Add structure to database
            sNew = db_template('Data');
            sNew.FileName = file_short(FileName);
            sNew.Comment  = sMat.Comment;
            iItem = length(sStudy.Data) + 1;
            sStudy.Data(iItem) = sNew;
            % Increment progress bar
            bst_progress('inc',1);
        end
    end
end
% Update database
bst_set('Study', iStudy, sStudy);
panel_protocols('UpdateNode', 'Study', iStudy);
db_reload_studies( h.bst.iStudy )
% Close progress bar
bst_progress('stop');


end

