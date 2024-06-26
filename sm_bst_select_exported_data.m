function sm_bst_select_exported_data(varargin)
global h

% Brainstorm's Input files for running analyses
h.bst.datafiles = [];

SubjectNames = fileparts(h.bst.Study.BrainStormSubject);   % returns Subject name
[~, h.bst.condition] = fileparts(h.bst.subj_data_dir);   % returns Subject name

%% Select trial data to compute covariances
h.bst.datafiles = bst_process('CallProcess', 'process_select_files_data', h.bst.datafiles, [], ...
    'subjectname',   SubjectNames, ...
    'condition',     h.bst.condition, ...
    'tag',           'sens_final', ...
    'includebad',    0, ...
    'includeintra',  0, ...
    'includecommon', 0);
    
%% Process: Select avg file if it exists because initial pre-processing already comlpeted. If not then select trial data and create averaged file
h.bst.avgfiles = bst_process('CallProcess', 'process_select_files_data', h.bst.datafiles, [], ...
    'subjectname',   SubjectNames, ...
    'condition',     h.bst.condition, ...
    'tag',           'Avg: 01-sens_final', ...
    'includebad',    0, ...
    'includeintra',  0, ...
    'includecommon', 0);

% h.bst.sourcefiles = bst_process('CallProcess', 'process_select_files_matrix', [], [], ...
%     'subjectname',   SubjectNames, ...
%     'condition',     h.bst.condition, ...
%     'tag',           'sig_final', ...
%     'includebad',    0, ...
%     'includeintra',  0, ...
%     'includecommon', 0);




