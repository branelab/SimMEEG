function [ResultsMat,opt] = sm_bst_source_modeling(varargin)
% This function runs source modeling directly within Brainstorm for files curently stored in the h.bst.Protocol folder that were exported from SimMEEG into Brainstrom
global h
%% Select trial data to compute covariances
sm_bst_select_exported_data;

%% Process: No Avg file so computing it and covariances (Noise & Data) based on intervals set within SimMEEG's act_int and ctrl_int
if isempty(h.bst.avgfiles)
    
    % %% SimMEEG simulated data are average referenced already because they are forward projection of source waveforms
    % %% Only implement the following code if real data was loaded into SimMEEG and then exported into Brainstorm
    %     if strcmpi(h.bst.datafiles(1).ChannelTypes,'EEG') % EEG average-ref files do not exists and thus need to create them
    %         %% Process: Apply montage AVERAGE reference for EEG
    %         h.bst.datafiles = bst_process('CallProcess', 'process_montage_apply', h.bst.datafiles, [], ...
    %             'montage',    'Average reference', ...
    %             'createchan', 0);
    %     end
    
    %% Process: Average: Everything
    h.bst.avgfiles = bst_process('CallProcess', 'process_average', h.bst.datafiles, [], ...
        'avgtype',       1, ...  % Everything
        'avg_func',      1, ...  % Arithmetic average:  mean(x)
        'weighted',      0, ...
        'keepevents',    0);
%     %% Process: Average: Everything
%     h.bst.sourcefiles = bst_process('CallProcess', 'process_average', h.bst.sourcefiles, [], ...
%         'avgtype',       1, ...  % Everything
%         'avg_func',      1, ...  % Arithmetic average:  mean(x)
%         'weighted',      0, ...
%         'keepevents',    0);
    
end


%% Process: Compute covariance (noise or data)
h.bst.datafiles = bst_process('CallProcess', 'process_noisecov', h.bst.datafiles, [], ...
    'baseline',       h.inv_soln(h.current_inv_soln).params.ctrl_int*1000, ...
    'datatimewindow', h.inv_soln(h.current_inv_soln).params.act_int*1000, ...
    'sensortypes',    'MEG, EEG, SEEG, ECOG', ...
    'target',         1, ...  % Noise covariance     (covariance over baseline time window)
    'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
    'identity',       0, ...
    'copycond',       0, ...
    'copysubj',       0, ...
    'copymatch',      0, ...
    'replacefile',    1);  % Replace

%% Process: Compute covariance (noise or data)

lat_multiplier = 1; % might need to change this to 1000 for MEG??
display(h.inv_soln(h.current_inv_soln).params.ctrl_int*lat_multiplier);

h.bst.datafiles = bst_process('CallProcess', 'process_noisecov', h.bst.datafiles, [], ...
    'baseline',       h.inv_soln(h.current_inv_soln).params.ctrl_int*lat_multiplier, ...
    'datatimewindow', h.inv_soln(h.current_inv_soln).params.act_int*lat_multiplier, ...
    'sensortypes',    'MEG, EEG, SEEG, ECOG', ...
    'target',         2, ...  % Data covariance      (covariance over data time window)
    'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
    'identity',       0, ...
    'copycond',       0, ...
    'copysubj',       0, ...
    'copymatch',      0, ...
    'replacefile',    1);  % Replace

inv_soln = h.menu_inv_soln.String{h.menu_inv_soln.Value};
%% creating options for inverse solutions
opt = struct('InverseMethod','','InverseMeasure','','SourceOrient','','DataTypes','','Comment','','DisplayMessages',1,...
    'ComputeKernel',1,'Loose',0.2,'UseDepth',[],'WeightExp',[],'WeightLimit',[],'NoiseMethod','','NoiseReg','','SnrMethod','',...
    'SnrRms','','SnrFixed','', 'ChannelTypes','','FunctionName','');
switch inv_soln
    case 'LCMV (BST)'
        opt.InverseMethod = 'lcmv';
        opt.InverseMeasure = 'nai'; % as of Oct-2022, Brainstorm only has the NAI option enabled
        opt.DataTypes = h.menu_sens_type.String(h.menu_sens_type.Value);
        opt.Comment = sprintf('PNAI: %s',opt.DataTypes{:});
        opt.UseDepth = 0;
        opt.WeightExp = 0;
        opt.WeightLimit = 0;
        opt.FunctionName = 'lcmvp';
        opt.SnrMethod = 'rms';
        opt.SnrRms = 1.0000e-06; % not sure where this comes from
        
    case {'MNE (BST)' 'sLORETA (BST)'}
        opt.InverseMethod = 'minnorm';
        opt.DataTypes = h.menu_sens_type.String(h.menu_sens_type.Value);
        switch inv_soln
            case 'MNE (BST)'
                opt.InverseMeasure = 'amplitude';
                opt.Comment = sprintf('MN: %s',opt.DataTypes{:});
                opt.UseDepth = h.radio_inv_bst_depth_weight.Value;
                opt.FunctionName = 'mn';
            case 'sLORETA (BST)'
                opt.InverseMeasure = 'sloreta';
                opt.Comment = sprintf('sLORETA: %s',opt.DataTypes{:});
                opt.UseDepth = 0;
                opt.FunctionName = 'sloreta';
        end
        opt.WeightExp = str2num(h.edit_inv_bst_depth_order.String);
        opt.WeightLimit = str2num(h.edit_inv_bst_depth_max.String);
        opt.SnrMethod = 'fixed';
        opt.SnrRms = 1.0000e-06; % not sure where this comes from
end
opt.SourceOrient = {'free'}; % SimMEEG currently only allows for this option
opt.NoiseMethod = h.cfg.study.bst_params.inv_NoiseMethod{h.menu_inv_bst_reg_noise_type.Value};
if h.menu_inv_bst_reg_noise_type.Value==1; opt.NoiseReg = str2double(h.edit_inv_bst_reg_noise_factor.String); else opt.NoiseReg = 0; end
opt.SnrFixed = 3;  % fixed to 3 dipole orientations




%% Process: Compute sources [2018]
h.bst.datafiles2 = bst_process('CallProcess', 'process_inverse_2018', h.bst.datafiles, [], ...
    'output',  1, ...  % Kernel only: shared
    'inverse', opt);

%% Get Imaging Kernel Data
h.bst.Study = bst_get('Study',h.bst.iStudy);    % updating study data
[tmp,tmp2,iResults] = bst_get('ResultsForDataFile',  h.bst.Study.Data(end).FileName,  h.bst.iStudy);
ResultsFile = h.bst.Study.Result(iResults(end)).FileName;
ResultsMat = in_bst_results(ResultsFile);   % output

%% deleting invSoln Results file so that database doesn't store unnessary files
% deleting datafiles
dfiles = file_fullpath({h.bst.datafiles.FileName}); 
isDeleted = file_delete( dfiles, 1, -1 );
% deleting avgfiles
dfiles = file_fullpath({h.bst.avgfiles.FileName}); 
isDeleted = file_delete( dfiles, 1, -1 );
% deleting noisecov and datacov files
dfiles = file_fullpath({h.bst.Study.NoiseCov.FileName}); 
isDeleted = file_delete( dfiles, 1, -1 );
% deleting result file
dfiles = file_fullpath(ResultsFile); 
isDeleted = file_delete( dfiles, 1, -1 );

% % delete result file 
% dfiles = file_fullpath(h.bst.datafiles2(1).FileName);
% isDeleted = file_delete( dfiles, 1, -1 );

% update 
% db_reload_studies(1:h.bst.iStudy)
db_reload_studies(h.bst.iStudy)

