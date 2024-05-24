function [data,lat] = sm_load_bst_trial_data(datadir)
% function [data,lat] = sm_load_bst_trial_data(datadir)
% This program will load the trials data store in the Brainstorm data directory (datadir) into a single matrix [samps x chans x trials] to be used within SimMEEG/BRANE lab software
% 

% datadir = 'F:\BackUp_Drive\brainstorm_db\SimMEEG_test\data\Subject01\concatenated_Average_reference';
cur_dir = pwd;
cd(datadir); 
xdir = dir('data*trial*.mat'); 
x = struct2cell(xdir); 
idx = find(contains(x(1,:),'data')); % finding only files with data listed 
cd(cur_dir);

for t = 1:length(idx)
    xd = load(fullfile(datadir,xdir(idx(t)).name));
    data(:,:,t) = xd.F';
end
lat = xd.Time;



