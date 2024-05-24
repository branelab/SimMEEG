function sm_load_source_locs(varargin)
% load source locations from .mat file with variables "source_locs" [sources x XYZ] and "source_ori" [sources x Xori, Yori, Zori]
global h
[xfile,fpath,indx]= uigetfile({'*.mat; *.mat', 'Matlab Files (*.mat)'},'','MultiSelect','off');

fname = fullfile(fpath,xfile); 
listOfVariables = who('-file', fname);

%% Source Locations
if ismember('source_locs', listOfVariables) % loading source locs
    load(fname,'source_locs');
    h.cfg.source.vx_locs = source_locs;
    for v=1:3; h.edit_source_locs(v).String = num2str(round(source_locs(v,:))); end
else
    warndlg(sprintf('Source locations do not exist\nFile: %s',fname));
end

%% Source Orientations
if ismember('source_ori', listOfVariables) % loading source locs
    load(fname,'source_ori');
    h.cfg.source.vx_ori = source_ori;
    
    [az,el] = cart2sph(source_ori(:,1),source_ori(:,2),source_ori(:,3))
    az = rad2deg(az); el = rad2deg(el);
    
    for v=1:3; h.edit_source_ori(v).String = sprintf('%.f %.f',az(v),el(v)); end
    
else
    warndlg(sprintf('Source orientations do not exist\nFile: %s',fname));
end
