function sm_load_sensors_excel(varargin)
global h
[xfile,fpath,indx]= uigetfile({'*.xls; *.xlsx', 'Excel Files (*.xls; *.xlsx)'},'','MultiSelect','off');
xtab = readcell(fullfile(fpath,xfile));
h.sens_table.Data = xtab;

end
