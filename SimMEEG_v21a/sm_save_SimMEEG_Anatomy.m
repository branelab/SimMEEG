function sm_save_SimMEEG_Anatomy(varargin)
global h

h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Saving SimMEEG Anatomy:'); drawnow;
[~,fname,~] = fileparts(h.anat_file);
fname = sprintf('%s_SimMEEG_Anatomy.mat',fname);
sname = fullfile(h.anat_path, fname);
[fname,fpath]=uiputfile(sname,'Save SimMEEG Anatomy');

if fname~=0
    h.anat_path=fpath; h.anat_file=fname;
    sxname = fullfile(fpath,fname);
    h.waitfor_txt.String = sprintf('Saving SimMEEG Anatomy:\n\n%s',sxname); drawnow;
    anat = h.anatomy;
    save(sxname,'-struct','anat');
else
end
update_anatomy_fields;
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');


