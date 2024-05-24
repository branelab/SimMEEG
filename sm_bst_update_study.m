function [h] = sm_bst_update_study(h)

%% find study with same sensors and headmodels - set same headmodel as current inverse solution type
chans = bst_get('ChannelForStudy',1:bst_get('StudyCount',h.bst.iProtocol));

for c=1:length(chans)
    if any(strcmpi([chans(c).Modalities],h.anatomy.sens.type)) && strcmpi(h.anatomy.sens.type,'MEG')
        cidx = c;
    elseif any(strcmpi([chans(c).Modalities],h.anatomy.sens.type)) && strcmpi(h.anatomy.sens.type,'EEG') && chans(c).nbChannels == length(h.anatomy.sens.label)
        cidx = c;
    end
end

%% get Study and Subj info
if isfield(h.bst,'subj_MriFile') % SimMEEG called from Brainstorm
    x = dir(sprintf('%s\\*.mat',h.bst.subj_data_dir));
    [h.bst.Study, h.bst.iStudy] = bst_get('AnyFile',fullfile(h.bst.subj_data_dir,x(1).name));
    SubjectNames = fileparts(h.bst.Study.BrainStormSubject);   % returns Subject name
    [h.bst.sSubject, h.bst.iSubject] = bst_get('Subject', SubjectNames);
else
    [h.bst.Study, h.bst.iStudy] = bst_get('AnyFile',chans(cidx).FileName);
    SubjectNames = fileparts(h.bst.Study.BrainStormSubject);   % returns Subject name
    [h.bst.sSubject, h.bst.iSubject] = bst_get('Subject', SubjectNames);
    
    ProtocolInfo = bst_get('ProtocolInfo');
    h.bst.subj_data_dir = bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(h.bst.Study.FileName));
end




%% Headmodel Type

[h.bst.sSubj, iSubj]  = bst_get('Study', h.bst.iSubject);
bst_hdm_types = {h.bst.sSubj.HeadModel.HeadModelType};

switch h.inv_soln(h.current_inv_soln).headmodel_type
    case 'Whole Brain'
        iHeadModel = find(strcmpi(bst_hdm_types,'volume'));
    case 'Cortical Surface'
        iHeadModel = find(strcmpi(bst_hdm_types,'surface'));
end
% h.bst.Study.iHeadModel = iHeadModel;
% bst_set('Study', h.bst.iStudy, h.bst.Study);
h.bst.sSubj.iHeadModel = 2;
bst_set('Study', 1, h.bst.sSubj);
% fprintf('Changed to BrainsStormSubject: %s \nChanged Headmodel to %s\n',h.bst.Study.BrainStormSubject, h.bst.sSubj.HeadModel(h.bst.sSubj.iHeadModel).HeadModelType);
% fprintf('Changed to BrainsStormSubject: %s \nChanged Headmodel to %s\n',h.bst.sSubj.BrainStormSubject, h.bst.sSubj.HeadModel(h.bst.sSubj.iHeadModel).HeadModelType);

