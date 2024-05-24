function [data,ref_data]=bl_reref(data,ref_chans,ex_chans,bl_flag)
% re-refereing by simple subtraction of ref_chan data from all other channel data.
% except those excluded in the bl_data.exclude_chans and bad channels
%
%
% INPUT
%   data = [samples x chans x  trials]
%   ref_chans = indices for ref channels
%   ex_chans = indices for excluding channels
%   bl_flag = (1) [default] use BRANE Lab's simple re-referencing (0) use EEGLab's re-referencing

if nargin<4
    bl_flag=1;
end
in_chans=setdiff(1:size(data,2),ex_chans);

if bl_flag==1
% simple re-referncing function --> not using other's programs
ref_data=mean(data(:,ref_chans,:),2);
% simple subtraction - This returns exactly the same as Field Trip's ft_preprocessing.m function.
data(:,in_chans,:)=bsxfun(@minus,data(:,in_chans,:),ref_data);
    
elseif bl_flag==0
    % EEG Lab's rereferencing calculations
    nchansin=length(in_chans);
    if ~isempty(ref_chans) % not average reference
        refmatrix = eye(length(in_chans)); % begin with identity matrix
        for index = 1:length(ref_chans)
            refmatrix(:,ref_chans(index)) = refmatrix(:,ref_chans(index))-1/length(ref_chans);
        end
    else % compute average reference
        refmatrix = eye(nchansin)-ones(nchansin)*1/nchansin;
    end
    chansout = in_chans;
    data2(in_chans,:) = refmatrix*data(in_chans,:);
    data=reshape(data2,[length(ref_chans) size(data,2) size(data,3)]);
    ref_data=nansum(reshape( ((1-refmatrix)*data(in_chans,:)), [length(ref_chans) size(data,2) size(data,3)]));
end


return;






%[bl_data.EEG.data,bl_data.EEG.chanlocs]=reref(bl_data.EEG_org.data, find(bl_data.ref_chans==1),'exclude', [find(bl_data.good_bad_chans==0) find(bl_data.exclude_chans==1)], 'keepref','on','elocs',bl_data.EEG.chanlocs);
%[EEG]=pop_reref(bl_data.EEG, find(bl_data.ref_chans==1),'exclude', [find(bl_data.good_bad_chans==0) find(bl_data.exclude_chans==1)]);
%  guidata(h_fig.BRANE_Lab, handles);
%  update_figure();    %updating figure 
