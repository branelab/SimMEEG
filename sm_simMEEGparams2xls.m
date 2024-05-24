function sm_simMEEGparams2xls(varargin)
global h
%% Create an Excel spread sheets with the following parameters in the header line and on separate pages named 'Study' and 'Source'

%% creating clmns names to designate which clmns in xls to store data
alpha = 'A':'Z';
clear clmn*
% iterating columns such that A:Z then AA:ZZ
c=0;
for a=1:length(alpha)
    clmns1(a,:) = sprintf('%s',alpha(a));
    for b=1:length(alpha)
        c = c+1;
        clmns2(c,:) = sprintf('%s%s',alpha(a),alpha(b));
    end
end
clmns = [cellstr(clmns1)' cellstr(clmns2)'];

% creating excel spreadsheets from default h.cfg while SimMEEG_GUI is open
sname = sprintf('%s.xlsx',h.cfg.study.study_name);
delete(sname); % deleting xlsx to overwrite it

sheet_names = {'study' 'source'};

for s = 1:length(sheet_names)
    if s==1; xstruct = h.cfg.study; elseif s==2;  xstruct = h.cfg.source; end
    
    fn = fieldnames(xstruct);
    
    xlswrite(sname,fn',sheet_names{s}); % Writing fieldnames in first row along columns
    
    %% writing out data
    hw = waitbar(0,sprintf('Saving fields for h.cfg.%s',sheet_names{s})); 
    for f=1:length(fn)
    waitbar(f/length(fn),hw,sprintf('Saving fields for h.cfg.%s',sheet_names{s})); 
        cell_start = sprintf('%s2',clmns{f});
        x = xstruct.(fn{f});
        if isnumeric(x)
            if all(max(max(diff(x)))==1) && length(x)>10 % numerical array is a list of sequential integer numbers like 256:512
                x = sprintf('%.f:%.f',x([1 end]));
            else
                x = num2str(x);
            end
        end
        try
            xlswrite(sname,cellstr(x),sheet_names{s},cell_start); % Writing fieldnames in first row along columns
        catch
            fprintf('Field not writeable to xls: %s \npossibly because field has subfield structure\n',fn{f});
            
        end
        
    end
    close(hw);
end


