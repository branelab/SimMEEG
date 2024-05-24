function sm_set_h_properties(h,set_props,h2)
% This program sets properties in "h" from "h2" for those listed set_props 
% It can be used to create a fake figure handle "h" for use with sm_batch scripts
% INPUT:
%   h = handle to figure
%   set_props = {'String' 'Value'};
%
% OUTPUT:
%   h2 = dummy variable that has properties of keep_props but nothing else.


% h2 = h.monte_params.h;

fn = fieldnames(h2);
for f=1:length(fn)
    for kp=1:length(set_props)
        if isfield(h,fn{f})
            if ~isempty( findobj(h.(fn{f}),'-property',set_props{kp},'-depth',0) )
%                 try
                for m = 1:length(h.(fn{f}))
                    eval(sprintf('h.%s(%.f).%s = h2.%s(%.f).%s;',(fn{f}),m,set_props{kp},(fn{f}),m,set_props{kp}));
                    try
                    if strcmp(h.(fn{f})(m).Style,'slider')     % need to make sure slider value is within limits of min_max so setting 
                        if h.(fn{f}).Value < h.(fn{f}).Min;  h.(fn{f}).Value = h.(fn{f}).Min; end
                        if h.(fn{f}).Value > h.(fn{f}).Max;  h.(fn{f}).Value = h.(fn{f}).Max; end
                    end
                    catch
%                     fprintf('"Style" property does not exist for h.%s\n',fn{f})
                    end
                end
%                 catch
%                     fprintf('Eval did not work for:\n%s\n',sprintf('h.%s(%.f).%s = h2.%s(%.f).%s;',(fn{f}),m,set_props{kp},(fn{f}),m,set_props{kp}));
%                 end
            end
        end
        
    end
end


