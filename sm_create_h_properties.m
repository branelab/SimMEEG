function h2 = sm_create_h_properties(h,keep_props)
% This program strips all properties in "h" while keeping on the ones listed in keep_props
% It can be used to create a fake figure handle "h" for use with sm_batch scripts
% INPUT:
%   h = handle to figure
%   keep_props = {'String' 'Value'};
%
% OUTPUT:
%   h2 = dummy variable that has properties of keep_props but nothing else.
h2 = struct('dummy_handle','');

%% Finding all fields with at least one of the keep_props
fn = fieldnames(h);
% xidx = zeros(length(fn),1);
for f=1:length(fn)
    if isobject(h.(fn{f}))
        for kp=1:length(keep_props)
            if ~isempty( findobj(h.(fn{f}),'-property',keep_props{kp},'-depth',0) )
                for m = 1:length(h.(fn{f}))
                    if isvalid(h.(fn{f})(m))
                        try
                            eval(sprintf('h2.%s(%.f).%s = h.%s(%.f).%s;',(fn{f}),m,keep_props{kp},(fn{f}),m,keep_props{kp}));
                        catch
                        
                        end

                    end
                end
            end
        end
    end
end




