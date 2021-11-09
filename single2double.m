function ds = single2double(ds)
% ds = data structure


fnames = fieldnames(ds);
for f=1:length(fnames)
    if eval(sprintf('isa(ds.%s,''single'')',fnames{f}))
        eval(sprintf('ds.%s = double(ds.%s);',fnames{f},fnames{f}))
    end
end

 
