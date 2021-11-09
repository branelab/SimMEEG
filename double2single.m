function ds = double2single(ds)
% ds = data structure


for n=1:length(ds)
 fnames = fieldnames(ds(n));
   for f=1:length(fnames)
        % diving down one more level of the structure 
        if eval(sprintf('isstruct(ds(%.f).%s)',n,fnames{f}))
            xds = eval(sprintf('ds(%.f).%s',n,fnames{f}));
            for xn=1:length(xds)
                xfnames = fieldnames(xds(xn));
                for xf=1:length(xfnames)
                    % diving down one more level of the structure
                    if eval(sprintf('isstruct(xds(%.f).%s)',xn,xfnames{xf}))        %% is struct
                        yds = eval(sprintf('xds(%.f).%s',xn,xfnames{xf}));
                        
                        for yn = 1:length(yds)
                            yfnames = fieldnames(yds(yn));
                            for yf=1:length(yfnames)
                                if eval(sprintf('isa(yds(%.f).%s,''double'')',yn,yfnames{yf}))
                                    eval(sprintf('yds(%.f).%s = single(yds(%.f).%s);',yn,yfnames{yf},yn,yfnames{yf}))
                                end
                            end
                        end
                        
                    elseif eval(sprintf('iscell(xds(%.f).%s)',xn,xfnames{xf}))      %% is cell
                        yds = eval(sprintf('xds(%.f).%s',xn,xfnames{xf}));
                        
                        for yn = 1:length(yds)
                           yds{yn} = single(yds{yn});
                        end
                        
                    else
                        
                        if eval(sprintf('isa(xds(%.f).%s,''double'')',xn,xfnames{xf}))
                            eval(sprintf('xds(%.f).%s = single(xds(%.f).%s);',xn,xfnames{xf},xn,xfnames{xf}))
                        end
                        
                    end
                end
            end
        elseif  eval(sprintf('iscell(ds(%.f).%s)',n,fnames{f}))      %% is cell
           
            xds = eval(sprintf('ds(%.f).%s',n,fnames{f}));
            
            for xn = 1:length(xds)
                xds{xn} = single(xds{xn});
            end
        else
            if eval(sprintf('isa(ds(%.f).%s,''double'')',n,fnames{f}))
                eval(sprintf('ds(%.f).%s = single(ds(%.f).%s);',n,fnames{f},n,fnames{f}))
            end
        end
    end
end

 
