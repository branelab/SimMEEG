function results = sm_batch_load_results(varargin)

data = varargin{1}; 
field_names = varargin{2};

for a = 1:length(data)
    for f=1:length(field_names)
        
        fname = field_names{f};
        try
            eval( sprintf('results(a).%s = data(a).%s;',fname,fname) );
        catch me
            fprintf('Error in loading results from %s\n%s\n',fname,me.message);
        end
    end
    
end

