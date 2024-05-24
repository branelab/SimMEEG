function [rw_idx,cl_idx] = sm_search_min_error_matrix(search_vals)

%% find first lowest value (min) then blank across rows (peaks) and columns (true sources), then search for min again until no more points to search .
x = search_vals; 

% starting point for search
rw_idx =[]; cl_idx=[];
v=0; 
while any(~isnan(x(:)))
v=v+1;
    [q] = nanmin(x(:)); 
    [row,col] = find(x==q);
    if length(row)>1; row = row(1); col = col(1); end % selecting first instance of minima when there are two equal minima like = 0 
    rw_idx(v) = row; cl_idx(v) = col; 
    % blanking out rows and columns for each min found
    x(row,:) = nan;
    x(:,col) = nan;
end


