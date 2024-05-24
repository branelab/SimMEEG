function [rw_idx,cl_idx,perm_mat] = sm_search_min_error_matrix(search_vals,non_val)

% search_vals(1,4) = .5; search_vals(3,4) = .25;
hit_idx = find(~isnan(search_vals)); [row_idx, clmn_idx] = ind2sub(size(search_vals),hit_idx);
%% test for repeats across rows and clmns
row_rep_flag = length(unique(row_idx))~=length(row_idx); % 1=repeats in rows in matrix
clmn_rep_flag = length(unique(clmn_idx))~=length(clmn_idx); % 1=repeats in rows in matrix
[~,w] = unique(row_idx,'stable'); row_rep_idx = setdiff(1:length(row_idx),w,'stable');
row_dup_idx = find(ismember(row_idx,row_idx(row_rep_idx))); % finding duplicates
row_single_idx = find(~ismember(row_idx,row_idx(row_rep_idx))); % finding duplicates

[~,w] = unique(clmn_idx,'stable'); clmn_rep_idx = setdiff(1:length(clmn_idx),w,'stable');
clmn_dup_idx = find(ismember(clmn_idx,clmn_idx(clmn_rep_idx))); % finding duplicates
clmn_single_idx = find(~ismember(clmn_idx,clmn_idx(clmn_rep_idx))); % finding duplicates

hit_single = intersect(row_single_idx, clmn_single_idx); % single hit indices along mat_idx
hit_dup = setdiff(1:length(row_idx),hit_single); % duplicate hit indices along mat_idx
mat_idx = [row_idx clmn_idx]; % row and clmn indices to find permutations with now row or clmn repeats
k=0; % iterating uniquely permutated search_vals with no duplicates to find minimum error.
if row_rep_flag==1 || clmn_rep_flag==1
    idx = sub2ind(size(search_vals), mat_idx(hit_single,1),mat_idx(hit_single,2));
    for rv = 1:length(hit_dup) %iteratively put in duplicate values
        accepted_idx = mat_idx(hit_single,:);    % start new mat_vals
        mat_vals = ones(size(search_vals))*non_val;
        mat_vals(idx) = search_vals(idx);
        didx = sub2ind(size(search_vals), mat_idx(hit_dup(rv),1),mat_idx(hit_dup(rv),2));
        mat_vals(didx) = search_vals(didx);
        accepted_idx = [accepted_idx; mat_idx(hit_dup(rv),:)]; % running indices that are accepted in the updating mat_vals matrix
        
        for dv = 1:length(hit_dup)  % running through duplicates once more to add to current matrix
            if sum(any(mat_idx(hit_dup(dv),:)==accepted_idx))>0 % finding the current duplicate's other duplicate index
%                 fprintf('duplication\n')
            else
                % no duplicates in currently accepted mat_vals matrix.
                didx = sub2ind(size(search_vals), mat_idx(hit_dup(dv),1),mat_idx(hit_dup(dv),2));
                mat_vals(didx) = search_vals(didx);
                accepted_idx = [accepted_idx; mat_idx(hit_dup(dv),:)];
            end
        end
        
        % store mat_vals matrix as a permutation
        k=k+1;
        perm_mat(:,:,k) = mat_vals;
        
    end
else
    perm_mat = search_vals;
    perm_mat(isnan(perm_mat)) = non_val;
end
%% find minimal mse in perm_mat
[~,p] = min(sum(sum(abs(perm_mat),1),2)); % p = index of perm_mat that had lowest mse
% find source index for under column to be true hits
pidx = find(perm_mat(:,:,p)~=non_val); [rw_idx,cl_idx]=ind2sub(size(search_vals),pidx);

