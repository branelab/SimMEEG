function sm_create_plv_contrasts(varargin)
global h
% 
% h.radio_3D_plot_plv_locs.Value = 0; sm_plot_plv_locs; % turning off leadfield and deleting from 3D MRI plot

%% no peak sources found
if all(isnan(h.inv_soln(h.current_inv_soln).classifier_metrics.Hits))
    h.inv_soln(h.current_inv_soln).plv_seed_idx = h.inv_soln(h.current_inv_soln).classifier_metrics.Hits;
    hm = warndlg(sprintf('No Peak Sources Found\nChange Image Scale --> Select Seed and Comparions Locations'),'No Sources Found '); WinOnTop(hm);
    hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
    return
end

if ~isfield(h.inv_soln(h.current_inv_soln),'plv_seed_idx') && ~isfield(h.inv_soln(h.current_inv_soln),'plv_comp_idx')
    hm = warndlg(sprintf('Please select Seed and Comparison Locations\n'),'Select Seed and Comparions Locations'); WinOnTop(hm);
    hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
    h.inv_soln(h.current_inv_soln).plv_contrasts = [];
else
    %% PLV contrasts are nchoosek comparisons among seed_idx and comp_idx
    t=0; plv_contrasts=[];
    seed_idx = h.inv_soln(h.current_inv_soln).plv_seed_idx;
    comp_idx = h.inv_soln(h.current_inv_soln).plv_comp_idx;
    
    num_contrasts = length(seed_idx)*length(comp_idx);
    if num_contrasts > 1e4
        answ = questdlg(sprintf('Number of PLV contrasts is %.f\n\nThis will take a long time for PLV to compute\nand may exceed computer''s memory capacity\n\nWould you like to continue?\n',num_contrasts),'Save SimMEEG Dataset?','Yes','No','No');
    else
        answ = 'Yes';
    end
    
    switch answ
        case 'Yes'
%             hm = msgbox(sprintf('Generating PLV contrasts.\n\nThis may take time if large Lead Field Grid selected.\n'),'Generating PLV/PLI contrasts'); WinOnTop(hm);
%             hm.Units = 'normalized';
%             hm.Position(1:2) = [h.main_fig.Position(1)+(sum(h.main_fig.Position([1 3]))/3) h.main_fig.Position(2)+(sum(h.main_fig.Position([2 4]))/2)];
%             hm.Units = 'normalized'; hm.Position = [sum(h.main_fig.Position([1 3]))/2 sum(h.main_fig.Position([2 4]))/2 .2 .1];
            
            
            if length(h.inv_soln(h.current_inv_soln).plv_seed_idx) == length(h.inv_soln(h.current_inv_soln).plv_comp_idx) ...
                    && isempty( setdiff(h.inv_soln(h.current_inv_soln).plv_seed_idx,h.inv_soln(h.current_inv_soln).plv_comp_idx) )
                % same voxels in seed and comps thus nchoose2 calculation
                cx_idx = nchoose2(1:length(h.inv_soln(h.current_inv_soln).plv_seed_idx));
                h.inv_soln(h.current_inv_soln).plv_contrast_idx = cx_idx;
                plv_contrasts = h.inv_soln(h.current_inv_soln).plv_seed_idx(cx_idx);
            else    % different voxels in seed and comparison lists
                % All seeds to all comps constrats
                plv_contrasts = []; cx_idx2=[];
                for s = 1:length(seed_idx)
                    % remove seed from comps
                    [cx_comp_idx, cc_idx] = setdiff(comp_idx,seed_idx(s));
                    
                    try          cx_idx = cat(2,repmat(seed_idx(s),length(cx_comp_idx),1), cx_comp_idx');
                    catch        cx_idx = cat(2,repmat(seed_idx(s),length(cx_comp_idx),1), cx_comp_idx);
                    end
                    % concatenate seed-to-comp contrasts
                    plv_contrasts = cat(1,plv_contrasts,cx_idx);
                    cx_idx2 = cat(1, cx_idx2, [repmat(s,length(cc_idx),1), cc_idx]);
                end
                % corresponds to indices of the plv_contrasts that are needed when calculating source waveforms (swf) in sm_calc_PLV_PLI.m
                % clmn 1 = seed_idx;  clmn 2 = comp_idx
                
                % remove duplicates [2 1], [3 1], ... that already exist as [1 2], [1 3],...  
                % remove same idx [1 1] [2 2] when there are misses
                rm_s = nan(length(cx_idx2),1);
                seed_vals = nchoose2(1:length(seed_idx));
                for s = 1:size(seed_vals,1)
                    if ~isempty(find(ismember( cx_idx2, fliplr(seed_vals(s,:)),'rows' ), 1))
                    rm_s(s) = find(ismember( cx_idx2, fliplr(seed_vals(s,:)),'rows' ));
                    end
                end
                rm_s(cx_idx2(:,1)-cx_idx2(:,2)==0)=1; % removing same [1 1] [2 2] when there are misses
                
                rm_s = rm_s(~isnan(rm_s)); % cx_idx2(rm_s,:)
                good_idx = setdiff(1:size(cx_idx2,1),rm_s);
                cx_idx2 = cx_idx2(good_idx,:); % final plv_contrast_idx with duplicates removed
                h.inv_soln(h.current_inv_soln).plv_contrast_idx = cx_idx2;
                plv_contrasts = plv_contrasts(good_idx,:);
                
            end
            h.inv_soln(h.current_inv_soln).plv_contrasts = plv_contrasts;
            h.listbox_plv_contrasts.String =  num2str(h.inv_soln(h.current_inv_soln).plv_contrasts);
            h.current_3D_plv_contrasts_listbox_order = 1:size(h.inv_soln(h.current_inv_soln).plv_contrasts,1);
            
            %% creating true contrasts
            chan_contrasts = nchoose2(1:length(h.sim_data.cfg.source.vx_idx));
            h.sim_data.cfg.source.plv_contrast_idx = chan_contrasts;
            h.sim_data.cfg.source.plv_contrasts = h.sim_data.cfg.source.vx_idx(chan_contrasts);
            h.listbox_true_plv_contrasts.String =  num2str(h.sim_data.cfg.source.plv_contrasts);
            update_listbox_plv_contrasts(); 
        case 'No'
    end
end

% if exist('hm','var')
%     delete(hm)
% end

