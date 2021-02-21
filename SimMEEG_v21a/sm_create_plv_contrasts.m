function sm_create_plv_contrasts(varargin)
global h

% h.radio_3D_plot_plv_locs.Value = 0; sm_plot_plv_locs; % turning off leadfield and deleting from 3D MRI plot

if ~isfield(h.inv_soln(h.current_inv_soln),'plv_seed_idx') && ~isfield(h.inv_soln(h.current_inv_soln),'plv_comp_idx')
    hm = warndlg(sprintf('Please select Seed and Comparison Locations\n'),'Select Seed and Comparions Locations'); WinOnTop(hm);
    hm.Position(3:4)=[350 80]; htext = findobj(hm, 'Type', 'Text'); htext.FontSize = 11; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
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
            hm = msgbox(sprintf('Generating PLV contrasts.\n\nThis may take time if large Lead Field Grid selected.\n'),'Generating PLV/PLI contrasts'); WinOnTop(hm);
            hm.Position(3:4)=[350 120]; htext = findobj(hm, 'Type', 'Text'); htext.FontSize = 11; htext.HorizontalAlignment = 'left'; drawnow;
            
            
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
                    [cx_comp_idx,cc_idx] = setdiff(comp_idx,seed_idx(s));
                    cx_idx = cat(2,repmat(seed_idx(s),length(cx_comp_idx),1), cx_comp_idx');
                    % concatenate seed-to-comp contrasts
                    plv_contrasts = cat(1,plv_contrasts,cx_idx);
                    cx_idx2 = cat(1, cx_idx2, [repmat(s,length(cc_idx),1), cc_idx]);    
                end
                % corresponds to indices of the plv_contrasts that are needed when calculating source waveforms (swf) in sm_calc_PLV_PLI.m
                % clmn 1 = seed_idx;  clmn 2 = comp_idx
                h.inv_soln(h.current_inv_soln).plv_contrast_idx = cx_idx2;
            end
            h.inv_soln(h.current_inv_soln).plv_contrasts = plv_contrasts;
        case 'No'
    end
end

if exist('hm','var')
    delete(hm)
end

