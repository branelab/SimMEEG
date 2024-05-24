function [norm_swf, norm_true_swf] = sm_calc_normalized_swf(voxel_idx)
global h




       %% New normalization procedure --> Normalize to +/1 max(abs) within act_samps
        act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
        max_p = nanmax(nanmax(abs(avg_data(act_samps,:))));
        avg_data = 100*(avg_data ./ max_p);

%% Normalize True swf
norm_true_swf=[];
for v=1:length(h.cfg.source.vx_idx)
    data = squeeze(h.sim_data.sig_final(:,v,:)); data = bsxfun(@minus,data,nanmean(data(h.sim_data.cfg.study.base_samps,:,:)))*h.sim_data.cfg.source.vx_amp(v);
    avg_data(:,v) = squeeze(nanmean(data,2));
    mu = repmat(nanmean(avg_data(h.sim_data.cfg.study.base_samps,v)),size(avg_data,1),1);
    sigma = repmat(nanstd(avg_data(h.sim_data.cfg.study.base_samps,v)),size(avg_data,1),1);
    if h.radio_normalize_swf.Value == 1 % plot normalized waves
        norm_true_swf(:,v) = (avg_data(:,v)-mu) ./ sigma;
    else
        norm_true_swf(:,v) = avg_data(:,v);
    end
    
    switch h.inv_soln(h.current_inv_soln).Type
        case {'SPA' 'SIA' 'MIA' 'LCMV'  'dics' 'pcc' 'sLORETA' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}     % scalar
        case {'eLORETA' 'MNE'}   % vector
            switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
                case 'RMS'
                    norm_true_swf(:,v) = abs(norm_true_swf(:,v));     % taking abs to match calculated swf
                case 'Max'    % vector dipole orientation with maximum swf power between active and control interval
                    norm_true_swf(:,v) = norm_true_swf(:,v);     % taking abs to match calculated swf
                case 'avg.pow'
                    norm_true_swf(:,v) = abs(norm_true_swf(:,v));     % taking abs to match calculated swf
            end
    end
    norm_true_swf = bsxfun(@minus,norm_true_swf,nanmean(norm_true_swf(h.cfg.study.base_samps,:)));
    
end

%% Normalize InvSoln swf
switch h.inv_soln(h.current_inv_soln).Type
    case {'SPA' 'LCMV (FT)' 'sLORETA (FT)'  'dics (FT)' 'pcc (FT)' 'SAM (FT)' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}    % scalar inverse solutions
        current_peak_swf = [h.inv_soln(h.current_inv_soln).soln.wts(:,voxel_idx)'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';
    case {'SIA' 'MIA'}                                  % multi-source scalar inverse solutions
        current_peak_swf = [h.inv_soln(h.current_inv_soln).soln.residual_wts(:,voxel_idx)'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';
    case {'eLORETA (FT)' 'MNE (FT)'}                    % Field Trip's vector inverse solutions
        act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
        ctrl_samps = h.inv_soln(h.current_inv_soln).params.ctrl_samps;
        clear swf swf_pwr
        for ox=1:3
            swf(:,ox,:)=squeeze(nanmean(h.sim_data.sens_final,3))*squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,voxel_idx,ox));
            swf_pwr(ox,:)=rms(swf(act_samps,ox,:),1)./rms(swf(ctrl_samps,ox,:),1);
        end
        current_peak_swf=[];
        switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
            case 'RMS'
                current_peak_swf = squeeze(rms(swf,2));
            case 'Max'    % vector dipole orientation with maximum swf power between active and control interval
                [~,max_ori]=max(swf_pwr);   % maximum orientation
                for v=1:size(swf,3); current_peak_swf(:,v) = squeeze(swf(:,max_ori(v),v)); end
            case 'avg.pow'
                switch h.inv_soln(h.current_inv_soln).Type
                    case {'eLORETA (FT)' 'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}
                        current_peak_swf = squeeze(nanmean(abs(swf),2));
                    case 'MNE (FT)'
                        switch h.inv_soln(h.current_inv_soln).headmodel_type
                            case 'Whole Brain'
                                current_peak_swf = h.inv_soln(h.current_inv_soln).soln.avg.pow(h.inv_soln(h.current_inv_soln).params.inside_idx(voxel_idx),:)';
                            case 'Cortical Surface'
                                current_peak_swf = h.inv_soln(h.current_inv_soln).soln.avg.pow(voxel_idx,:)';
                        end
                end
                
        end
   case {'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}      % Brainstorm's vector inverse solutions
        swf = h.inv_soln(h.current_inv_soln).soln.ImagingKernel * squeeze(nanmean(h.sim_data.sens_final,3))';
       switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
            case 'RMS'
                iVertSource = 1:size(swf,1);
                idx1 = 1:3:length(iVertSource); idx2 = 2:3:length(iVertSource); idx3 = 3:3:length(iVertSource);
                swf_rms = squeeze(sqrt(swf(idx1,:).^2 + swf(idx2,:).^2 + swf(idx3,:).^2));
                current_peak_swf = swf_rms(voxel_idx,:)';
            case 'Max'    % vector dipole orientation with maximum img power between active and control interval
                iVertSource = 1:size(swf,1);
                idx1 = 1:3:length(iVertSource); idx2 = 2:3:length(iVertSource); idx3 = 3:3:length(iVertSource);
                [swf_max, max_ori] = max(abs(squeeze(cat(3, swf(idx1,:), swf(idx2,:), swf(idx3,:)))), [], 3);
                current_peak_swf = swf_max(voxel_idx,:)';
        end
        
end


% baselining
current_peak_swf = bsxfun(@minus,current_peak_swf,nanmean(current_peak_swf(h.cfg.study.base_samps,:)));

%% normalizing z-transform relative to standard deviation of the baseline
mu_base = repmat(nanmean(current_peak_swf(h.sim_data.cfg.study.base_samps,:)),size(current_peak_swf,1),1);
sigma_base = repmat(nanstd(current_peak_swf(h.sim_data.cfg.study.base_samps,:)),size(current_peak_swf,1),1);

norm_swf = (current_peak_swf-mu_base) ./ sigma_base;    % Z-score
