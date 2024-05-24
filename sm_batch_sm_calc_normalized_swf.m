function [norm_swf, norm_true_swf] = sm_batch_sm_calc_normalized_swf(varargin)
h = varargin{1};
voxel_idx = varargin{end};


%% Normalize True swf
%% Plot TRUE source waves
if nargout>=2
    norm_true_swf=[];
    for v=1:length(h.monte_params.cfg.source.vx_idx)
        data = squeeze(h.sim_data.sig_final(:,v,:)); data = bsxfun(@minus,data,nanmean(data(h.sim_data.cfg.study.base_samps,:,:)))*h.sim_data.cfg.source.vx_amp(v);
        avg_data(:,v) = squeeze(nanmean(data,2));
    end
    if h.radio_normalize_swf.Value == 1
        %         mu = repmat(nanmean(avg_data(h.sim_data.cfg.study.base_samps,v)),size(avg_data,1),1);
        %         sigma = repmat(nanstd(avg_data(h.sim_data.cfg.study.base_samps,v)),size(avg_data,1),1);
        %
        %% New normalization procedure --> Normalize to +/1 max(abs) within act_samps
        act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
        max_p = nanmax(nanmax(abs(avg_data(act_samps,:))));
        avg_data = 100*(avg_data ./ max_p);
    end
    
    for v=1:length(h.monte_params.cfg.source.vx_idx)
        if h.radio_normalize_swf.Value == 1 % plot normalized waves
            %             norm_true_swf(:,v) = (avg_data(:,v)-mu) ./ sigma;
            norm_true_swf(:,v) = avg_data(:,v);
        else
            norm_true_swf(:,v) = avg_data(:,v);
        end
        
        switch h.inv_soln(h.current_inv_soln).Type
            case {'SPA' 'SIA' 'MIA' 'LCMV (FT)' 'sLORETA (FT)'  'dics (FT)' 'pcc (FT)' 'SAM (FT)' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}     % scalar
            case {'eLORETA (FT)' 'MNE (FT)' 'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}   % vector
                switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
                    case {'RMS' 'rms'}
                        norm_true_swf(:,v) = abs(norm_true_swf(:,v));     % taking abs to match calculated swf
                    case {'Max' 'Max'}    % vector dipole orientation with maximum swf power between active and control interval
                        norm_true_swf(:,v) = norm_true_swf(:,v);     % taking abs to match calculated swf
                    case 'avg.pow'
                        norm_true_swf(:,v) = abs(norm_true_swf(:,v));     % taking abs to match calculated swf
                end
        end
        norm_true_swf = bsxfun(@minus,norm_true_swf,nanmean(norm_true_swf(h.monte_params.cfg.study.base_samps,:)));
    end
end


%% Normalize InvSoln swf at voxel_idx
if ~isempty(voxel_idx)
    p_idx = ~isnan(voxel_idx);
    norm_swf = nan(size(h.sim_data.sens_final,1),length(p_idx));
    act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
    ctrl_samps = h.inv_soln(h.current_inv_soln).params.ctrl_samps;
    
    switch h.inv_soln(h.current_inv_soln).Type
        case {'SPA' 'LCMV (FT)' 'sLORETA (FT)'  'dics (FT)' 'pcc (FT)' 'SAM (FT)' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}    % BRANE Lab beamformers
            norm_swf(:,p_idx) = [h.inv_soln(h.current_inv_soln).soln.wts(:,voxel_idx(p_idx))'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';
        case {'SIA' 'MIA'}    % BRANE Lab beamformers
            norm_swf(:,p_idx) = [h.inv_soln(h.current_inv_soln).soln.residual_wts(:,voxel_idx(p_idx))'*squeeze(nanmean(h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:),3))']';
        case {'eLORETA (FT)' 'MNE (FT)'}    % Field Trips Vector inverse solutions
            clear swf swf_pwr
            for ox=1:3
                swf(:,ox,:)=squeeze(nanmean(h.sim_data.sens_final,3))*squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,voxel_idx(p_idx),ox));
                swf_pwr(ox,:)=rms(swf(act_samps,ox,:),1)./rms(swf(ctrl_samps,ox,:),1);
            end
            %        norm_swf=[];
            switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
                case 'RMS'
                    norm_swf(:,p_idx) = squeeze(rms(swf,2));
                case 'Max'    % vector dipole orientation with maximum swf power between active and control interval
                    [~,max_ori]=max(swf_pwr);   % maximum orientation
                    for v=1:size(swf,3); norm_swf(:,v) = squeeze(swf(:,max_ori(v),v)); end
                case 'avg.pow'
                    switch h.inv_soln(h.current_inv_soln).Type
                        case 'eLORETA (FT)'
                            norm_swf(:,find(p_idx==1)) = squeeze(nanmean(abs(swf),2));
                        case 'MNE (FT)'
                            switch h.inv_soln(h.current_inv_soln).headmodel_type
                                case 'Whole Brain'
                                    norm_swf(:,find(p_idx==1)) = h.inv_soln(h.current_inv_soln).soln.avg.pow(h.inv_soln(h.current_inv_soln).params.inside_idx(find(p_idx==1)),:)';
                                case 'Cortical Surface'
                                    norm_swf(:,find(p_idx==1)) = h.inv_soln(h.current_inv_soln).soln.avg.pow(voxel_idx,:)';
                            end
                    end
                    
            end
        case {'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}    % Brainstorm's vector inverse solutions
            swf = h.inv_soln(h.current_inv_soln).soln.ImagingKernel * squeeze(nanmean(h.sim_data.sens_final,3))';
            switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
                case 'RMS'
                    iVertSource = 1:size(swf,1);
                    idx1 = 1:3:length(iVertSource); idx2 = 2:3:length(iVertSource); idx3 = 3:3:length(iVertSource);
                    swf_rms = squeeze(sqrt(swf(idx1,:).^2 + swf(idx2,:).^2 + swf(idx3,:).^2));
                    norm_swf(:,p_idx) = swf_rms(voxel_idx(p_idx),:)';
                case 'Max'    % selecting waveform for maximum ori "max_ori" vector already calculatd for the inv_soln
                    dims = size(swf);
                    swf = permute(reshape(swf,[3 dims(1)/3 dims(2)]),[3 2 1]); % [3dipole_vecotrs x grid locs]
                    vx_idx = voxel_idx(p_idx);
                    swf2 = [];
                    for v=1:length(vx_idx); swf2(:,v) = swf(:,vx_idx(v),h.inv_soln(h.current_inv_soln).soln.max_ori(vx_idx(v))); end
                    %                 iVertSource = 1:size(swf,1);
                    %                 idx1 = 1:3:length(iVertSource); idx2 = 2:3:length(iVertSource); idx3 = 3:3:length(iVertSource);
                    %                 [swf_max, max_ori] = max(abs(squeeze(cat(3, swf(idx1,:), swf(idx2,:), swf(idx3,:)))), [], 3);
                    %                 norm_swf(:,p_idx) = swf_max(voxel_idx(p_idx),:)';
                    norm_swf(:,p_idx) = swf2;
            end
    end
    
    
    % baselining
    norm_swf = bsxfun(@minus,norm_swf,nanmean(norm_swf(h.monte_params.cfg.study.base_samps,:)));
    
    
    %% New normalization procedure --> Normalize to +/1 max(abs) within act_samps
    if isfield(h.inv_soln(h.current_inv_soln),'classifier_metrics')
        if isempty(h.inv_soln(h.current_inv_soln).classifier_metrics) % Hits will be empty when called in from thread starting with run_soruce_modeling.m
            hit_idx = [];
        else
            hit_idx = h.inv_soln(h.current_inv_soln).classifier_metrics.Hits;
        end
    else
         hit_idx = [];
    end
    % hit_idx mathing with voxel_idx
    if any(ismember(voxel_idx,hit_idx))    % no hits within voxel_idx - defaulting to max of swf for all voxel_idx
        max_p = nanmax(nanmax(abs(norm_swf(act_samps,ismember(voxel_idx,hit_idx)))));   % getting max_p for Hits
    else
        max_p = nanmax(nanmax(abs(norm_swf(act_samps,:))));
    end
    norm_swf = 100* (norm_swf ./ max_p);
end



