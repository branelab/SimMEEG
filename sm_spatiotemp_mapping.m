function sm_spatiotemp_mapping(varargin)
global h


if ~isfield(h.inv_soln(h.current_inv_soln).soln.P,'img_org')     % saving original image from initial invSoln
    h.inv_soln(h.current_inv_soln).soln.P.img_org = h.inv_soln(h.current_inv_soln).soln.P.img;
end

% s=h.sim_data.cfg.study.act_samps;
% s_ctrl=h.sim_data.cfg.study.ctrl_samps;
s=h.cfg.study.bl_bmf.act_samps;
s_ctrl=h.cfg.study.bl_bmf.ctrl_samps;

ssp = nanmean(h.sim_data.sig_final(s,:,:),3);
ssp_ctrl = nanmean(h.sim_data.sig_final(s_ctrl,:,:),3);
norm_ssp = bsxfun(@rdivide, ssp, nanstd(ssp_ctrl,[],1));
norm_ssp = bsxfun(@minus, norm_ssp, nanmean(norm_ssp,1));   % baselining across act_int
vox_pos = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos; % h.anatomy.leadfield.voxel_pos;
st_idx =[]; % spatiotemporal indices of found peak sources
search_thresh = str2num(h.edit_inv_peak_spread.String);  % search distance to find peaks
dist_thresh = str2num(h.edit_inv_dist_search_thresh.String);   % search dist around true sources to find hits

%% Spatiotemporal search for finding peak locations for single-source inv solns

switch h.inv_soln(h.current_inv_soln).Type
    case {'SPA' 'LCMV (FT)' 'sLORETA (FT)' 'dics (FT)' 'pcc (FT)' 'SAM (FT)' 'SIA' 'MIA' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}    % BRANE Lab beamformers
        swf = h.inv_soln(h.current_inv_soln).soln.wts' * squeeze(nanmean(h.sim_data.sens_final(s,h.anatomy.sens.good_sensors,:),3))';
        swf_ctrl = h.inv_soln(h.current_inv_soln).soln.wts' * squeeze(nanmean(h.sim_data.sens_final(s_ctrl,h.anatomy.sens.good_sensors,:),3))';
        if h.radio_normalize_swf.Value == 1; swf_base = abs(h.inv_soln(h.current_inv_soln).soln.wts' * squeeze(nanmean(h.sim_data.sens_final(h.sim_data.cfg.study.base_samps,h.anatomy.sens.good_sensors,:),3))'); end
    case {'eLORETA (FT)' }    % Field Trips inverse solutions
        swf=[]; swf_ctrl=[];
        for ox = 1:size(h.inv_soln(h.current_inv_soln).soln.wts,3)
            swf(:,ox,:)=squeeze(nanmean(h.sim_data.sens_final(s,h.anatomy.sens.good_sensors,:),3))*squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,:,ox));
            swf_ctrl(:,ox,:)=squeeze(nanmean(h.sim_data.sens_final(s,h.anatomy.sens.good_sensors,:),3))*squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,:,ox));
        end
        switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
            case 'RMS'
                swf = squeeze(rms(swf,2))'; % taking RMS of waveforms across orientations
                swf_ctrl = squeeze(rms(swf_ctrl,2))'; % taking RMS of waveforms across orientations
            case 'Max'    % vector dipole orientation with maximum swf power between active and control interval
                swf = squeeze(max(swf,[],2))'; % taking RMS of waveforms across orientations
                swf_ctrl = squeeze(max(swf_ctrl,[],2))'; % taking RMS of waveforms across orientations
            case 'avg.pow'
                swf = squeeze(rms(swf,2))'; % taking RMS of waveforms across orientations
                swf_ctrl = squeeze(rms(swf_ctrl,2))'; % taking RMS of waveforms across orientations
        end
    case {'MNE (FT)' 'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}    % vector inverse solutions
        swf=[]; swf_ctrl=[];
        for ox = 1:size(h.inv_soln(h.current_inv_soln).soln.wts,3)
            swf(:,ox,:)=squeeze(nanmean(h.sim_data.sens_final(s,h.anatomy.sens.good_sensors,:),3))*squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,:,ox));
            swf_ctrl(:,ox,:)=squeeze(nanmean(h.sim_data.sens_final(s,h.anatomy.sens.good_sensors,:),3))*squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,:,ox));
        end
        
        switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
            case 'RMS'
                swf = squeeze(rms(swf,2))'; % taking RMS of waveforms across orientations
                swf_ctrl = squeeze(rms(swf_ctrl,2))'; % taking RMS of waveforms across orientations
            case 'Max'    % vector dipole orientation with maximum swf power between active and control interval
                swf = squeeze(max(swf,[],2))'; % taking RMS of waveforms across orientations
                swf_ctrl = squeeze(max(swf_ctrl,[],2))'; % taking RMS of waveforms across orientations
            case 'avg.pow'
                swf = squeeze(h.inv_soln(h.current_inv_soln).soln.avg.pow(:,s)); % taking RMS of waveforms across orientations
                swf_ctrl = squeeze(h.inv_soln(h.current_inv_soln).soln.avg.pow(:,s)); % taking RMS of waveforms across orientations
        end
        
end

%% find peaks in active interval
fprintf('Running Spatiotemporal mapping between %.f - %.f ms\n',h.sim_data.cfg.study.lat_sim(s([1 end]))*1000);
[~, xs] = findpeaks(max(abs(swf)));
h.inv_soln(h.current_inv_soln).soln.time_peak_samps = h.inv_soln(h.current_inv_soln).params.act_samps(xs);
h.inv_soln(h.current_inv_soln).soln.time_peak_lat = h.cfg.study.lat_sim(h.inv_soln(h.current_inv_soln).soln.time_peak_samps);
fprintf('Peak Latencies found (ms): '); fprintf('%.f ', h.inv_soln(h.current_inv_soln).soln.time_peak_lat*1000); fprintf('\n');

%% find all local maxima in abs(swf) for each peak latency found
p_idx = [];
for ss=1:length(xs)
    %         fprintf('%.f ',h.inv_soln(h.current_inv_soln).soln.time_peak_lat(ss)*1000)
    voxel_vals=[vox_pos, abs(swf(:,xs(ss)))];
    null_dist = sort(reshape(abs(swf_ctrl),[numel(swf_ctrl), 1])); thresh_val = null_dist(ceil(length(null_dist)*.95)); % nanmedian(abs(swf(:,ss)));  % finding median value to speed up;
    [peak_voxel,pidx]=BRANELab_find_peak_voxel_thresh(voxel_vals,thresh_val,search_thresh);   % searches for all peaks
    p_idx = unique([p_idx; pidx]);
end

img = nanmax(abs(swf(:,xs)),[],2);

h.current_3D_peak_voxels = [vox_pos(p_idx,:), img(p_idx), p_idx]; h.current_3D_peak_idx = p_idx;
h.inv_soln(h.current_inv_soln).peak_idx = p_idx;
h.inv_soln(h.current_inv_soln).peak_voxels = h.current_3D_peak_voxels;
h.inv_soln(h.current_inv_soln).soln.P.img = img;    % overwriting with combined spatiotemporal maps using rms(swf)

h.inv_soln(h.current_inv_soln).soln.plot_min_max = [min(img) max(img)];



