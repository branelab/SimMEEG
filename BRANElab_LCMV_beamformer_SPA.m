function [SPA]=BRANElab_LCMV_beamformer_SPA(H,data,act_samps,ctrl_samps,loc_flag,noise_alpha)
%function [wts,Hscalar,ori,P,RNcov,null_thresh]=BRANElab_LCMV_beamformer_SPA(H,data,act_samps,ctrl_samps,loc_flag)
% This program will find voxels with maximal amplitudes from scalar beamformer
%
% INPUT:
%   H =  lead field  [channels x 3 orientations x voxels]
%   data = channel data used to find voxels with maximum signal-to-noise. [samples x channels x trials] .
%   act_samps = samples in active interval used for calculating signal-to-noise ratio.
%   ctrl_samps = samples in control interval used for calculating signal-to-noise ratio.
%   loc_flag =  (0) returns 'pseudoZ' for single-source or 'MPZ' for multisource.
%                  (1) returns 'pseudoZ_er' for single-source or 'MER' for multisource.
%                  (2) returns 'pseudoZ_rer' for single-source or 'rMER' for multisource.
%                  (3) returns 'activity_index' for single-source or 'MAI' for multisource.
%           Note: single-source localizer returned if Href=[]; or single-source localizer returned if Href=[channels x voxels].  
%
% OUTPUT:
%   wts = source weights [channels x voxels];
%   Hscalar = scalar lead-field matrix based on single- or mult-source
%   ori = orientation of voxels [voxels x (X,Y,Z)];
%   P.type = possible localizers -  formulae from Moiseev et al., 2011 .
%           pseudoZ       = (H'*Rinv*H)/(H'*Rinv*N*Rinv*H);
%           pseudoZ_er   = (H'*Rinv*Rbar*Rinv*H)/(H'*Rinv*N*Rinv*H)
%           pseudoZ_rer   = (H'*Rinv*Rbar*Rinv*H)/(H'*Rinv*H)
%           activity_index  = (H'*Ninv*H)/(H'*Rinv*H)
%           MPZ             = trace((H'*Rinv*H)/(H'*Rinv*N*Rinv*H))-n,    where n = number of multiple sources.
%           MER             =  trace((H'*Rinv*Rbar*Rinv*H)/(H'*Rinv*N*Rinv*H))-n,    where n = number of multiple sources.
%           rMER             =  trace((H'*Rinv*Rbar*Rinv*H)/(HF'*Rinv*HF))-n,    where n = number of multiple sources.
%           MAI             =  trace((H'*Ninv*H)/(HF'*Rinv*HF))-n,    where n = number of multiple sources.
%   P.img = image for localizer values for the active interval [voxels x 1]
%   P_ctrl.img = image for localizer values for the control interval [voxels x 1]
%   P_log.img = normalized power image relative to ctrl_int: P_log.img=20*log(P.img./P_ctrl.img);
%
%   RNcov = R, Rinv, N, Ninv, Rbar, Nbar covariance matrices.
%   null_thresh = localizer threshold based on LCMV of ctrl interval by splitting ctrl samps into half active and half control.
%
%   written  by A. Herdman Jun 13, 2015
%

%% beamformer based on following equation
%   W = R^-1 * H * ( H' * R^-1 * H)^-1    , where R = covariance matrix [sensor x sensor] and H = Ledfield [sensor x voxel].

fprintf('Performing Single-Source beamforming\n');
SPA = [];

%% LCMV
% calculating covarince matrices for signal+noise (Rinv) and noise (Ninv
xs=round(size(ctrl_samps,2)/2);

%% checking rank of data. If it doesn't match the num_chans then white noise will be added until rank sufficient.
R=[]; t=0;
[R,~,~,~,~,~]=BRANELab_calc_cov(data,act_samps,ctrl_samps(xs+1:end));
r2=rank(R);
xstd=min(min(std(data)))/size(data,2);
while r2<size(H,1)
    t=t+1;
    data=data+(xstd*randn(size(data)));
    [R,~,~,~,~,~]=BRANELab_calc_cov(data,act_samps,ctrl_samps(xs+1:end));
    r2=rank(R);
    fprintf('Rank = %.f\tChannel Count=%.f\n',r2,size(data,2));
    if 100*(t*xstd)/max(max(std(data)))>10; fprintf('WARNING! Insufficient Rank in Data or Data is empty\n'); MIA=[]; return; end % stopping infinite loop
    % fprintf('Rank = %.f\n',r2);
    fprintf('Percent white noise added for sufficient rank = %.3f %%\n\n',100*(t*xstd)/max(max(std(data))));
end

%% Data Covariance 
%[R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(data,act_samps,ctrl_samps);
[R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(data,act_samps,ctrl_samps(xs+1:end));
RNcov.R=R; RNcov.Rinv=Rinv; RNcov.Rbar=Rbar; RNcov.N=N; RNcov.Ninv=Ninv; RNcov.Nbar=Nbar; 

H_idx=1:size(H,3); [wts,Hscalar,ori,P]=bl_lcmv_scalar_v5(H,H_idx,[],RNcov,[],loc_flag);

%%  Noise estimate from control samples
% splitting ctrl interval in half to calculate null distribution of noise
%xs=round(size(ctrl_samps,2)/2);
[nR,nN,nRbar,nNbar,nRinv,nNinv]=BRANELab_calc_cov(data,ctrl_samps(1:xs),ctrl_samps(xs+1:end));
nRNcov.R=nR; nRNcov.Rinv=nRinv; nRNcov.Rbar=nRbar; nRNcov.N=nN; nRNcov.Ninv=nNinv; nRNcov.Nbar=nNbar;

% localizer image for ctrl_samps.
if loc_flag==0; Rx=nR; Nx=nN;  bm_type2='pseudoZ'; % pseudoZ
elseif loc_flag==1; Rx=nRbar; Nx=nN; bm_type2='pseudoZ_er'; % pseudoZ_er
elseif loc_flag==2; Rx=nRbar; Nx=nR; bm_type2='pseudoZ_rer'; % pseudoZ_er
elseif loc_flag==3; Rx=nN; Nx=nR; bm_type2='acitvity_index'; % pseudoZ_er2
end

% First calculation of noise - like SPA
H_idx=1:size(H,3); [~,~,~,nP]=bl_lcmv_scalar_v5(H,H_idx,[],nRNcov,[],loc_flag);
nd=sort(nP.img); null_thresh=nd(ceil(length(nd)*noise_alpha));

% OUTPUTS
SPA.wts=wts;
SPA.Hscalar=Hscalar;
SPA.ori=ori;
SPA.P=P; 
SPA.P_ctrl=nP; 
SPA.RNcov=RNcov; 
SPA.null_thresh=null_thresh;
% SPA.P_log.type='Log normalized Map';
% SPA.P_log.img=20*log(P.img./nP.img);
% SPA.P_SNR.img=SPA.P.img./SPA.P_ctrl.img;


