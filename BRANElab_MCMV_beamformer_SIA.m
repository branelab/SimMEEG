function [SIA]=BRANElab_MCMV_beamformer_SIA(H,mHref,mref,data,act_samps,ctrl_samps,vx_locs,loc_flag,plot_flag,perc_crit,noise_alpha,max_sources)
%function [wts,ori,P,RNcov,MCMV_idx,mHref,null_thresh]=BRANElab_MCMV_beamformer_SIA(H,mHref,mref,data,act_samps,ctrl_samps,vx_locs,loc_flag,plot_flag,perc_crit)
% This program will find voxels iteratively based on procedures in Moiseev et al., 2011 NeuroImage
%
% INPUT:
%   H =  lead field  [channels x 3 orientations x voxels]
%   mHref = lead field for voxels to place nulls (can use wts.scalar_LeadFields from "BRANElab_LCMV_beamformer_singlesource_v3.m".).
%   mref = indices of the mHref nulls -- to start with single-source serach mHref=[]; mref=[];.
%   data = channel data used to find voxels with maximum signal-to-noise. [samples x channels x trials] .
%   act_samps = samples in active interval used for calculating signal-to-noise ratio.
%   ctrl_samps = samples in control interval used for calculating signal-to-noise ratio.
%   loc_flag =  (0) returns 'pseudoZ' for single-source or 'MPZ' for multisource.
%                  (1) returns 'pseudoZ_er' for single-source or 'MER' for multisource.
%                  (2) returns 'pseudoZ_rer' for single-source or 'rMER' for multisource.
%                  (3) returns 'activity_index' for single-source or 'MAI' for multisource.
%           Note: single-source localizer returned if Href=[]; or single-source localizer returned if Href=[channels x voxels].
%   plot_flag = (0) do not plot found voxel locations (1) plot found voxel locations
%   perc_crit = percent of improved value in max(localizer.img) from previous search to be used as a stopping criterion.
%                   The search stops when percent improved < perc_crit
%
% OUTPUT:
%   wts  = MCMV beamformer weights[channels x 3 x voxels]
%   ori = scalar orientation
%   P.type =     possible localizers - only found sources have values, all other voxels are set to zero -- formulae from Moiseev et al., 2011 .
%           pseudoZ       = (H'*Rinv*H)/(H'*Rinv*N*Rinv*H);
%           pseudoZ_er   = (H'*Rinv*Rbar*Rinv*H)/(H'*Rinv*N*Rinv*H)
%           pseudoZ_rer   = (H'*Rinv*Rbar*Rinv*H)/(H'*Rinv*H)
%           activity_index  = (H'*Ninv*H)/(H'*Rinv*H)
%           MPZ             = trace((H'*Rinv*H)/(H'*Rinv*N*Rinv*H))-n,    where n = number of multiple sources.
%           MER             =  trace((H'*Rinv*Rbar*Rinv*H)/(H'*Rinv*N*Rinv*H))-n,    where n = number of multiple sources.
%           rMER             =  trace((H'*Rinv*Rbar*Rinv*H)/(HF'*Rinv*HF))-n,    where n = number of multiple sources.
%           MAI             =  trace((H'*Ninv*H)/(HF'*Rinv*HF))-n,    where n = number of multiple sources.
%   P.img = image for localizer value [voxels x 1]
%
%   RNcov = R, Rinv, N, Ninv, Rbar, Nbar covariance matrices.
%   MCMV_idx = indices for found multi-sources
%
%
%   written  by A. Herdman Jun 13, 2015
%       updates:
%           Aug 15, 2015 - added MAI and correct bug in "bl_lcmv_scalar.m" - now performs eigen solution for localizer speficied.
%
%% NOTES:
%   Tony - Tested on simulated and real data:
%        - added stoping criteria, works with simulated data
%        - PseudoZ_ER threshold = based only on ctrl interval ==> split ctrl_samps in half to get act_samps and ctrl_samps calculate covariance.

%% beamformer based on following equation
%   W = R^-1 * H * ( H' * R^-1 * H)^-1    , where R = covariance matrix [sensor x sensor] and H = Ledfield [sensor x voxel].


if nargin<12
    max_sources = floor(size(H,1)/5);   % 1 dipole source has 5 degrees of freedom (3 spatial and 2 orientations), thus dividing number of channels/sensors by 5
end


[num_chan,~,num_voxels]=size(H);
% making sure that zeros in H matrix are really zero
H(abs(H)<eps)=eps;

wts_MCMV=[];
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
    if 100*(t*xstd)/max(max(std(data)))>10; fprintf('WARNING! Insufficient Rank in Data or Data is empty\n'); SIA=[]; return; end % stopping infinite loop
    % fprintf('Rank = %.f\n',r2);
    fprintf('Percent white noise added for sufficient rank = %.3f %%\n\n',100*(t*xstd)/max(max(std(data))));
end

%% calc covariance
%[R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(data,act_samps,ctrl_samps);
[R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(data,act_samps,ctrl_samps(xs+1:end));
RNcov.R=R; RNcov.Rinv=Rinv; RNcov.Rbar=Rbar; RNcov.N=N; RNcov.Ninv=Ninv; RNcov.Nbar=Nbar;

% splitting ctrl interval in half to calculate null distribution of noise
xs=round(size(ctrl_samps,2)/2);
[nR,nN,nRbar,nNbar,nRinv,nNinv]=BRANELab_calc_cov(data,ctrl_samps(1:xs),ctrl_samps(xs+1:end));
nRNcov.R=nR; nRNcov.Rinv=nRinv; nRNcov.Rbar=nRbar; nRNcov.N=nN; nRNcov.Ninv=nNinv; nRNcov.Nbar=nNbar;

% % localizer image for ctrl_samps --> null_dist for stopping criteria MCMV searches.
% [wn,~,~,~]=bl_lcmv_scalar(H,[],RNcov,[],loc_flag);
if loc_flag==0; Rx=nR; Nx=nN;  bm_type2='pseudoZ'; % pseudoZ
elseif loc_flag==1; Rx=nRbar; Nx=nN; bm_type2='pseudoZ_er'; % pseudoZ_er
elseif loc_flag==2; Rx=nRbar; Nx=nR; bm_type2='pseudoZ_rer'; % pseudoZ_er
elseif loc_flag==3; Rx=nN; Nx=nR; bm_type2='acitvity_index'; % pseudoZ_er2
end
% parfor v=1:size(wn,2)
%     null_dist(v)=sqrt((wn(:,v)'*Rx*wn(:,v))/(wn(:,v)'*Nx*wn(:,v)));    % this is the null distribution for each voxel's estimate of the single-source pseudoZ function
% %    null_dist(v)=(wn(:,v)'*Rx*wn(:,v))/(wn(:,v)'*Nx*wn(:,v));    % this is the null distribution for each voxel's estimate of the single-source pseudoZ function
% end
% % if length(ctrl_samps)<=length(act_samps)*2;
% %     null_dist=sort(abs(reshape(null_dist,[numel(null_dist) 1]))/sqrt(length(ctrl_samps)/(length(act_samps)*2))); % because halfing the samples in ctrl interval the threshold should actually be sqrt(2).
% % else
% %     null_dist=sort(abs(reshape(null_dist,[numel(null_dist) 1])));
% % end
% %null_dist=sort(abs(reshape(null_dist,[numel(null_dist) 1]))/sqrt(length(ctrl_samps)/(length(act_samps)*2))); % because halfing the samples in ctrl interval the threshold should actually be sqrt(2).
% null_dist=sort(abs(reshape(null_dist,[numel(null_dist) 1])));
% alpha_level=.99;
% %null_thresh=2*null_dist(ceil( (alpha_level)*size(null_dist,1)));  % because halfing the samples in ctrl interval the threshold should actually be sqrt(2).
% null_thresh=null_dist(ceil( (alpha_level)*size(null_dist,1)));  % because halfing the samples in ctrl interval the threshold should actually be sqrt(2).
% %keyboard;

% First calculation of noise - like SPA
H_idx=1:size(H,3); [~,~,~,nP]=bl_lcmv_scalar_v4(H,H_idx,[],nRNcov,[],loc_flag);
null_thresh=max(nP.img);

mm=1;  %num_loops=30; % number of allowable iterations to find source.
zz=1; sn=1; nd=[];
while mm~=0
    % breaking loop after finding maximum # sources
    if length(mref)>=max_sources;  mm=0; fprintf('Found maximum number of sources set by user: %.f\n',length(mref)); break; end % breaking loop because maximum number of sources found as set by the user
 
    % loop through "bl_lcmv_scalar" fucntion until maximal signal across all peak voxels are maximum???.
    H_idx=1:size(H,3); [~,Hscalar2,ori,P2]=bl_lcmv_scalar_v4(H,H_idx,mHref,RNcov,mref,loc_flag);
    
    % updating of noise to include mHref
%     if ~isempty(mref)
%         % iteratively updating noise to include Href
%         H_idx=1: size(H,3); [~,~,~,nP]=bl_lcmv_scalar_v4(H,H_idx,mHref,nRNcov,mref,loc_flag);
%         null_thresh=max(abs(nP.img));
%         nd=sort(nP.img); nd=nd(nd~=0); null_thresh=nd(ceil(length(nd)*noise_alpha));
%     end
    
    bm_type=P2.type;
    img=abs(P2.img);
    voxel_vals=[vx_locs, img];
    %   nd2=sort(P2.img); nd2=nd2(nd2~=0); thresh_val=nd2(ceil(length(nd2)*noise_alpha)); % used to find next peak
    thresh_val=.99*max(img);
    thresh_limit=15;    % 15 mm
    [peak_voxel,p_idx]=BRANELab_find_peak_voxel_thresh(voxel_vals,thresh_val,thresh_limit);
    peak_voxels=flipud(sortrows([peak_voxel p_idx],4));
    mref=[mref peak_voxels(1,5)];
    mHref=[mHref Hscalar2(:,peak_voxels(1,5))];     % nullying best found ori
    
    %w=wts2(:,peak_voxels(1,5));
    w2=(Rinv*mHref)/(mHref'*Rinv*mHref); % adjusting weights to include all found sources - this should speed things up
    w=w2(:,end);
    
    %     % monitoring pseudoZ functions for stopping criteria
    %     if loc_flag==0; z(sn)=sqrt((w'*R*w)/(w'*N*w));  bm_type2='pseudoZ';% pseudoZ
    %     elseif loc_flag==1; z(sn)=sqrt((w'*Rbar*w)/(w'*N*w));  bm_type2='pseudoZ_ER'; % pseudoZ_ER
    %     elseif loc_flag==2; z(sn)=sqrt((w'*Rbar*w)/(w'*R*w));  bm_type2='pseudoZ_rER';% pseudoZ_rER
    %     elseif loc_flag==3; z(sn)=sqrt((w'*N*w)/(w'*R*w));  bm_type2='activity_index';% activity_index
    %     end
    %
    % monitoring pseudoZ functions for stopping criteria
    if loc_flag==0; z(sn)=(w'*R*w)/(w'*N*w);  bm_type2='pseudoZ';% pseudoZ
    elseif loc_flag==1; z(sn)=(w'*Rbar*w)/(w'*N*w);  bm_type2='pseudoZ_ER'; % pseudoZ_ER
    elseif loc_flag==2; z(sn)=(w'*Rbar*w)/(w'*R*w);  bm_type2='pseudoZ_rER';% pseudoZ_rER
    elseif loc_flag==3; z(sn)=(w'*N*w)/(w'*R*w);  bm_type2='activity_index';% activity_index
    end
    
    %z(sn)=img(peak_voxels(1,5));
    pmax(sn)=abs(img(peak_voxels(1,5)));
    xref(sn)=peak_voxels(1,5);
    fprintf('Sources Found = '); fprintf('%.f ',mref);
    %        fprintf('\nIteration = %.f\tPeak Voxel#:%.f (ori=%.3f %.3f %.3f) \t%s=%.2f\t%s Threshold =%.2f\t SNR = %.2f\n',zz,xref(sn),ori(xref(sn),1),ori(xref(sn),2),ori(xref(sn),3),bm_type,z(sn),bm_type2,null_thresh,z(sn)/null_thresh);
    %        fprintf('\nIteration = %.f\tPeak Voxel#:%.f (ori=%.3f %.3f %.3f) \t%s=%.2f\t%s Threshold =%.2f\t pseudoSNR = %.2f\n',zz,xref(sn),ori(xref(sn),1),ori(xref(sn),2),ori(xref(sn),3),bm_type,z(sn),bm_type,null_thresh,z(sn)/null_thresh);
    
    %%% new stoping criteria - when pmax doesn't improve more than 'perc_crit' then stop .
    if sn>=2; perc_imp=abs(100*(1-(pmax(end-1)/pmax(end)))); else perc_imp=100; end
    
    
    %    if ~isempty(nd); keyboard; end
    
    if plot_flag==1;
        %  bl_plot_lcmv_locs(peak_voxels(1,5),vx_locs,abs(img),hot(255),[min(abs(img)) max(abs(img))]); colorbar;
        figure(999); clf; set(gcf,'color',[1 1 1]*.4); %,'position',[556   390   796   293]);
        bl_plot_lcmv_locs(peak_voxels(1,5),vx_locs,abs(img),hot(255),[min(abs(img(img~=0))) max(abs(img))]); colorbar;
        title('Signal Map');
        subplot(1,4,4); cla; plot(pmax); title('Signal Map'); %plot([pmax; z]');
        if sn>1; legend({bm_type2},'location','SouthEast'); end   % if sn>1; legend({'pmax','z'},'location','NorthWest'); end
        if ~isempty(nd)
            figure(1999); clf; set(gcf,'color',[1 1 1]*.4); %,'position',[556   390   796   293]);
            bl_plot_lcmv_locs(mref,vx_locs,abs(nP.img),hot(255),[min(abs(img(img~=0))) max(abs(img))]); colorbar;
            title('Noise Map'); subplot(1,4,4); cla; histfit(nd); %plot([pmax; z]');
        end
        drawnow
    end
    fprintf('\nIteration = %.f\tPeak Voxel#:%.f (ori=%.3f %.3f %.3f) \n%s=%.4f   %s=%.4f   Threshold =%.4f   pseudoSNR = %.4f   Percent improved=%.2f\n',zz,xref(sn),ori(xref(sn),1),ori(xref(sn),2),ori(xref(sn),3),bm_type,pmax(sn),bm_type2,z(sn),null_thresh,z(sn)/null_thresh,perc_imp);
    
    %         %%% new stoping criteria - when pmax doesn't improve more than 'perc_crit' then stop and take max of last 2.
    if sn>=2 && perc_imp<perc_crit; zz=zz+1;
        %          keyboard;
    elseif sn>=2 && perc_imp>perc_crit; kk=0;   % reset if % improvment gets larger than cirterion (i.e., rising again not plateauing)
    end
    sn=sn+1;
    
    %%% stopping criterion #2- if sources still being found after num_loops.
    zz=zz+1;    % if zz>30 then stopping loop
%     if zz>num_loops; mm=0; end
    %%% finding best Lead-field when accounting for all found sources.
    if length(mref)>1;
        [~,mH]=bl_lcmv_scalar_v4(H,mref,[],RNcov,[],loc_flag);
        mHref2=mH(:,mref);
        %wts_MCMV=zeros(size(wts2));
        for v=1:size(mref,2);
            H_idx=mref(v);
            r_idx=setdiff(1:size(mref,2),v); %selecting all other leadfield for nulling
            [~,mH,mori,P3]=bl_lcmv_scalar_v4(H,H_idx,mHref2(:,r_idx),RNcov,mref(r_idx),loc_flag);
            mHscalar=mH(:,H_idx); ori_MCMV(v,:)=mori(H_idx,:);
            mHref(:,v)=mHscalar;
        end
    elseif  length(mref)==1;
        [~,mH,mori,~]=bl_lcmv_scalar_v4(H,mref,[],RNcov,[],loc_flag);
        mHref(:,1)=mH(:,mref); ori_MCMV(1,:)=mori(mref,:);
    elseif  isempty(mref);
        fprintf('WARNING! No sources found!\n');
        mref=[];  mm=0; break;
    end
    
    % finding if there are sources above null_thresh after accoutning for all multisources
    wts=zeros(num_chan,num_voxels);
    wts(:,mref)=(Rinv*mHref)/(mHref'*Rinv*mHref); % Alex's suggested script - returns same as nutMEG script
    for v=1:size(mref,2)
        mpz(v)=(wts(:,mref(v))'*R*wts(:,mref(v)))/(wts(:,mref(v))'*N*wts(:,mref(v))); % MPZ
        mer(v)=(wts(:,mref(v))'*Rbar*wts(:,mref(v)))/(wts(:,mref(v))'*N*wts(:,mref(v))); % MER
        rmer(v)=(wts(:,mref(v))'*Rbar*wts(:,mref(v)))/(wts(:,mref(v))'*R*wts(:,mref(v))); % rMER
        mai(v)=(wts(:,mref(v))'*Ninv*wts(:,mref(v)))/(wts(:,mref(v))'*Rinv*wts(:,mref(v))); % MAI
    end
    
    if loc_flag==0; r_idx=find(mpz>null_thresh); % MPZ
    elseif loc_flag==1; r_idx=find(mer>null_thresh); % MER
    elseif loc_flag==2; r_idx=find(rmer>null_thresh); % rMER
    elseif loc_flag==3; r_idx=find(mai>null_thresh); % MAI
    end
    fprintf('Sources = '); fprintf('%.f ',mref); fprintf('\t');
    fprintf('Multi Sources = '); fprintf('%.f ',mref(r_idx)); fprintf('\n');
    % stopping criteria
    if length(r_idx)<size(mref,2) && ~isempty(r_idx); mm=0; end % stops loop when no more sources can be found above null_thresh after acounting for all multi-sources.
    if max(z)<=null_thresh ; mm=0; end % stops when the pseudoZ_er value goes below 1 (meaning SNR<1);
            

end
fprintf('Stopping loop - No more sources!\n');

%% thresholding and calculate final weights for MCMV
m=0;
while m==0;
    
    if length(mref)>1;
        [~,mH]=bl_lcmv_scalar_v4(H,mref,[],RNcov,[],loc_flag);
        mHref2=mH(:,mref);
        %wts_MCMV=zeros(size(wts2));
        for v=1:size(mref,2);
            H_idx=mref(v);
            r_idx=setdiff(1:size(mref,2),v); %selecting all other leadfield for nulling
            [~,mH,mori,P3]=bl_lcmv_scalar_v4(H,H_idx,mHref2(:,r_idx),RNcov,mref(r_idx),loc_flag);
            mHscalar=mH(:,H_idx); ori_MCMV(v,:)=mori(H_idx,:);
            mHref(:,v)=mHscalar;
        end
    elseif  length(mref)==1;
        [~,mH,mori,~]=bl_lcmv_scalar_v4(H,mref,[],RNcov,[],loc_flag);
        mHref(:,1)=mH(:,mref); ori_MCMV(1,:)=mori(mref,:);
    elseif  isempty(mref);
        fprintf('WARNING! No sources found!\n');
        mref=[];
        MCMV_idx=[];
        P.type='';
        P.img=[];
        ori=[];
        mHref=[];
        mm=99;
    end
    
    
    if ~isempty(mref);
        % recalculating localizer value based on only including voxels not below noise threshold.
        wts(:,mref)=(Rinv*mHref)/(mHref'*Rinv*mHref); % Alex's suggested script - returns same as nutMEG script
        pimg=zeros(1,length(mref));
        for v=1:size(mref,2)
            if loc_flag==0;  pimg(v)=(wts(:,mref(v))'*R*wts(:,mref(v)))/(wts(:,mref(v))'*N*wts(:,mref(v)));   % MPZ
            elseif loc_flag==1; pimg(v)=(wts(:,mref(v))'*Rbar*wts(:,mref(v)))/(wts(:,mref(v))'*N*wts(:,mref(v)));    %MER
            elseif loc_flag==2; pimg(v)=(wts(:,mref(v))'*Rbar*wts(:,mref(v)))/(wts(:,mref(v))'*R*wts(:,mref(v)));   % rMER
            elseif loc_flag==3; pimg(v)=(wts(:,mref(v))'*Ninv*wts(:,mref(v)))/(wts(:,mref(v))'*Rinv*wts(:,mref(v)));    % MAI
            end
        end
        
%         % Final updating of noise to include mHref
%         if ~isempty(mref)
%             H_idx=1: size(H,3); [~,~,~,nP]=bl_lcmv_scalar_v4(H,H_idx,mHref,nRNcov,mref,loc_flag);      % this is conservative.
%             nd=sort(nP.img); nd=nd(nd~=0); null_thresh=nd(ceil(length(nd)*noise_alpha));
%         end
        
        if isempty(find(pimg<=null_thresh, 1)); m=99;   % breaking loop because all voxels are above threshold now
            m=99;
        elseif isempty(find(pimg>null_thresh,1));       % no significant voxels after thresholding
            fprintf('WARNING! No sources found!\n');
            mref=[];pimg=[]; ori_MCMV=[]; MCMV_idx=[]; P=[]; P.type=''; P.img=[]; ori=[]; mHref=[]; wts=[];
            m=99;
        else
            r_idx=find(pimg>null_thresh);
            mref=mref(r_idx);
            ori_MCMV=ori_MCMV(r_idx,:);
            mHref=mHref(:,r_idx);
        end
    else
        fprintf('WARNING! No sources found!\n');
        mref=[];pimg=[]; ori_MCMV=[]; MCMV_idx=[]; P=[]; P.type=''; P.img=[]; ori=[]; mHref=[]; wts=[];
        m=99;
    end
    
end

MCMV_idx=mref;
P.img=zeros(num_voxels,1);

if loc_flag==0; % MPZ
    P.type='MPZ';
    P.img(MCMV_idx)=pimg; %(wts(:,v)'*R*wts(:,v))/(wts(:,v)'*N*wts(:,v));
elseif loc_flag==1;  %MER
    P.type='MER';
    P.img(MCMV_idx)=pimg; %(wts(:,v)'*Rbar*wts(:,v))/(wts(:,v)'*N*wts(:,v));
elseif loc_flag==2; % rMER
    P.type='rMER';
    P.img(MCMV_idx)=pimg; %(wts(:,v)'*Rbar*wts(:,v))/(wts(:,v)'*N*wts(:,v));
elseif loc_flag==3; % MAI
    P.type='MAI';
    P.img(MCMV_idx)=pimg; %(wts(:,v)'*Rbar*wts(:,v))/(wts(:,v)'*N*wts(:,v));
end
    
P.img(P.img<=null_thresh)=0; % double checking to remove all voxel values below null_thresh
ori=ori_MCMV; 

fprintf('Final Sources =\t Value'); fprintf('\n');
for v=1:length(MCMV_idx); fprintf('%.f \t %.3f\n',MCMV_idx(v),P.img(MCMV_idx(v))); end
fprintf('Noise \t %.3f\n',null_thresh);

%% generate noise img by running MCMV beamformer wts for all voxels but the MCMV_idx are the ref voxels.
[~,nHscalar2,nori,nP]=bl_lcmv_scalar_v4(H,H_idx,mHref,RNcov,MCMV_idx,loc_flag);

%% Residual maps
cov_flag = 1;   % using Rinv for weights
[residual_wts,residual_Hscalar,residual_ori]=calc_nulled_wts(H,mHref,RNcov,mref,loc_flag,ori_MCMV,cov_flag);
[pimg3,npimg3,~]=calc_localizer(residual_wts,mref,RNcov,nRNcov,loc_flag);
nd=sort(abs(npimg3(~isnan(npimg3)))); nd=nd(nd~=0);  residual_thresh=nd(ceil(length(nd)*(1-noise_alpha)));

%% % Final image and weights
cov_flag = 2;   % using Ninv for weights
[nulled_wts,nulled_Hscalar,nulled_ori]=calc_nulled_wts(H,mHref,RNcov,mref,loc_flag,ori_MCMV,cov_flag);
[pimg2,npimg2,~]=calc_localizer(nulled_wts,mref,RNcov,nRNcov,loc_flag);
nd=sort(abs(npimg2(~isnan(npimg2)))); nd=nd(nd~=0);  nulled_thresh=nd(ceil(length(nd)*(1-noise_alpha)));


%% OUTPUTS
SIA.wts=wts;
SIA.P=P;
SIA.RNcov=RNcov;
SIA.MCMV_idx=MCMV_idx;
SIA.mHref=mHref;
SIA.null_thresh=null_thresh;
SIA.ori=nan(size(SIA.wts,2),3); 
SIA.ori(SIA.MCMV_idx,:)=ori;
SIA.P_log.img=20*log(P.img/SIA.null_thresh);
SIA.P_SNR.img=SIA.P.img/SIA.null_thresh;    % MER SNR = MER(act_int) / pseudoZ(ctrl_int)
SIA.P_noise.img=nP.img;
SIA.P_noise.ori=nori;
SIA.P_noise.Hscalar=nHscalar2;

% MCMV + Residual maps --> using Rinv to create maps
SIA.residual_wts=residual_wts;
SIA.residual_Hscalar=residual_Hscalar;
SIA.P.residual_img=pimg3;
SIA.P_ctrl.residual_img=npimg3;
SIA.residual_thresh=residual_thresh;
SIA.residual_ori = residual_ori;

% nulled --> using Ninv to create maps 
SIA.nulled_wts=nulled_wts;
SIA.nulled_Hscalar=nulled_Hscalar;
SIA.P.nulled_img=pimg2;
SIA.P_ctrl.nulled_img=npimg2;
SIA.nulled_thresh=nulled_thresh;
SIA.nulled_ori = nulled_ori;


function [nulled_wts,nulled_Hscalar,vx_ori]=calc_nulled_wts(H,mHref,RNcov,mref,loc_flag,mref_ori,cov_flag)
% H_idx=1:size(H,3); [~,~,vx_ori,~]=bl_lcmv_scalar_v4(H,H_idx,[],RNcov,[],loc_flag);
H_idx=1:size(H,3); [~,~,vx_ori,~]=bl_lcmv_scalar_v4(H,H_idx,mHref,RNcov,mref,loc_flag);
if ~isempty(mref) && ~isempty(mref_ori)
    vx_ori(mref,:)=mref_ori; % replacing SPA_ori with SIA_ori
end
% nullying all voxels using mHref for the found mref locations
parfor v=1:size(H,3)
    h = H(:,:,v)*vx_ori(v,:)';     % lead field for dipole with best orientation (yields highest signal).
    h = h/norm(h);              % normalizes the leadfield across channels for each voxel (Needed to get proper weights)
    
    if ismember(v,mref)
        vnx=find(mref==v);
        nx=1:length(mref);
        null_vx=setdiff(nx,vnx);
    else
        null_vx=1:length(mref);
    end
    %       null_vx=setdiff(mref,v);
    %    if isempty(find(mref==v, 1))
    %         null_vx=1:length(mref); % using all mref for nullying
    %     else
    %         vx=find(mref==v,1); null_vx=setdiff(1:length(mref),vx); % removing the current voxel in mref when nullying
    %     end
    HF=[h mHref(:,null_vx)]; % inlcuding ref voxels in "nulled_idx" for nullying
    if cov_flag==1      % use Rinv for weights
        v_wts=(RNcov.Rinv*HF)/(HF'*RNcov.Rinv*HF); % Alex's suggested script - returns same as nutMEG script
    elseif cov_flag==2  % use Rinv for weights
        v_wts=(RNcov.Ninv*HF)/(HF'*RNcov.Ninv*HF); % using Ninv Noise Covarianace to get better supression of noise because it doesn't include possible active source as in Rinv
    end
    nulled_wts(:,v)=v_wts(:,1);
    nulled_Hscalar(:,v)=HF(:,1);
end


function [pimg,npimg,bmf_type]=calc_localizer(wts,mref,RNcov,nRNcov,loc_flag)

%% Final image calculation after final nullying
if loc_flag==0
    bmf_type='MPZ';% pseudoZ for active interval
    pimg=(wts'*RNcov.R*wts)./(wts'*RNcov.N*wts);  pimg=diag(pimg);
    npimg=(wts'*nRNcov.R*wts)./(wts'*nRNcov.N*wts);   npimg=diag(npimg);
elseif loc_flag==1
    bmf_type='MER'; % pseudoZ_ER for active interval
    if ~isempty(mref)
        N2=RNcov.N/length(mref); nN2=nRNcov.N/length(mref);
    else
        N2=RNcov.N; nN2=nRNcov.N;
    end
    pimg=(wts'*RNcov.Rbar*wts)./(wts'*N2*wts);  pimg=diag(pimg); % all ref_idx
    npimg=(wts'*nRNcov.Rbar*wts)./(wts'*nN2*wts);  npimg=diag(npimg); % all ref_idx
elseif loc_flag==2
    bmf_type='rMER';% pseudoZ_rER for active interval
    pimg=(wts'*RNcov.Rbar*wts)./(wts'*RNcov.R*wts);    pimg=diag(pimg);  % all ref_idx
    npimg=(wts'*nRNcov.Rbar*wts)./(wts'*nRNcov.R*wts);   npimg=diag(npimg); % all ref_idx
elseif loc_flag==3
    bmf_type='MAI';% activity_index for active interval
    pimg=(wts'*RNcov.N*wts)./(wts'*RNcov.R*wts);  pimg=diag(pimg); % all ref_idx
    npimg=(wts'*nRNcov.N*wts)./(wts'*nRNcov.R*wts);  npimg=diag(npimg); % all ref_idx
end

