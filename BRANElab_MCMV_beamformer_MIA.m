function [MIA]=BRANElab_MCMV_beamformer_MIA(H,mHref,mref,data,act_samps,ctrl_samps,vx_locs,loc_flag,plot_flag,perc_crit,noise_alpha,text_flag,anatomy,max_sources)
%function [MIA]=BRANElab_MCMV_beamformer_MIA(H,mHref,mref,data,act_samps,ctrl_samps,vx_locs,loc_flag,plot_flag,perc_crit);
% This program will find voxels with maximal amplitudes from scalar beamformer and iteratively null these voxels up to the max num_sources.
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
%   plot_flag = (0) do NOT plot images online (1)  plot images online
%   perc_crit = percent of improved value in max(localizer.img) from previous search to be used as a stopping criterion.
%                   The search stops when percent improved < perc_crit
%   text_flag = (0) print out online final MIA information (1) print out information for each loop
%   anatomy = anatomy stucture
%             .vol.bnd.tri = triangles (aka faces) for volumes generated using Field Trip's code
%             .vol.bnd.pos = node position (aka vertices) for volumes generated using Field Trip's code
%             .leadfield.pos = position of voxels [M x 3]
%             .leadfield.inside = indices for voxels within the brain volume or cortical surfaces.
%             .leadfield.voxel_pos = brain volume voxel locations (inside voxel locations)
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
%        - real data showed a lot of low signal localized sources, most likely becuase noise thresh is not updated.
%           Thus, placing null_thresh within iterative loop to update when each source is found.
%
% Update Sept 28, 2018:
%   Tony
%       Major - found bug in previous versions whereby the pseudoZer Noise (N) was not being divided by
%               the number of ref sources. Previous MIA versions were thus
%               too conservative and had more misses than what whould have
%               occured. Now tested with 9 simulated source with high SNR
%               and it finds all of them.
%       Minor - added new images for online plotting - requires "anatomy" structure
%

%% beamformer based on following equation
%   W = R^-1 * H * ( H' * R^-1 * H)^-1    , where R = covariance matrix [sensor x sensor] and H = Ledfield [sensor x voxel].

fprintf('Performing MCMV MIA beamforming\n');

if nargin<12
    text_flag=0;
elseif nargin<13 && plot_flag==1
        text_flag=0;plot_flag=0;
        fprintf('Not plotting image results online becuase input variable "anatomy" is not provided\n');
elseif nargin<14
    text_flag=0; plot_flag=0;
    max_sources = floor(size(H,1)/5);   % 1 dipole source has 5 degrees of freedom (3 spatial and 2 orientations), thus dividing number of channels/sensors by 5
end





[num_chan,~,num_voxels]=size(H);
% making sure that zeros in H matrix are really zero
H(abs(H)<eps)=eps;

% wts_MCMV=[];
ori_MCMV=[];
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

%% calc covariance
%[R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(data,act_samps,ctrl_samps);
[R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(data,act_samps,ctrl_samps(xs+1:end));
RNcov.R=R; RNcov.Rinv=Rinv; RNcov.Rbar=Rbar; RNcov.N=N; RNcov.Ninv=Ninv; RNcov.Nbar=Nbar;

% splitting ctrl interval in half to calculate null distribution of noise
% xs=round(size(ctrl_samps,2)/2);
[nR,nN,nRbar,nNbar,nRinv,nNinv]=BRANELab_calc_cov(data,ctrl_samps(1:xs),ctrl_samps(xs+1:end));
nRNcov.R=nR; nRNcov.Rinv=nRinv; nRNcov.Rbar=nRbar; nRNcov.N=nN; nRNcov.Ninv=nNinv; nRNcov.Nbar=nNbar;

% % localizer image for ctrl_samps --> null_dist for stopping criteria MCMV searches.
% [wn,~,~,~]=bl_lcmv_scalar(H,[],nRNcov,[],loc_flag);
% if loc_flag==0; Rx=nR; Nx=nN;  bm_type2='pseudoZ'; % pseudoZ
% elseif loc_flag==1; Rx=nRbar; Nx=nN; bm_type2='pseudoZ_er'; % pseudoZ_er
% elseif loc_flag==2; Rx=nRbar; Nx=nR; bm_type2='pseudoZ_rer'; % pseudoZ_er
% elseif loc_flag==3; Rx=nN; Nx=nR; bm_type2='acitvity_index'; % pseudoZ_er2
% end
% [wn,nHscalar,nori,nP]=bl_lcmv_scalar(H,[],nRNcov,[],loc_flag);
% parfor v=1:size(wn,2)
%     if isnan(sqrt((wn(:,v)'*Rx*wn(:,v))/(wn(:,v)'*Nx*wn(:,v))))
%         fprintf('%.f\n',v);
%     else
%         null_dist(v)=sqrt((wn(:,v)'*Rx*wn(:,v))/(wn(:,v)'*Nx*wn(:,v)));    % this is the null distribution for each voxel's estimate of the single-source pseudoZ function
%     end
%     %    null_dist(v)=(wn(:,v)'*Rx*wn(:,v))/(wn(:,v)'*Nx*wn(:,v));    % this is the null distribution for each voxel's estimate of the single-source pseudoZ function
% %     sd(v)=wn(:,v)'*Rx*wn(:,v);
% %     nd(v)=(wn(:,v)'*Nx*wn(:,v));
% end
% null_dist=sort(abs(reshape(null_dist,[numel(null_dist) 1])));
% alpha_level=.99;
% null_thresh=null_dist(ceil( (alpha_level)*size(null_dist,1)));

% H_idx=1: size(H,3); [~,~,~,nP]=bl_lcmv_scalar_v5(H,H_idx,[],nRNcov,[],loc_flag);
% nd=sort(abs(nP.img)); null_thresh=nd(ceil(length(nd)*noise_alpha));

fprintf('Running iterations. Please be patient ...\n');

%% mref=[]; mHref=[]; wts_MCMV=[]; ori_MCMV=[];
mm=1;  num_loops=3; % number of reversal for two sources to be found.
nqloop=0; null_thresh2=[]; z9=[]; nz9=[]; max_iter=20;
while mm~=0
    
    nqloop=nqloop+1;
    kk=0; sn=0; clear pmax pmean npmax npmean zk nzk Zmax nZmax xref z ref_ori Href2 null_thresh; Href=mHref;ref_idx=mref; zz=1; k2=1; ref_ori=[];
    
    %% INNER LOOP
    while kk<num_loops % reversals of same indices found
        sn=sn+1;
        %  for sn=1:100  % seraching through top 2 sources
        
        % loop through "bl_lcmv_scalar" fucntion until maximal signal across all peak voxels are maximum???.
        %         keyboard;
        H_idx=1: size(H,3); [~,Hscalar2,ori,P2]=bl_lcmv_scalar_v5(H,H_idx,Href,RNcov,ref_idx,loc_flag);
        
        % iteratively updating noise to include Href
        H_idx=1: size(H,3); [~,~,~,nP]=bl_lcmv_scalar_v5(H,H_idx,Href,nRNcov,ref_idx,loc_flag);
        %         nd=sort(nP.img(~isnan(nP.img)));  null_thresh(sn)=nd(ceil(length(nd)*noise_alpha));
        nd=sort(abs(nP.img(~isnan(nP.img)))); nd=nd(nd~=0);  null_thresh(sn)=nd(ceil(length(nd)*(1-noise_alpha)));
        
        bm_type=P2.type;
        %         img=abs(P2.img);
        voxel_vals=[vx_locs, abs(P2.img)];
        
        %         thresh_val=.95*max(abs(P2.img));
        %         thresh_limit=15;    % 15 mm
        %         [peak_voxel,p_idx]=BRANELab_find_peak_voxel_thresh(voxel_vals,thresh_val,thresh_limit);
        %         peak_voxels=flipud(sortrows([peak_voxel p_idx],4));
        %         ref_idx=[mref peak_voxels(1,5)];
        %         Href=[mHref Hscalar2(:,peak_voxels(1,5))];     % nullying best found ori
        %         Href2(sn,:)=Hscalar2(:,peak_voxels(1,5)); % tacking Href per inner loop
        %         ref_ori(sn,:)=ori(peak_voxels(1,5),:);
        
        [~,peak_voxels(1,5)]=max(abs(P2.img)); % only finding next largest peak to speed up processing;
        ref_idx=[mref peak_voxels(1,5)];
        Href=[mHref Hscalar2(:,peak_voxels(1,5))];     % nullying best found ori
        Href2(sn,:)=Hscalar2(:,peak_voxels(1,5)); % tacking Href per inner loop
        ref_ori(sn,:)=ori(peak_voxels(1,5),:);
        
        %% According to our 2011 paper (Moiseev et al., 2011 NeuroImage), monitoring the pseudoZ is
        %  the more appropriate and practical means for stopping the iterative loops. When this saturates then most sources should have been found.
        %w=wts2(:,peak_voxels(1,5));
        w=(Rinv*Href)/(Href'*Rinv*Href); % adjusting weights to include all found sources - this should speed things up
        w3=(Rinv*Href)/(Href'*Rinv*Href); % adjusting weights to include all found sources - this should speed things up
        wk=w3(:,end);
        
        % monitoring pseudoZ functions for stopping criteria
        if loc_flag==0
            z=(w'*R*w)./(w'*N*w);  z=diag(z); Zmax(sn)=max(z); % active interval all ref_idx
            nz=(w'*nR*w)./(w'*nN*w);   nz=diag(nz); nZmax(sn)=max(nz); % control interval all ref_idx
            zk(sn)=(wk'*R*wk)/(wk'*N*wk); % current/last ref_idx
            nzk(sn)=(wk'*nR*wk)/(wk'*nN*wk); % current/last ref_idx
            bm_type2='pseudoZ';% pseudoZ
        elseif loc_flag==1
            N2=N/size(Href,2); nN2=nN/size(Href,2);
            z=(w'*Rbar*w)./(w'*N2*w);  z=diag(z); Zmax(sn)=max(z); % active interval all ref_idx
            nz=(w'*nRbar*w)./(w'*nN2*w);  nz=diag(nz); nZmax(sn)=max(nz); % control interval all ref_idx
            zk(sn)=(wk'*Rbar*wk)/(wk'*N2*wk);
            nzk(sn)=(wk'*nRbar*wk)/(wk'*nN2*wk);
            bm_type2='pseudoZ_ER'; % pseudoZ_ER
        elseif loc_flag==2
            z=(w'*Rbar*w)./(w'*R*w);    z=diag(z); Zmax(sn)=max(z); % active interval all ref_idx
            nz=(w'*nRbar*w)./(w'*nR*w);   nz=diag(nz); nZmax(sn)=max(nz); % control interval all ref_idx
            zk(sn)=(wk'*Rbar*wk)/(wk'*R*wk);
            nzk(sn)=(wk'*nRbar*wk)/(wk'*nR*wk);
            bm_type2='pseudoZ_rER';% pseudoZ_rER
        elseif loc_flag==3
            z=(w'*N*w)./(w'*R*w);  z=diag(z); Zmax(sn)=max(z); % active interval all ref_idx
            nz=(w'*nN*w)./(w'*nR*w);  nz=diag(nz); nZmax(sn)=max(nz); % control interval all ref_idx
            zk(sn)=(wk'*N*wk)/(wk'*R*wk);
            nzk(sn)=(wk'*nN*wk)/(wk'*nR*wk);
            bm_type2='activity_index';% activity_index
        end
        
        %z(sn)=img(peak_voxels(1,5));
        pmax(sn)=abs(P2.img(peak_voxels(1,5)));
        pmean(sn)=nanmean(abs(P2.img));
        npmax(sn)=abs(nP.img(peak_voxels(1,5)));
        npmean(sn)=nanmean(abs(nP.img));
        xref(sn)=peak_voxels(1,5);
        if text_flag==1;   fprintf('Sources Found = '); fprintf('%.f ',mref); end
        %        fprintf('\nIteration = %.f\tPeak Voxel#:%.f (ori=%.3f %.3f %.3f) \t%s=%.2f\t%s Threshold =%.2f\t SNR = %.2f\n',zz,xref(sn),ori(xref(sn),1),ori(xref(sn),2),ori(xref(sn),3),bm_type,z(sn),bm_type2,null_thresh,z(sn)/null_thresh);
        %        fprintf('\nIteration = %.f\tPeak Voxel#:%.f (ori=%.3f %.3f %.3f) \t%s=%.2f\t%s Threshold =%.2f\t pseudoSNR = %.2f\n',zz,xref(sn),ori(xref(sn),1),ori(xref(sn),2),ori(xref(sn),3),bm_type,z(sn),bm_type,null_thresh,z(sn)/null_thresh);
        
        %%% new stoping criteria - when pmax doesn't improve more than 'perc_crit' then stop .
        if sn>=2; perc_imp=abs(100*(1-(pmax(end-1)/pmax(end)))); else; perc_imp=100; end
        %% plotting
        if plot_flag==1 && ~isempty(anatomy)
            ln_clr=[0 0.6 0; .8 .4 0; 0 0 0]; % for active, control, diff plots
            grid_locs=anatomy.leadfield.pos;
            slice_orient=[1 1 1]; sFaceAlpha=.5; vx_locs=anatomy.leadfield.voxel_pos;
            inside_idx=find(anatomy.leadfield.inside==1);
            vol_types=3; grid_locs=anatomy.leadfield.pos;
            vol=anatomy.headmodel.bnd(end); opt.vol_nums=1; vol.FaceColor=[.6 .6 .6]; vol.FaceAlpha=0.3; vol.EdgeColor='none';
            %             bl_plot_mesh(vol,opt); bl_plot_mesh(vol)
            
            slice_orient=[1 1 1]; sFaceAlpha=.5; vx_locs=anatomy.leadfield.voxel_pos;
            cmap=jet(255);
            h999=figure(999); clf; set(h999,'position',[20 200 1200 680]); clf;
            img1=P2.img; img1(img1==0)=nan; %min_max1=[min(img1) max(img1)];    % active image
            img2=nP.img; img2(img2==0)=nan; %min_max2=[min(img2) max(img2)];    % control image
            img3=img1-img2;
            min_max1=[min(img1) max(img1)]; min_max2=[min(img2) min(img2)+range(img1)]; min_max3=[min(img3) max(img3)];
            
            % Active Int map
            nthresh=(1-1e-10)*max(abs(img1));   % only finding next largest peak to speed up plotting;
            subplot(3,4,1); cla; haxis = gca; bl_plot_lcmv_peak_img_FT_new(haxis, img1,ori,nthresh,15,vx_locs,cmap,min_max1,vol,[],[0 90],1,vol_types,grid_locs,inside_idx,slice_orient,sFaceAlpha);
            mrk_size=150; s1=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'k+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'ks','sizedata',mrk_size,'linewidth',4);
            mrk_size=150; s1=scatter3(vx_locs(xref(sn),1),vx_locs(xref(sn),2),vx_locs(xref(sn),3),'g+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(xref(sn),1),vx_locs(xref(sn),2),vx_locs(xref(sn),3),'gs','sizedata',mrk_size,'linewidth',4);
            if sn>2; mrk_size=150; s1=scatter3(vx_locs(xref(sn-1),1),vx_locs(xref(sn-1),2),vx_locs(xref(sn-1),3),'+','sizedata',mrk_size,'linewidth',1,'markeredgecolor',[1 0 1]); s2=scatter3(vx_locs(xref(sn-1),1),vx_locs(xref(sn-1),2),vx_locs(xref(sn-1),3),'s','sizedata',mrk_size,'linewidth',4,'markeredgecolor',[1 0 1]); end
            px=get(gca,'position'); a1=axes('position',[px(1)-.05 px(2) .015 px(4)]); plot_colorbar(a1,min_max1,cmap,0)
            subplot(3,4,2); cla; haxis = gca; bl_plot_lcmv_peak_img_FT_new(haxis, img1,ori,nthresh,15,vx_locs,cmap,min_max1,vol,[],[0 0],1,vol_types,grid_locs,inside_idx,slice_orient,sFaceAlpha);
            mrk_size=150; s1=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'k+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'ks','sizedata',mrk_size,'linewidth',4);
            mrk_size=150; s1=scatter3(vx_locs(xref(sn),1),vx_locs(xref(sn),2),vx_locs(xref(sn),3),'g+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(xref(sn),1),vx_locs(xref(sn),2),vx_locs(xref(sn),3),'gs','sizedata',mrk_size,'linewidth',4);
            if sn>2; mrk_size=150; s1=scatter3(vx_locs(xref(sn-1),1),vx_locs(xref(sn-1),2),vx_locs(xref(sn-1),3),'+','sizedata',mrk_size,'linewidth',1,'markeredgecolor',[1 0 1]); s2=scatter3(vx_locs(xref(sn-1),1),vx_locs(xref(sn-1),2),vx_locs(xref(sn-1),3),'s','sizedata',mrk_size,'linewidth',4,'markeredgecolor',[1 0 1]); end
            title(sprintf('Active Interval map (Current Peak Voxel=%.f Value=%.3f)',xref(sn),img1(xref(sn))),'color',[0 .6 0]);
            subplot(3,4,3); cla; haxis = gca; bl_plot_lcmv_peak_img_FT_new(haxis, img1,ori,nthresh,15,vx_locs,cmap,min_max1,vol,[],[-90 0],1,vol_types,grid_locs,inside_idx,slice_orient,sFaceAlpha);
            mrk_size=150; s1=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'k+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'ks','sizedata',mrk_size,'linewidth',4);
            mrk_size=150; s1=scatter3(vx_locs(xref(sn),1),vx_locs(xref(sn),2),vx_locs(xref(sn),3),'g+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(xref(sn),1),vx_locs(xref(sn),2),vx_locs(xref(sn),3),'gs','sizedata',mrk_size,'linewidth',4);
            if sn>2; mrk_size=150; s1=scatter3(vx_locs(xref(sn-1),1),vx_locs(xref(sn-1),2),vx_locs(xref(sn-1),3),'+','sizedata',mrk_size,'linewidth',1,'markeredgecolor',[1 0 1]); s2=scatter3(vx_locs(xref(sn-1),1),vx_locs(xref(sn-1),2),vx_locs(xref(sn-1),3),'s','sizedata',mrk_size,'linewidth',4,'markeredgecolor',[1 0 1]); end
            subplot(3,4,4); cla; hold on; plot(pmax); plot(pmean); plot(npmax); plot(npmean); axis([0 max_iter 0 nanmax(pmax)*1.2]); legend({'Peak Act' 'Mean Act' 'Peak Ctrl' 'Mean Ctrl'},'Location','best'); title('Active PseudoZ values');
            
            % Control Int map
            nthresh=(1-1e-10)*max(abs(img2));
            subplot(3,4,5); cla; haxis = gca; [~,pidx2]=bl_plot_lcmv_peak_img_FT_new(haxis, img2,ori,nthresh,15,vx_locs,cmap,min_max2,vol,[],[0 90],1,vol_types,grid_locs,inside_idx,slice_orient,sFaceAlpha);
            mrk_size=150; s1=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'k+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'ks','sizedata',mrk_size,'linewidth',4);
            mrk_size=150; s1=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'+','sizedata',mrk_size,'linewidth',1,'markeredgecolor',[1 0 0]); s2=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'s','sizedata',mrk_size,'linewidth',4,'markeredgecolor',[1 0 0]);
            px=get(gca,'position'); a1=axes('position',[px(1)-.05 px(2) .015 px(4)]); plot_colorbar(a1,min_max2,cmap,0)
            subplot(3,4,6); cla; haxis = gca; bl_plot_lcmv_peak_img_FT_new(haxis, img2,ori,nthresh,15,vx_locs,cmap,min_max2,vol,[],[0 0],1,vol_types,grid_locs,inside_idx,slice_orient,sFaceAlpha);
            mrk_size=150; s1=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'k+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'ks','sizedata',mrk_size,'linewidth',4);
            mrk_size=150; s1=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'+','sizedata',mrk_size,'linewidth',1,'markeredgecolor',[1 0 0]); s2=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'s','sizedata',mrk_size,'linewidth',4,'markeredgecolor',[1 0 0]);
            title(sprintf('Control Interval map (Current Peak Voxel=%.f Value=%.3f)',pidx2(1),img2(pidx2(1))),'color',[1 0 0]);
            subplot(3,4,7); cla; haxis = gca; bl_plot_lcmv_peak_img_FT_new(haxis, img2,ori,nthresh,15,vx_locs,cmap,min_max2,vol,[],[-90 0],1,vol_types,grid_locs,inside_idx,slice_orient,sFaceAlpha);
            mrk_size=150; s1=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'k+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'ks','sizedata',mrk_size,'linewidth',4);
            mrk_size=150; s1=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'+','sizedata',mrk_size,'linewidth',1,'markeredgecolor',[1 0 0]); s2=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'s','sizedata',mrk_size,'linewidth',4,'markeredgecolor',[1 0 0]);
            if length(z)>1; subplot(3,4,8); cla; hold on; plot(z); plot(nz); axis([1 length(z) 0 max(z)*1.2]); plot([1 length(z)],[null_thresh(sn) null_thresh(sn)],'k--'); legend({'Act Z' 'Ctrl Z' 'Null Thresh'},'Location','best'); title('Inner Loop PseudoZ for refs'); end
            
            
            % SNR map
            nthresh=(1-1e-10)*max(abs(img3));
            subplot(3,4,9); cla; haxis = gca; [~,pidx2]=bl_plot_lcmv_peak_img_FT_new(haxis, img3,ori,nthresh,15,vx_locs,cmap,min_max3,vol,[],[0 90],1,vol_types,grid_locs,inside_idx,slice_orient,sFaceAlpha);
            mrk_size=150; s1=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'k+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'ks','sizedata',mrk_size,'linewidth',4);
            mrk_size=150; s1=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'w+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'ws','sizedata',mrk_size,'linewidth',4);
            px=get(gca,'position'); a1=axes('position',[px(1)-.05 px(2) .015 px(4)]); plot_colorbar(a1,min_max3,cmap,0)
            subplot(3,4,10); cla; haxis = gca; bl_plot_lcmv_peak_img_FT_new(haxis, img3,ori,nthresh,15,vx_locs,cmap,min_max3,vol,[],[0 0],1,vol_types,grid_locs,inside_idx,slice_orient,sFaceAlpha);
            mrk_size=150; s1=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'k+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'ks','sizedata',mrk_size,'linewidth',4);
            mrk_size=150; s1=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'w+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'ws','sizedata',mrk_size,'linewidth',4);
            title('Active-Control map','color',[1 1 1]*0);
            subplot(3,4,11); cla; haxis = gca; bl_plot_lcmv_peak_img_FT_new(haxis, img3,ori,nthresh,15,vx_locs,cmap,min_max3,vol,[],[-90 0],1,vol_types,grid_locs,inside_idx,slice_orient,sFaceAlpha);
            mrk_size=150; s1=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'k+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(mref,1),vx_locs(mref,2),vx_locs(mref,3),'ks','sizedata',mrk_size,'linewidth',4);
            mrk_size=150; s1=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'w+','sizedata',mrk_size,'linewidth',1); s2=scatter3(vx_locs(pidx2(1),1),vx_locs(pidx2(1),2),vx_locs(pidx2(1),3),'ws','sizedata',mrk_size,'linewidth',4);
            if nqloop>1 && length(z9)>1;  subplot(3,4,12); cla; hold on; plot(z9); plot(nz9); axis([1 length(z) 0 max(z9)*1.2]);             plot([1 length(z9)],[null_thresh2(end) null_thresh2(end)],'k--'); legend({'Act Z' 'Ctrl Z' 'Null Thresh'},'Location','best'); title('Outer Loop PseudoZ for refs'); end
        drawnow
        end
        %%
        if text_flag==1; fprintf('\nIteration = %.f\tPeak Voxel#:%.f (ori=%.3f %.3f %.3f) \n%s=%.4f Active=%.4f Theshold=%.4f Control =%.4f pseudoSNR = %.4f Percent improved=%.2f\n',zz,xref(sn),ori(xref(sn),1),ori(xref(sn),2),ori(xref(sn),3),bm_type,pmax(sn),zk(sn),null_thresh(sn),nzk(sn),zk(sn)/nzk(sn),perc_imp); end
        
        %%% new stoping criteria - when pmax doesn't improve more than 'perc_crit' then stop and take max of last 2.
        if sn>=2 && perc_imp<perc_crit; kk=kk+1; nquit=0;
            %          keyboard;
        elseif sn>=2 && perc_imp>perc_crit; kk=0;   % reset if % improvment gets larger than cirterion (i.e., rising again not plateauing)
        end
        
        %%% stopping criterion #1 - reversals of same locations found
        %       if sn>2 && peak_voxels(1,5)==xref(sn-2); kk=kk+1; else kk=0; end
        
        %%% stopping criterion #1- source found more than 5 times.
        %        %mx_ref=max(hist(xref,1:num_voxels));
        %        if max(hist(xref,1:num_voxels))>=num_find; kk=99; end
        
        %%% stopping criterion #2- if sources constantly changing then stop after 20 iterations.
        zz=zz+1;    % if zz>30 then stopping loop
        if zz>max_iter; kk=99; end
        
        %%% stopping criteriaon #3 - for inner loop if there are more than "num_loops" consecutive peak pseudoZ vlaues < null_thresh
        %         if z(sn)<null_thresh; zk=zk+1; else zk=0; end
        %         if zk(sn)<nzk(sn); k2=k2+1; else k2=0; end
        if zk(sn)<null_thresh(sn); k2=k2+1; else k2=0; end
        %         if zk(sn)./nzk(sn)<1; k2=k2+1; else k2=0; end   % stop loop when SNR is less than 1.0 for num_loops
        % This stops finding and including noise sources in the MCMV calculations that could result in canceling "real" sources.
        if k2>num_loops; nquit=1; if text_flag==1; fprintf('No more sources found\n'); end; kk=num_loops; else nquit=0;end
        %         if sum(z./nz<1.0)>0; nquit=1; fprintf('No more sources found\n'); kk=num_loops; else nquit=0;end
    end
    %% OUTER LOOP
    if nquit==0
        %         zk(1:end-num_loops)=0;   % only taking last "num_loops" to find largest peak voxel - more stable estimate of the largest peak.
        %         [~,vv]=max(abs(zk));
        %         mref=[mref xref(vv)];
        %
        zk(1:end-num_loops)=0;   % only taking last "num_loops" to find largest peak voxel - more stable estimate of the largest peak.
        [~,vv]=max(zk); % only taking last "num_loops" to find largest peak voxel - more stable estimate of the largest peak.
        mref=[mref xref(vv)];
        mHref=[mHref Href2(vv,:)']; % adding found peak source leadfields to ref
        ori_MCMV=cat(1,ori_MCMV,ref_ori(vv,:));
    elseif nquit==1 && nqloop>1
        %         fprintf('One source still might be present\n');
        %         [~,vv]=max(abs(z)); % finding largest in all possible sources
        %         mref=[mref xref(vv)];
        if text_flag==1; fprintf('Stopping loop - No more sources!\n'); end
        mm=0;
    end
    
    %     %%% finding best Lead-field when accounting for all found sources.
    %     if length(mref)>1;
    %         [~,mH]=bl_lcmv_scalar_v5(H,mref,[],RNcov,[],loc_flag);
    %         mHref2=mH(:,mref);
    %         %wts_MCMV=zeros(size(wts2));
    %         for v=1:size(mref,2);
    %              H_idx=mref(v);
    %             r_idx=setdiff(1:size(mref,2),v); %selecting all other leadfield for nulling
    %             [~,mH,mori,P3]=bl_lcmv_scalar_v5(H,H_idx,mHref2(:,r_idx),RNcov,mref(r_idx),loc_flag);
    %              mHscalar=mH(:,H_idx); ori_MCMV(v,:)=mori(H_idx,:);
    %             mHref(:,v)=mHscalar;
    %         end
    %     elseif  length(mref)==1;
    %          [~,mH,mori,~]=bl_lcmv_scalar_v5(H,mref,[],RNcov,[],loc_flag);
    %      mHref(:,1)=mH(:,mref); ori_MCMV(1,:)=mori(mref,:);
    %     elseif  isempty(mref);
    %         fprintf('WARNING! No sources found!\n');
    %         mref=[];  mm=0; break;
    %     end
    if  isempty(mref)
        fprintf('WARNING! No sources found!\n');
        mref=[];  mm=0; break;
    elseif length(mref)>=max_sources              % breaking loop after finding maximum # sources
        fprintf('Found maximum number of sources set by user: %.f\n',length(mref)); 
        mm=0; break; % breaking loop because maximum number of sources found as set by the user
    end
    
    % finding if there are sources above null_thresh after accounting for all multisources
    wts=zeros(num_chan,num_voxels);
    wts(:,mref)=(Rinv*mHref)/(mHref'*Rinv*mHref); % Alex's suggested script - returns same as nutMEG script
    for v=1:size(mref,2)
        mpz(v)=(wts(:,mref(v))'*R*wts(:,mref(v)))/(wts(:,mref(v))'*N*wts(:,mref(v))); % MPZ
        mer(v)=(wts(:,mref(v))'*Rbar*wts(:,mref(v)))/(wts(:,mref(v))'*N*wts(:,mref(v))); % MER
        rmer(v)=(wts(:,mref(v))'*Rbar*wts(:,mref(v)))/(wts(:,mref(v))'*R*wts(:,mref(v))); % rMER
        mai(v)=(wts(:,mref(v))'*Ninv*wts(:,mref(v)))/(wts(:,mref(v))'*Rinv*wts(:,mref(v))); % MAI
    end
    wts=zeros(num_chan,num_voxels);
    wts(:,mref)=(Rinv*mHref)/(mHref'*Rinv*mHref); % Alex's suggested script - returns same as nutMEG script
    w=wts(:,mref);
    % monitoring pseudoZ functions for stopping criteria
    if loc_flag==0
        z9=(w'*R*w)./(w'*N*w);  z9=diag(z9); Zmax9(nqloop)=max(z9); % all ref_idx
        nz9=(w'*nR*w)./(w'*nN*w);   nz9=diag(nz9); nZmax9(nqloop)=max(nz9); % all ref_idx
        bm_type2='pseudoZ';% pseudoZ
    elseif loc_flag==1
        N2=N/length(mref); nN2=nN/length(mref);
        z9=(w'*Rbar*w)./(w'*N2*w);  z9=diag(z9); % all ref_idx
        nz9=(w'*nRbar*w)./(w'*nN2*w);  nz9=diag(nz9); % all ref_idx
        bm_type2='pseudoZ_ER'; % pseudoZ_ER
    elseif loc_flag==2
        z9=(w'*Rbar*w)./(w'*R*w);    z9=diag(z9); Zmax9(nqloop)=max(z9); % all ref_idx
        nz9=(w'*nRbar*w)./(w'*nR*w);   nz9=diag(nz9); nZmax9(nqloop)=max(nz9); % all ref_idx
        bm_type2='pseudoZ_rER';% pseudoZ_rER
    elseif loc_flag==3
        z9=(w'*N*w)./(w'*R*w);  z9=diag(z9); Zmax9(nqloop)=max(z9); % all ref_idx
        nz9=(w'*nN*w)./(w'*nR*w);  nz9=diag(nz9); nZmax9(nqloop)=max(nz9); % all ref_idx
        bm_type2='activity_index';% activity_index
    end
    
    % iteratively updating noise to include mHref
    H_idx=1: size(H,3); [~,~,~,nP]=bl_lcmv_scalar_v5(H,H_idx,mHref,nRNcov,mref,loc_flag);
    nd=sort(abs(nP.img(~isnan(nP.img)))); nd=nd(nd~=0); null_thresh2(nqloop)=nd(ceil(length(nd)*(1-noise_alpha)));
    
    r_idx=find(z9>null_thresh2(nqloop)); % MPZ
    %     r_idx=find(z>nz==1); % MPZ
    
    if text_flag==1;  fprintf('Sources = '); fprintf('%.f ',mref); fprintf('\t'); fprintf('Multi Sources = '); fprintf('%.f ',mref(r_idx)); fprintf('\n'); end
    % stopping criteria
    if length(r_idx)<size(mref,2) && ~isempty(r_idx); mm=0; end % stops loop when no more sources can be found above null_thresh after acounting for all multi-sources.
    if max(z9)<=null_thresh2(nqloop);   mm=0; end % stops when the pseudoZ_er value goes below 1 (meaning SNR<1);
end

%% thresholding and calculate final weights for MCMV
m=0;
while m==0
    
    if length(mref)>1
        [~,mH]=bl_lcmv_scalar_v5(H,mref,[],RNcov,[],loc_flag);
        mHref2=mH(:,mref);
        %wts_MCMV=zeros(size(wts2));
        for v=1:size(mref,2)
            H_idx=mref(v);
            r_idx=setdiff(1:size(mref,2),v); %selecting all other leadfield for nulling
            [~,mH,mori,~]=bl_lcmv_scalar_v5(H,H_idx,mHref2(:,r_idx),RNcov,mref(r_idx),loc_flag);
            mHscalar=mH(:,H_idx); ori_MCMV(v,:)=mori(H_idx,:);
            mHref(:,v)=mHscalar;
        end
    elseif  length(mref)==1
        [~,mH,mori,~]=bl_lcmv_scalar_v5(H,mref,[],RNcov,[],loc_flag);
        mHref(:,1)=mH(:,mref); ori_MCMV(1,:)=mori(mref,:);
    elseif  isempty(mref)
        fprintf('WARNING! No sources found!\n');
        mref=[];
        MCMV_idx=[];
        P.type='';
        P.img=[];
        ori=[];
        mHref=[];
        m=99;
    end
    
    
    if ~isempty(mref)
        %         % recalculating localizer value based on only including voxels not below noise threshold.
        %         wts(:,mref)=(Rinv*mHref)/(mHref'*Rinv*mHref); % Alex's suggested script - returns same as nutMEG script
        %         pimg=zeros(1,length(mref));
        %         for v=1:size(mref,2)
        %             if loc_flag==0;  pimg(v)=(wts(:,mref(v))'*R*wts(:,mref(v)))/(wts(:,mref(v))'*N*wts(:,mref(v)));   % MPZ
        %             elseif loc_flag==1; N2=N/length(mref); pimg(v)=(wts(:,mref(v))'*Rbar*wts(:,mref(v)))/(wts(:,mref(v))'*N2*wts(:,mref(v)));    %MER
        %             elseif loc_flag==2; pimg(v)=(wts(:,mref(v))'*Rbar*wts(:,mref(v)))/(wts(:,mref(v))'*R*wts(:,mref(v)));   % rMER
        %             elseif loc_flag==3; pimg(v)=(wts(:,mref(v))'*Ninv*wts(:,mref(v)))/(wts(:,mref(v))'*Rinv*wts(:,mref(v)));    % MAI
        %             end
        %         end
        
        % recalculating localizer value based on only including voxels not below noise threshold.
        wts=zeros(num_chan,num_voxels);
        wts(:,mref)=(Rinv*mHref)/(mHref'*Rinv*mHref); % Alex's suggested script - returns same as nutMEG script
        w=wts(:,mref);
        % monitoring pseudoZ functions for stopping criteria
        if loc_flag==0
            pimg=(w'*R*w)./(w'*N*w);  pimg=diag(pimg);
            npimg=(w'*nR*w)./(w'*nN*w);   npimg=diag(npimg);
            bm_type2='pseudoZ';% pseudoZ
        elseif loc_flag==1
            N2=N/length(mref); nN2=nN/length(mref);
            pimg=(w'*Rbar*w)./(w'*N2*w);  pimg=diag(pimg); % all ref_idx
            npimg=(w'*nRbar*w)./(w'*nN2*w);  npimg=diag(npimg); % all ref_idx
            bm_type2='pseudoZ_ER'; % pseudoZ_ER
        elseif loc_flag==2
            pimg=(w'*Rbar*w)./(w'*R*w);    pimg=diag(pimg); Zmax9(nqloop)=max(pimg); % all ref_idx
            npimg=(w'*nRbar*w)./(w'*nR*w);   npimg=diag(npimg); nZmax9(nqloop)=max(npimg); % all ref_idx
            bm_type2='pseudoZ_rER';% pseudoZ_rER
        elseif loc_flag==3
            pimg=(w'*N*w)./(w'*R*w);  pimg=diag(pimg); Zmax9(nqloop)=max(pimg); % all ref_idx
            npimg=(w'*nN*w)./(w'*nR*w);  npimg=diag(npimg); nZmax9(nqloop)=max(npimg); % all ref_idx
            bm_type2='activity_index';% activity_index
        end
        
        % Final updating of noise to include mHref
        if ~isempty(mref)
            %             H_idx=1: size(H,3); [~,~,~,nP]=bl_lcmv_scalar_v5(H,H_idx,mHref,nRNcov,mref,loc_flag);      % this is conservative.
            %             nd=sort(abs(nP.img(~isnan(nP.img)))); nd=nd(nd~=0); null_thresh=nd(ceil(length(nd)*(1-noise_alpha)));
            
            cov_flag = 1;   % using Rinv for weights
            [nulled_wts]=calc_nulled_wts(H,mHref,RNcov,mref,loc_flag,ori_MCMV,cov_flag);
            [pimg2,npimg2]=calc_localizer(nulled_wts,mref,RNcov,nRNcov,loc_flag);
            nd=sort(abs(npimg2(~isnan(npimg2)))); nd=nd(nd~=0);  nulled_thresh=nd(ceil(length(nd)*(1-noise_alpha)));
            
            
        end
        
        if isempty(find(pimg2(mref)<=nulled_thresh, 1))    % breaking outer loop because all voxels are above threshold now
            m=99;
        elseif isempty(find(pimg2(mref)>nulled_thresh,1))       % no significant voxels after thresholding
            fprintf('WARNING! No sources found!\n');
            mref=[];pimg=[]; ori_MCMV=[]; MCMV_idx=[]; P=[]; P.type=''; P.img=[]; ori=[]; mHref=[]; wts=[];
            m=99;
        else
            r_idx=find(pimg2(mref)>nulled_thresh);
            mref=mref(r_idx);
            ori_MCMV=ori_MCMV(r_idx,:);
            mHref=mHref(:,r_idx);
        end
    else
        fprintf('WARNING! No sources found!\n');
        mref=[];pimg=[];npimg=[]; ori_MCMV=[]; MCMV_idx=[]; P=[]; P.type=''; P.img=[]; ori=[]; mHref=[]; wts=[];
        m=99;
    end
    
end

P.img=zeros(num_voxels,1);
nP.img=zeros(num_voxels,1);
MCMV_idx=mref;

%% Residual maps
cov_flag = 1;   % using Rinv for weights
[residual_wts,residual_Hscalar,residual_ori]=calc_nulled_wts(H,mHref,RNcov,mref,loc_flag,ori_MCMV,cov_flag);
[pimg3,npimg3,bmf_type]=calc_localizer(residual_wts,mref,RNcov,nRNcov,loc_flag);
nd=sort(abs(npimg3(~isnan(npimg3)))); nd=nd(nd~=0);  residual_thresh=nd(ceil(length(nd)*(1-noise_alpha)));

%% % Final image and weights
cov_flag = 2;   % using Ninv for weights
[nulled_wts,nulled_Hscalar,nulled_ori]=calc_nulled_wts(H,mHref,RNcov,mref,loc_flag,ori_MCMV,cov_flag);
[pimg2,npimg2,bmf_type]=calc_localizer(nulled_wts,mref,RNcov,nRNcov,loc_flag);
nd=sort(abs(npimg2(~isnan(npimg2)))); nd=nd(nd~=0);  nulled_thresh=nd(ceil(length(nd)*(1-noise_alpha)));


%%
%  H_idx=1: size(H,3); [f_wts,f_Hscalar,f_ori,f_P]=bl_lcmv_scalar_v5(H,H_idx,mHref,RNcov,mref,loc_flag);      % this is conservative.
%  f_wts(:,mref)=wts(:,mref);
% %  f_Hscalar(:,mref)=mHref;
%  f_ori(mref,:)=ori_MCMV;
%  f_P.img(mref)=pimg;

% if loc_flag==0; % MPZ
%     P.type='MPZ';
%     P.img(MCMV_idx)=pimg; %(wts(:,v)'*R*wts(:,v))/(wts(:,v)'*N*wts(:,v));
% elseif loc_flag==1;  %MER
%     P.type='MER';
%     P.img(MCMV_idx)=pimg; %(wts(:,v)'*Rbar*wts(:,v))/(wts(:,v)'*N*wts(:,v));
% elseif loc_flag==2; % rMER
%     P.type='rMER';
%     P.img(MCMV_idx)=pimg; %(wts(:,v)'*Rbar*wts(:,v))/(wts(:,v)'*N*wts(:,v));
% elseif loc_flag==3; % MAI
%     P.type='MAI';
%     P.img(MCMV_idx)=pimg; %(wts(:,v)'*Rbar*wts(:,v))/(wts(:,v)'*N*wts(:,v));
% end

if ~isempty(pimg)
    P.img(MCMV_idx)=pimg;
    P.img(pimg<=null_thresh)=0; % double checking to remove all voxel values below null_thresh
    nP.img(MCMV_idx)=npimg;
    nP.img(pimg<=null_thresh)=0; % double checking to remove all voxel values below null_thresh
    % ori=ori_MCMV;
    
    fprintf('Final Results\nSources \t Active \t Control \t SNR'); fprintf('\n');
    for v=1:length(MCMV_idx); fprintf('%.f \t\t %.3f \t\t %.3f \t\t %.3f\n',MCMV_idx(v),pimg2(MCMV_idx(v)),npimg2(MCMV_idx(v)),pimg2(MCMV_idx(v))./npimg2(MCMV_idx(v)) ); end
    fprintf('Threshold \t %.3f\n',nulled_thresh);
    
    fprintf('Final Results\nSources \t Active \t Control \t SNR'); fprintf('\n');
    for v=1:length(MCMV_idx); fprintf('%.f \t\t %.3f \t\t %.3f \t\t %.3f\n',MCMV_idx(v),pimg(v),npimg(v),pimg(v)./npimg(v) ); end
    fprintf('Threshold \t %.3f\n',null_thresh);
else
    fprintf('No sources were found in Final Solution!\n');
end


%% OUTPUTS
MIA.RNcov=RNcov;
MIA.MCMV_idx=MCMV_idx;
MIA.MCMV_ori=ori_MCMV;
MIA.wts=wts;
MIA.mHref=mHref;
MIA.ori=ori_MCMV;
MIA.P.type=bmf_type;
MIA.P.img=P.img;
MIA.P_ctrl.type=bmf_type;
MIA.P_ctrl.img=nP.img;
MIA.null_thresh=null_thresh;

% MCMV + Residual maps --> using Rinv to create maps
MIA.residual_wts=residual_wts;
MIA.residual_Hscalar=residual_Hscalar;
MIA.P.residual_img=pimg3;
MIA.P_ctrl.residual_img=npimg3;
MIA.residual_thresh=residual_thresh;
MIA.residual_ori = residual_ori;

% nulled --> using Ninv to create maps
MIA.nulled_wts=nulled_wts;
MIA.nulled_Hscalar=nulled_Hscalar;
MIA.P.nulled_img=pimg2;
MIA.P_ctrl.nulled_img=npimg2;
MIA.nulled_thresh=nulled_thresh;
MIA.nulled_ori = nulled_ori;

% [nulled_wts,nulled_Hscalar]=calc_nulled_wts(H,mHref,RNcov,mref,loc_flag,mref_ori);
% [pimg,npimg,bmf_type]=calc_localizer(nulled_wts,mref,RNcov,nRNcov,loc_flag);

function [nulled_wts,nulled_Hscalar,vx_ori]=calc_nulled_wts(H,mHref,RNcov,mref,loc_flag,mref_ori,cov_flag)
% H_idx=1:size(H,3); [~,~,vx_ori,~]=bl_lcmv_scalar_v5(H,H_idx,[],RNcov,[],loc_flag);
H_idx=1:size(H,3); [~,~,vx_ori,~]=bl_lcmv_scalar_v5(H,H_idx,mHref,RNcov,mref,loc_flag);
if ~isempty(mref) && ~isempty(mref_ori)
    vx_ori(mref,:)=mref_ori; % replacing SPA_ori with MIA_ori
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



