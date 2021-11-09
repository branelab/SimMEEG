% Field Trip Tutorial on Beamforming Motor Evoked Fields
%% Defining Pre and Post Stimulus interval using ctrl_int and act_int
    cfg = [];
    cfg.toilim = ctrl_int;
    datapre = ft_redefinetrial(cfg,ft_data);
    cfg.toilim = act_int;
    datapost = ft_redefinetrial(cfg, ft_data);
    
    cfg = [];
    cfg.covariance='yes';
    cfg.covariancewindow = [ctrl_int(1) act_int(2)]; % calculating covariance across entire ctrl and act interval just as stated in Field Trip's Tutorial
    avg_data = ft_timelockanalysis(cfg,ft_data);
    
    cfg = [];
    cfg.covariance='yes';
    avgpre = ft_timelockanalysis(cfg,datapre);
    avgpst = ft_timelockanalysis(cfg,datapost);
    
    
    %% first call to ft_sourceanalysis keeping the spatial filters
    cfg=[];
    cfg.grad = h.inv_soln(h.next_inv_soln).elec_good;
    cfg.method = 'lcmv';
    cfg.projectnoise = 'yes'; %     = project noise estimate through filter,           can be 'yes' or 'no'
    cfg.grid   = h.inv_soln(h.next_inv_soln).leadfield;
    cfg.headmodel    = h.inv_soln(h.next_inv_soln).headmodel;
    cfg.lcmv.keepfilter='yes';
    cfg.lcmv.lambda = '5%';
    LCMV = ft_sourceanalysis(cfg, avg_data);
    
    %% second and third call to ft_sourceanalysis now applying the precomputed filters to pre and %post intervals
    cfg=[];
    cfg.method='lcmv';
    cfg.grad = h.inv_soln(h.next_inv_soln).elec_good;
    cfg.grid   = h.inv_soln(h.next_inv_soln).leadfield;
    cfg.sourcemodel.filter = LCMV.avg.filter;
    cfg.headmodel    = h.inv_soln(h.next_inv_soln).headmodel;
    sourcepre = ft_sourceanalysis(cfg, avgpre);
    sourcepst = ft_sourceanalysis(cfg, avgpst);
    
    img = (sourcepst.avg.pow-sourcepre.avg.pow)./sourcepre.avg.pow;
    LCMV.P.img = img(inside_idx);
    
    y=cell2mat(LCMV.avg.filter); wts=permute(reshape(y,[3 size(y,1)/3 size(y,2)]),[2 1 3]);
    LCMV.swf=[]; LCMV.swf_snr=[]; LCMV.swf_pwr=[]; LCMV.noise_pwr=[];
    for ox=1:3
        LCMV.swf(:,ox,:)=avg_data.avg'*squeeze(wts(:,ox,:))';
        %    swf_snr(:,ox,:)=20*log10(bsxfun(@rdivide,swf(:,ox,:),nanmean(swf(h.cfg.study.base_samps,ox,:),1)));
        %    LCMV.swf_snr(:,ox,:)=(bsxfun(@rdivide,LCMV.swf(:,ox,:),nanmean(LCMV.swf(h.cfg.study.base_samps,ox,:),1))); % SNR in percent from baseline
        LCMV.swf_snr(:,ox,:)=(bsxfun(@rdivide,LCMV.swf(:,ox,:),nanstd(LCMV.swf(h.cfg.study.base_samps,ox,:),1))); % normlaized to baseline
        LCMV.swf_pwr(ox,:)=rms(LCMV.swf(h.cfg.study.bl_bmf.act_samps,ox,:),1)./rms(LCMV.swf(h.cfg.study.bl_bmf.ctrl_samps,ox,:),1);
        LCMV.noise_pwr(ox,:)=rms(LCMV.swf(h.cfg.study.bl_bmf.ctrl_samps,ox,:),1)./rms(LCMV.swf(h.cfg.study.bl_bmf.ctrl_samps,ox,:),1);
    end
    LCMV.P.img2= LCMV.avg.pow(inside_idx)./LCMV.avg.noise(inside_idx);   % averaging across orientations
    nd=sort(LCMV.P.img); LCMV.null_thresh=nd(ceil(length(nd)*h.cfg.study.bl_bmf.noise_alpha)); %img(img<pthresh)=0;
    % Approximating orientations for plotting purposes only.
    x=LCMV.swf_pwr(1,:)./max(LCMV.swf_pwr,[],1);
    y=LCMV.swf_pwr(2,:)./max(LCMV.swf_pwr,[],1);
    z=LCMV.swf_pwr(3,:)./max(LCMV.swf_pwr,[],1);
    LCMV.ori=[x;y;z]';
    
    
    
      img=LCMV.P.img; ori=LCMV.ori; null_thresh=max(img)*.75;
        min_max=[min(img) max(img)*.95];  vol_types=1;
        seed_idx = 1:3;  ln_wdth = 1; ln_wdth2 = 2;
        figure(1004); clf; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT_new(img,ori,null_thresh,15,h.anatomy.leadfield.voxel_pos,jet(255),min_max,h.anatomy.vol.bnd(3),[],...
            h.cfg.study.bl_bmf.vw_angle,1,vol_types,h.anatomy.leadfield.pos,inside_idx,h.cfg.study.bl_bmf.slice_orient,h.cfg.study.bl_bmf.sFaceAlpha);
        hold on; mrk_size=150; s1=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'k+','sizedata',mrk_size,'linewidth',ln_wdth); s2=scatter3(h.cfg.source.vx_locs(seed_idx,1),h.cfg.source.vx_locs(seed_idx,2),h.cfg.source.vx_locs(seed_idx,3),'ks','sizedata',mrk_size,'linewidth',ln_wdth2);
        title('FT LCMV'); colorbar; caxis(min_max);
    
    
    
        