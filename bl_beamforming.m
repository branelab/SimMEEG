function bl_beamforming(data, anat, beamtype, mrifile, grid_res);

% INPUT
%   beamtype =  beamformer type
%                   (0) = LCMV SPA
%                   (1) = MCMV SIA
%                   (2) = MCMV MIA
%   T1file = MRI template file for warping subject MRI into template space
%            (e.g., T1file='C:\usr\local\matlab_progs\fieldtrip-20161103\external\spm8\templates\T1.nii').
%   mrifile = subject MRI file name (e.g., 'S01.mri' - CTF format);
%   grid_res = voxel grid resolution (mm) basd on FT's source templates, thus must have one of ollowing values = 4, 5, 6, 7, 8, or 10 mm.
%   anat = structure of anatomy read using FieldTrip's programs (if empty then program creates this variable).
%           anat.mri        = ft_read_mri(mrifile);
%           anat.sens       = ft_read_sens(dsfile);
%           anat.seg        = ft_volumesegment(cfg, anat.mri);
%           anat.headmodel  = ft_prepare_headmodel(cfg, anat.seg);
%           anat.grid       = ft_prepare_sourcemodel(cfg);
%           anat.grid_lf    = ft_prepare_leadfield(cfg);
%   

%% INPUT variables
%addpath(genpath('C:\usr\local\matlab_progs\fieldtrip-20161103\'))
%rmpath(genpath('C:\usr\local\matlab_progs\fieldtrip-20161103\external\'));

%T1file='C:\usr\local\matlab_progs\fieldtrip-20161103\external\spm8\templates\T1.nii';
mrifile ='\\PONYO\Data\SESAEF\SESAEF_Adult\G003\G003.mri';
dsfile='\\PONYO\Data\SESAEF\SESAEF_Adult\G003\G003_SESAEF_SingleTrial.ds';
fdata=epoch_data(1).data(:,ds.MEG_chan_idx,epoch_data(1).good_trial);
fdata=filter_data(fdata,600,[1 20],'butter','bandpass',2);

ctrl_samps=960:1200;    % -400 to 0 ms
ctrl_samps=1080:1200;    % -400 to 0 ms
act_samps=1231:1351;    % 50 to 250 ms
loc_flag=1; % 1= zER/MER evoked localizer,  0 = pseduoZ/MPZ
plot_flag=0;    % plot results
perc_crit=2;    % percent improved to stop iterative steps used in SIA and MIA


if isempty(anat);
%% load Anatomy
anat.mri = ft_read_mri(mrifile);
anat.sens = ft_read_sens(dsfile);
anat.sens = ft_convert_units(anat.sens,'mm');

% segmenting brain, skull, and scalp
cfg        = [];
cfg.output = {'brain', 'skull', 'scalp'};
anat.seg        = ft_volumesegment(cfg, anat.mri);

% % getting meshes for skull & brain
% cfg = [];
% cfg.method = 'projectmesh';
% cfg.tissue = 'brain';
% cfg.numvertices = 3000;
% anat.mesh_brain = ft_prepare_mesh(cfg, anat.seg);
% 
% cfg = [];
% cfg.method = 'projectmesh';
% cfg.tissue = 'skull';
% cfg.numvertices = 2000;
% anat.mesh_skull = ft_prepare_mesh(cfg, anat.seg);
% 
% cfg = [];
% cfg.method = 'projectmesh';
% cfg.tissue = 'scalp';
% cfg.numvertices = 1000;
% anat.mesh_scalp = ft_prepare_mesh(cfg, anat.seg);

% % plotting
% figure
% ft_plot_mesh(mesh_scalp, 'edgecolor', 'none', 'facecolor', 'skin')
% material dull
% camlight
% lighting phong
% print -dpng natmeg_dip_scalp.png


% construct the volume conductor model (i.e. head model) for each subject
cfg        = [];
cfg.method = 'singleshell'; %'openmeeg';
anat.headmodel  = ft_prepare_headmodel(cfg, anat.seg);

% getting template grid
% ftpath   = 'C:\usr\local\matlab_progs\fieldtrip-20161103\'; % this is the path to fieldtrip at Donders
% load(fullfile(ftpath, '\template\sourcemodel\standard_sourcemodel3d5mm'));
% template_grid = sourcemodel;
% clear sourcemodel;

% create the subject specific grid, using the template grid that has just been created
cfg                = [];
cfg.grid.warpmni   = 'yes';
%cfg.grid.template  = template_grid;
cfg.grid.resolution = grid_res;    % in mm
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = anat.mri;
cfg.grid.unit      ='mm';
anat.grid          = ft_prepare_sourcemodel(cfg);

%% calculate LeadField matrix
% Create leadfield grid
meg_idx=strmatch('meggrad',anat.sens.chantype);
cfg                 = [];
cfg.grad            = anat.sens;
cfg.vol             = anat.headmodel;
cfg.grid            = anat.grid;
%cfg.method = 'bem_openmeeg';
%cfg.reducerank = 2; % default for MEG is 2, for EEG is 3
%cfg.normalize       = 'yes';
%cfg.normalizeparam  = 0.5; % depth normalization parameter (default = 0.5)
%cfg.grid.resolution = 0.5;   % use a 3-D grid with a 1 cm resolution
%cfg.grid.unit       = 'cm';
%cfg.grid.tight      = 'yes';
[anat.grid_lf] = ft_prepare_leadfield(cfg);
end

%% make a figure of the single subject headmodel, and grid positions
figure(999); clf; set(gcf,'color',[.6 .6 .6]);  hold on;
ft_plot_vol(anat.headmodel, 'facecolor',[1 .6 1],'edgecolor', 'none', 'facealpha', 0.4);
ft_plot_mesh(anat.grid.pos(anat.grid.inside,:),'facecolor','none','edgecolor', 'none', 'vertexcolor',[0 0 0]);
ft_plot_sens(anat.sens, 'facecolor','none','edgecolor', [1 1 0], 'facealpha', 0.4); axis([-200 200 -200 200 -50 200]); axis off;

%% restructuring LeadField for MCMV
anat.H=nan(length(meg_idx),3,size(anat.grid_lf.leadfield,2));
anat.H_idx=nan(size(anat.grid_lf.leadfield,2),1);
for vx=1:size(anat.grid_lf.leadfield,2)
    if ~isempty(anat.grid_lf.leadfield{vx})
        x=anat.grid_lf.leadfield{vx};
        anat.H(:,:,vx)=x(meg_idx,:);
        anat.H_idx(vx)=vx;
    end
end
%vx_idx=find(~isnan(H_idx)); % voxels within brain that have leadfields
anat.vx = ft_inside_vol(anat.grid_lf.pos, anat.headmodel); % voxels within brain that have leadfields
anat.vx_idx=find(anat.vx==1);


%% calculagting Beamformers

rmpath(genpath('C:\usr\local\matlab_progs\fieldtrip-20161103\'));

mHref=[];mref=[];
[LCMV.wts,LCMV.Hscalar,LCMV.ori,LCMV.P,LCMV.RNcov,LCMV.null_thresh]=BRANElab_LCMV_beamformer_SPA(anat.H(:,:,anat.vx_idx),fdata,act_samps,ctrl_samps,loc_flag);
[MIA]=BRANElab_MCMV_beamformer_MIA(anat.H(:,:,anat.vx_idx),mHref,mref,fdata,act_samps,ctrl_samps,anat.grid_lf.pos(anat.vx_idx,:),loc_flag,plot_flag,perc_crit);

addpath(genpath('C:\usr\local\matlab_progs\fieldtrip-20161103\'))
rmpath(genpath('C:\usr\local\matlab_progs\fieldtrip-20161103\external\'));
paxis=[-500 800 -50 80];
%paxis=[-500 800 -20 20];
figure(991); clf; set(gcf,'color','w');


%% plot MCMV
MIA_img=zeros(length(anat.vx_idx),1);
MIA_img(MIA.MCMV_idx)=1;
subplot(2,2,2); cla; bl_plot_lcmv_peak_img_FT(MIA_img,0,1,anat.grid_lf.pos(anat.vx_idx,:),lines(256),[0 1],anat.headmodel,anat.sens,1); view(0,0);
title('MCMV (MIA)','color','k');
MCMV_swf=[MIA.wts(:,MIA.MCMV_idx)'*squeeze(nanmean(fdata,3))']*1e14;
subplot(2,2,4);cla;
plot(epoch_data(1).lat*1000,squeeze(nanmean(fdata,3))*1e14,'color',[1 1 1]*.6);
hold on;
plot(epoch_data(1).lat*1000,squeeze(MCMV_swf)); axis(paxis)

%% plot LCMV
s=sort(abs(LCMV.P.img));
img_thresh=s(round(length(LCMV.P.img)*.95));
%img_thresh=LCMV.null_thresh/sqrt(2);
subplot(2,2,1); cla; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT(LCMV.P.img,img_thresh,15,anat.grid_lf.pos(anat.vx_idx,:),jet(256),[-(max(abs(LCMV.P.img))) (max(abs(LCMV.P.img)))],anat.headmodel,anat.sens,1); view(0,0);
%subplot(2,2,2); cla; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT(LCMV.P.img,img_thresh,15,grid_lf.pos(vx_idx,:),jet(256),[-(max(abs(LCMV.P.img))) (max(abs(LCMV.P.img)))],headmodel,sens,1); view(180,0);
%subplot(2,4,2); cla; [peak_voxels,p_idx]=bl_plot_lcmv_peak_img_FT(LCMV.P.img,LCMV.null_thresh,15,grid_lf.pos(vx_idx,:),jet(256),[-(max(abs(LCMV.P.img))) (max(abs(LCMV.P.img)))],headmodel,sens,0); view(0,0);
title('LCMV','color','k');
LCMV_swf=[LCMV.wts(:,p_idx)'*squeeze(nanmean(fdata,3))']*1e14;
subplot(2,2,3);cla;
plot(epoch_data(1).lat*1000,squeeze(nanmean(fdata,3))*1e14,'color',[1 1 1]*.6);
hold on;
plot(epoch_data(1).lat*1000,squeeze(LCMV_swf)); 
axis(paxis);

return


