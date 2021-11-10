
%% INFANT Anatomy
load('C:\BRANELab\matlab_progs\general_progs\EEG_sim\SimSignals_GUI\anatomy\ANATOMY_DEFAULT_BRANELab_EEG_BioSemi_66ch_MRIvolume_1yr_infant.mat');

a=fieldnames(anatomy);
for k=1:length(a)
    eval(sprintf('%s = anatomy.%s;',a{k},a{k}))
end
save('C:\BRANELab\matlab_progs\general_progs\EEG_sim\SimSignals_GUI\anatomy\ANATOMY_DEFAULT_BRANELab_MEEG_1yr_infant.mat','mri', 'vol', 'elec', 'please_read', 'pial_mesh', 'sens_adult_dewar', 'sens_infant_dewar', 'sens_adult_dewar_supine', 'mesh_hd', 'bnd', 'elec_good', 'chan_idx', 'vol_leadfield', 'mesh_sourcemodel', 'mesh_leadfield', 'leadfield', 'vol2', 'vx_locs', 'eeg_vol_leadfield', 'vol_infant_dewar', 'vol_adult_dewar', 'vol_adult_dewar_supine', 'meg_vol_leadfield', 'vol_eeg');

% adult dewar sitting
anatomy.sens_adult_dewar.coilori = anatomy.sens_adult_dewar.coilori(1:size(anatomy.sens_adult_dewar.label,1),:);
anatomy.sens_adult_dewar.coilpos = anatomy.sens_adult_dewar.coilpos(1:size(anatomy.sens_adult_dewar.label,1),:);
anatomy.sens_adult_dewar.tra     = anatomy.sens_adult_dewar.tra(1:size(anatomy.sens_adult_dewar.label,1),1:size(anatomy.sens_adult_dewar.label,1));
% adult dewar supine
anatomy.sens_adult_dewar_supine.coilori = anatomy.sens_adult_dewar_supine.coilori(1:size(anatomy.sens_adult_dewar_supine.label,1),:);
anatomy.sens_adult_dewar_supine.coilpos = anatomy.sens_adult_dewar_supine.coilpos(1:size(anatomy.sens_adult_dewar_supine.label,1),:);
anatomy.sens_adult_dewar_supine.tra     = anatomy.sens_adult_dewar_supine.tra(1:size(anatomy.sens_adult_dewar_supine.label,1),1:size(anatomy.sens_adult_dewar_supine.label,1));
% infant dewar
anatomy.sens_infant_dewar.coilori = anatomy.sens_infant_dewar.coilori(1:size(anatomy.sens_infant_dewar.label,1),:);
anatomy.sens_infant_dewar.coilpos = anatomy.sens_infant_dewar.coilpos(1:size(anatomy.sens_infant_dewar.label,1),:);
anatomy.sens_infant_dewar.tra     = anatomy.sens_infant_dewar.tra(1:size(anatomy.sens_infant_dewar.label,1),1:size(anatomy.sens_infant_dewar.label,1));


%% ADULT Anatomy
sens_adult = anatomy.sens_adult_dewar; % first load in infant anatomy above 
anatomy2 = load('C:\BRANELab\matlab_progs\general_progs\EEG_sim\SimSignals_GUI\anatomy\ANATOMY_DEFAULT_Biosemi_72ch_Tujillo.mat');
% % re-aligning sensors to adult dewar <-- but no need to do this because they were already aligned
% cfg           = [];
% cfg.method    = 'interactive';
% cfg.elec    = sens_adult;
% % cfg.channel      = sens_idx;
% cfg.headshape = anatomy_adult.scalp.bnd(1);
% sens_aligned  = ft_electroderealign(cfg,sens_adult);
anatomy2.sens_adult=sens_adult;
a=fieldnames(anatomy2);
for k=1:length(a)
    eval(sprintf('%s = anatomy2.%s;',a{k},a{k}))
end
save('C:\BRANELab\matlab_progs\general_progs\EEG_sim\SimSignals_GUI\anatomy\ANATOMY_DEFAULT_MEEG_adult_Tujillo.mat','mri', 'vol', 'elec', 'please_read', 'pial_mesh', 'sens_adult_dewar', 'sens_infant_dewar', 'sens_adult_dewar_supine', 'mesh_hd', 'bnd', 'elec_good', 'chan_idx', 'vol_leadfield', 'mesh_sourcemodel', 'mesh_leadfield', 'leadfield', 'vol2', 'vx_locs', 'eeg_vol_leadfield', 'vol_infant_dewar', 'vol_adult_dewar', 'vol_adult_dewar_supine', 'meg_vol_leadfield', 'vol_eeg');


%% Plot anatomy
figure(1); clf;

% Infant adult dewar supine
subplot(2,2,1); cla; 
ft_plot_mesh(anatomy.vol.bnd(3),'facecolor',[1 .6 .4]*.75,'edgecolor',[1 .6 .4]*.85,'facealpha',.4,'edgealpha',.4);
ft_plot_mesh(anatomy.pial_mesh,'facecolor',[1 1 1]*.8,'edgecolor','none','facealpha',1,'edgealpha',1);
ft_plot_sens(anatomy.sens_adult_dewar_supine,'label','off', 'chantype','meg','edgecolor','k'); %,'label','off');
ft_plot_sens(anatomy.elec,'label','off', 'chantype','eeg','facecolor','b','edgecolor','b'); %,'label','off');
view(-90,0); c=camlight; lighting gouraud; material dull; axis tight; p_axis = [get(gca,'XLim') get(gca,'YLim')];

% Infant adult dewar sitting
subplot(2,2,2); cla; 
ft_plot_mesh(anatomy.vol.bnd(3),'facecolor',[1 .6 .4]*.75,'edgecolor',[1 .6 .4]*.85,'facealpha',.4,'edgealpha',.4);
ft_plot_mesh(anatomy.pial_mesh,'facecolor',[1 1 1]*.8,'edgecolor','none','facealpha',1,'edgealpha',1);
ft_plot_sens(anatomy.sens_adult_dewar,'label','off', 'chantype','meg','edgecolor','k'); %,'label','off');
ft_plot_sens(anatomy.elec,'label','off', 'chantype','eeg','facecolor','b','edgecolor','b'); %,'label','off');
view(-90,0); c=camlight; lighting gouraud; material dull;  axis tight; %axis( p_axis )

% Infant infant dewar
subplot(2,2,3); cla; 
ft_plot_mesh(anatomy.vol.bnd(3),'facecolor',[1 .6 .4]*.75,'edgecolor',[1 .6 .4]*.85,'facealpha',.4,'edgealpha',.4);
ft_plot_mesh(anatomy.pial_mesh,'facecolor',[1 1 1]*.8,'edgecolor','none','facealpha',1,'edgealpha',1);
ft_plot_sens(anatomy.sens_infant_dewar,'label','off', 'chantype','meg','edgecolor','k'); %,'label','off');
ft_plot_sens(anatomy.elec,'label','off', 'chantype','eeg','facecolor','b','edgecolor','b'); %,'label','off');
view(-90,0); c=camlight; lighting gouraud; material dull;  axis( p_axis )

% Adult dewar
subplot(2,2,4); cla; 
ft_plot_mesh(anatomy2.scalp.bnd,'facecolor',[1 .6 .4]*.75,'edgecolor',[1 .6 .4]*.85,'facealpha',.4,'edgealpha',.4);
ft_plot_mesh(anatomy2.pial_mesh,'facecolor',[1 1 1]*.8,'edgecolor','none','facealpha',1,'edgealpha',1);
ft_plot_sens(anatomy2.sens_adult,'label','off', 'chantype','meg','edgecolor','k'); %,'label','off');
ft_plot_sens(anatomy2.elec,'label','off', 'chantype','eeg','facecolor','b','edgecolor','b'); %,'label','off');
view(-90,0); c=camlight; lighting gouraud; material dull;  axis tight

