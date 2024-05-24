function bl_create_hdm_leadfields()


HeadModelMat = load('D:\brainstorm_db\SimMEEG_Anatomy\data\chan66\NewCondition\headmodel_vol_openmeeg.mat');
ChannelMat = load('D:\brainstorm_db\SimMEEG_Anatomy\data\chan66\NewCondition\channel_BioSemi_64_A01_02.mat');
[elec, grad] = out_fieldtrip_channel(ChannelMat, 1);


[ftHeadmodel, ftLeadfield, iChannels] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, 1:66);


mesh = ftHeadmodel.bnd;
hdm = ft_headmodel_concentricspheres(mesh);
hdm.skinsurfaces = 3;
hdm.innerskullsurface = 1; 
hdm.unit = 'm';
hdm.label = ftHeadmodel.label;

leadfieldopt = {'reducerank', [3], 'backproject', [], 'normalize', [], 'normalizeparam', [], 'weight', []};
lf = ft_compute_leadfield([5 -30 20], elec, hdm, leadfieldopt{:});
