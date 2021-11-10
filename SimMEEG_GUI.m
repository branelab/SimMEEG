function SimMEEG_GUI(varargin)
%% GUI for user-friendly selection of parameters needed to run SimSignals.m
%    varagin = bst structure called directly from Brainstorm see "sm_bst2ft_anatomy_from_bst_files" for input structure
%
% This version is for general distribution
%
% Updates since v19:
%   - Major Bug Fixed - Source modeling was setting "Simulated" forward leadfields to the "inverse Head Model" selection --> this yields wrong location errors
%   - added button "Plot SpatioTemporal" to plot spatiotemporal mapping of inverse solutions - can then click on source waveform axis to move time cursor
%   - fixed monte-carlo bug that allows it to run even when location range dimenions were different across X, Y, Z edit boxes
%   - Monte-Carlo Tab: a new look --> reconfigured Monte-Carlo panels
%   - Mote-Carlo tab: added head model selection for forward and inverse modeling
%   - fixed source waveform calculations in bs_calc_FC.m for MNE and eLORETA by changing (permute) leadfields to needed dimensions
%   - fixed "Orientation Constraint" so that user can select "Random" for "Cortical Surface" HeadModels
%   - added button "plot LeadField" to plot dipole locations for leadfields and transparency slider for them
%   - added a radio button to plot raw or normalized source waveforms from inverse solutions
%   - fixed "Noise Projection" menu sub-selection items in panel "Generator Sensor Noise" on the "Simulate M/EEG" Panel
%   - added ability to plot evoked as butterfly plot by selecting multiple channels in the list box on the "Simulate M/EEG" Panel
%   - added "Real Sensor" in "Noise Projection" on the "Simulate M/EEG" Panel - this allows users to "Load Real Noise" from a datafile.mat with variable data=[samples x channels x trials].
%   - added button that with randomize the "Real Sensor" noise loaded by apply Alex Moiseev's PCA-FFT phase randomization procedure --> this removes sensor-to-sensor phase coherences and thus PLV/PLI interactions.
%   - change "run_source_modeling.m" to only include sensors (h.anatomy.sens.good_sensors) in lead fields when calculating inverse solutions
%   - added "Spatial Temporal" Source searching to find peaks voxels in space with lowest residual variance (best fit - after waveform normalization) to the true source waveforms
%   - added text data points to error plots on "Source Modeling Tab"
%   - source waveform error is now "% Residual Variance" = sum of squares of differences between peak source and true source waveforms
%   - changed source plot colors to peak "Hits" = true source colors and all other false positives as transparent orange
%   - when selecting from "Peaks Found" listbox, sources are highlighted in purple and waveforms become solid
%   - main figure on startup --> adjusts screensize to fit within monitor's resolution
%   - created "Aanlyses" panel on "Source Modeling" tab
%   - added "Seeded Functional Connectivity" Analyses with plotting TFRs and 3D connectivity grap; includes Surrogate Statistics
%   - added ability to inverse project sensor noise to Seed and comparison locations to use as baseline statistics similar to surrogate statistics for PLV, null distributions are
%       based on reshaped matrix of within the active interval of the inverse solution
%   - fixed "plot normalized waves" so that "true source waves" are also normalized just liked the "peak source waves" so that the "% Residual Variance" is more accurate measure.
%   - fixed dB calculation for wavelet TFR in plot_SimSignals.m --> dB = 10*( log10(active) - log10(baseline) );
%
%
% Bug or Need-to-Fix report
%   - Vector inverse solution --> need to add option for selecting maximal dipole orientation <-- currently set to rms of 3 dipole orientations to yield final maps and waveforms
%   - allow use to define "dip_thresh" search distances for spatiotemporal searches --> line 44 in "sm_spatiotemp_mapping.m"
%   - add in tooltips for all buttons, edit boxes, etc.
%   - Time-Frequency Analysis --> change to or add in "ft_freqanalysis.m" because it has a lot more functionality and evidence for use with M/EEG
%
% Modules to add
%   - inverse project sens_noise to true source dipoles using leadfield weights to yield "source_noise" --> then project sens_noise using inv_soln weights and compare to source_noise
%   	Aslo, project sim_data.noise_wav to sensor space to yield "source_noise"
%

version_num = 2.1;

clear h
global h
warning('off','all');
% opengl software;

%% Finding varargin that has anatomy variables called from Brainstorm
for vi = 1:length(varargin)
    if isfield(varargin{1},'subj_MriFile')
        bst = varargin{vi};
        h.bst_called_flag = 1;   % SimMEEG was called from Brainstorm
    end
end

% btn = questdlg(sprintf('SimMEEG Software\n\nCopyright (C) 2020  \nDr. Anthony Thomas Herdman and The University of British Columbia\n\nThis program comes with ABSOLUTELY NO WARRANTY\n\nThis program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation\n\nPlease kindly acknowledge Dr. Anthony Herdman at UBC when presenting work from SimMEEG and/or redistributing this software.\n\nThis software is not licensed for clinical use or training.\n\n Liability: The developers, distrubtors, and licensors of this product are not liable in any form or in any way for losses that you, the licensee, may incur from use or re-distribution of this software.\n\n Click "Accept" if you agree to these licensing terms, otherwise click "Decline"'),'License Agreement','Accept','Decline','Accept');
h.license_flag = 0;
sm_popup_license_terms;

if h.license_flag==1
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Adding paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h.pwd=pwd;
    %% User-defined or Brainstorm paths
    if exist('bst','var')   % Brainstorm preset paths - if SimMEEG called directly from BrainStorm
        if isfield(bst,'subj_anat_dir') && isfield(bst,'subj_data_dir') && isfield(bst,'FieldTrip_dir')
            h.data_dir = bst.subj_data_dir;
            h.anat_path = bst.subj_anat_dir;
            h.FieldTrip_dir = bst.FieldTrip_dir;
        else                                                            % User-defined paths
            h.data_dir = uigetdir(h.pwd,'Set Data Directory');
            h.anat_path = uigetdir(h.pwd,'Set Anatomy Directory');
            h.FieldTrip_dir = uigetdir(h.pwd,'Open Field Trip Directory');
        end
    else                                                                % User-defined paths
        h.data_dir = uigetdir(h.pwd,'Set Data Directory');
        h.anat_path = uigetdir(h.pwd,'Set Anatomy Directory');
        h.FieldTrip_dir = uigetdir(h.pwd,'Open Field Trip Directory');
        h.bst_subj_data_dir ='';
    end
    %% Adding FieldTrip paths
    cd(h.FieldTrip_dir); %C:\BRANELab\Software\BRANE_Lab_v3\fieldtrip-20200519\;
    ft_defaults;
    cd(h.pwd);
    % Field Trip's external directory with filtfilt messes up BRANE Lab's filter_data_new.m function so need to remove this directory for filtering function
    h.FieldTrip_dir_external = genpath([ h.FieldTrip_dir,'\external\']);
    rmpath(h.FieldTrip_dir_external);
    
    %% Initializing UserData
    h.UserData.bkg_clr = [1 1 1];    % background color
    % h.cfg.study.source_locs = [91 48 78; 35 109 78; 151 109 78]; % [primary vis, Left Aud, Right Aud];
    h.cfg.study.source_locs = [86 26 77; 37 114 77; 150 110 79]; % [primary vis, Left Aud, Right Aud];
    % h.tfr_ROI(1)=''; h.tfr_ROI(2)=''; h.tfr_ROI(3)=''; % time-frequency ROIs for each source
    h.src_clr=[0 .6 0; 0 0 1; 1 0 0];     %[.8 .8 1; 1 .8 .8; .8 1 .8];
    h.src_clr2=[.7 .9 .7; .8 .8 1; 1 .8 .8];
    h.plv_clr = [.7 0 .9; 1 0 1; 1 .6 0];
    h.plv_clr2 = [.7 0 .9; 1 0 1; 1 .6 0]*.5;
    % h.src_clr=[0 0 0; .8 0 .8; 1 .5 0];     %[.8 .8 1; 1 .8 .8; .8 1 .8];
    % h.src_clr2=[.6 .6 .6; 1 .8 1; 1 .8 .7];
    
    h.chan_clr = [1 0 .5]*.7; %[1 0 1];
    h.sens_clr = [0 0 0];
    h.brain_clr = [1 1 1]*.5; % [.6 .8 1];  %
    h.scalp_clr = [178 116 62]/255; %[1 .8 .8]; %[.6 .8 1];
    
    h.trial_wave_clr = [1 1 1]*.8;
    
    % h.evk_clr=[0 0 1; 1 0 0; 0 1 0];     % Evoked not determined for this GUI. Evoked waves are calculated while running the simulation.
    h.freq_vals=1:200;
    h.caxis_power=[-100 100]; h.caxis_PLV=[-100 100]; h.caxis_PLI=[-100 100];
    h.num_sig_freqs=0;
    h.current_source_num=1; h.current_tfr_roi_idx=1;
    h.topo_lat_samp =[];
    h.inv_soln = [];
    h.current_3D_thresh = .5;
    h.current_filt_freqs = [0 0 0 0];
    h.current_inv_soln = [];
    h.plv_contrasts = [1 2; 1 3; 2 3];
    h.monte_carlo_flag = 0;
    h.new_study_flag  = 0;
    h.start_flag = 1;
    % h.bst_subj_anat_dir = h.data_dir;
    % h.bst_subj_data_dir = h.data_dir;
    h.brainstorm_db = h.data_dir;
    h.bst_subj_anat_dir = h.data_dir;
    h.anatomy = struct('mri','','sens_eeg','','sens_meg','','mesh_volumes',[],'headmodel_eeg_vol','','headmodel_meg_vol','','headmodel_eeg_cortex','','headmodel_meg_cortex','','leadfield_eeg_vol','','leadfield_eeg_cortex','','leadfield_meg_vol','','leadfield_meg_cortex','','sens','','headmodel','','leadfield','','mesh_cortex','');
    h.real_datadir = ' ';
    h.font_size_warndlg = 11;
    h.find_spatiotemp_peaks_flag = 0; % (0) regular peak search inv_soln map (1) spatiotemporal search for peaks in inv_soln source waveforms
    h.false_positive_FaceAlpha = .5;   % FaceAlpha transparency for plotting source locations for false positives
    h.false_positive_lineAlpha = .25;   % lineAlpha transparency for plotting source waves for false positives
    h.cfg.source.TFR_results = [];   % True Source Time-freq results
    h.current_inv_tfr_time_point = []; % time-freq point on inverse-source modleing Tab's TFR plots
    h.current_inv_true_fc_data = []; %
    h.current_inv_peak_fc_data = []; %
    h.fieldtrip_tfr_data = [];
    h.current_fieldtrip_tfr_data = [];
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h.main_fig = figure; set(h.main_fig,'Units','normalized','Position',[.05 .05 .85 .85],'Name',sprintf('SimMEEG v%.1f',version_num),'NumberTitle','off');
    addToolbarExplorationButtons(h.main_fig) % Adds zoom, rotate, ... buttons to figure toolbar
    h.main_fig.CloseRequestFcn = @closereq_SimMEEG;
    
    %% %%%%% Panel "Study Parameters" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h.study_panel = uipanel(h.main_fig,'Title','Study Parameters','FontSize',12,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[.6 0 1],...
        'Units','normalize','Position',[.005 .865 .86 .13],'FontWeight','bold');
    h.txt_length   = .12;
    h.ui_length    = .05;
    h.ui_height    = .2;
    h.txt_clmn_pos = .01:.185:.99; %.01:.205:.99;
    h.ui_clmn_pos  = .135:.185:.99; %.135:.205:.99;
    h.txt_row_pos  = fliplr(.1:h.ui_height+.025:1-h.ui_height+.02);
    %% srate
    h.edit_srate_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(1) h.txt_row_pos(1) h.txt_length h.ui_height],'FontSize',10,'HorizontalAlignment','left','String','Sample Rate');
    h.edit_srate = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Tag','edit_srate','Position',[h.ui_clmn_pos(1) h.txt_row_pos(1) h.ui_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','center','String','256','Callback',@update_study_cfg);
    %% Trial Duration
    h.edit_dur_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(1) h.txt_row_pos(2) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Trial Duration (sec)');
    h.edit_dur = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Tag','edit_dur','Position',[h.ui_clmn_pos(1) h.txt_row_pos(2) h.ui_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','center','String','-1 1','Callback',@update_study_cfg);
    % update_study_cfg;
    % h.num_samps_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
    %     'Position',[h.edit_dur.Position(1)+h.edit_dur.Position(3)+.01 1-h.study_panel.Position(4)-.01 .1 h.study_panel.Position(4)],...
    %     'FontSize',10,'HorizontalAlignment','left','String',sprintf('Samples = %.f',h.cfg.study.num_samps));
    %% Numbers of Trials
    h.edit_num_trials_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(1) h.txt_row_pos(3) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Number Trials');
    h.edit_num_trials = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(1) h.txt_row_pos(3) h.ui_length h.ui_height],...
        'FontSize',8,'HorizontalAlignment','center','String','90','Callback',@update_study_cfg);
    %% Signal-to-Noise Ratio (SNR)
    h.edit_sens_SNR_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(1) h.txt_row_pos(4) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','SNR (dB)');
    h.edit_sens_SNR = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(1) h.txt_row_pos(4) h.ui_length h.ui_height],...
        'FontSize',8,'HorizontalAlignment','center','String','-7','Value',1);
    %% Noise Amplitude Percent
    h.edit_noise_amp_perc_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(2) h.txt_row_pos(2) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Noise Ampltidue (%)');
    h.edit_noise_amp_perc = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(2) h.txt_row_pos(2) h.ui_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','center','String','10','Callback',@update_study_cfg);
    %% Noise Frequncy Band
    h.edit_noise_freqs_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(2) h.txt_row_pos(1) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Noise Band (Hz)');
    h.edit_noise_freqs = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(2) h.txt_row_pos(1) h.ui_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','center','String','1 100','Callback',@update_study_cfg);
    %% Noise Type
    h.edit_noise_flag_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(2) h.txt_row_pos(3) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Noise Type');
    h.edit_noise_flag = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[h.ui_clmn_pos(2)-.05 h.txt_row_pos(3) h.ui_length+.05 h.ui_height],...
        'FontSize',8,'HorizontalAlignment','center','String',{'Broad Band' 'Narrow Band' 'Notched Band' 'Pink/Brown'},'Value',2,'Callback',@update_study_cfg);
    h.edit_pink_noise_slope_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(2) h.txt_row_pos(4)-.05 h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Pink Noise Slope');
    h.edit_pink_noise_slope = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(2) h.txt_row_pos(4)-.05 h.ui_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','center','String','1 100','Callback',@update_study_cfg);
    %% Baseline Interval
    h.edit_base_int_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(3) h.txt_row_pos(1) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Baseline Interval');
    h.edit_base_int = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(3) h.txt_row_pos(1) h.ui_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','center','String','-0.5 0','Callback',@update_study_cfg);
    %% Post-Event Interval
    h.edit_poststim_int_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(3) h.txt_row_pos(2) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Post-Event Interval','Tooltip','(seconds) Used for calculating SNR');
    h.edit_poststim_int = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(3) h.txt_row_pos(2) h.ui_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','center','String','0 0.5','Callback',@update_study_cfg,'Tooltip','(seconds) Used for calculating SNR');
    %% PLV Target Threshold
    h.edit_plv_thresh_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(3) h.txt_row_pos(3) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','PLV Target Threshold');
    h.edit_plv_thresh = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(3) h.txt_row_pos(3) h.ui_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','center','String','0.05','Callback',@update_study_cfg);
    %% plots Simulated Data pop-ups
    % h.radio_plot_sim_flag = uicontrol(h.study_panel,'Style','radiobutton', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
    %     'Position',[.89 h.txt_row_pos(1) h.txt_length h.ui_height],...
    %     'FontSize',10,'HorizontalAlignment','left','String','Plot Sim Data');
    %% Num PLV Perms
    h.edit_max_perm_plv_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(3) h.txt_row_pos(4) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','PLV Permutations');
    h.edit_max_perm_plv = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[h.ui_clmn_pos(3) h.txt_row_pos(4) h.ui_length h.ui_height],...
        'FontSize',8,'HorizontalAlignment','center','String',{'7' '8' '9' '10'},'Value',3,'Callback',@update_study_cfg);
    %% Plot Time Interval
    h.edit_plot_time_int_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(4) h.txt_row_pos(1) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Plot Interval (sec)');
    h.edit_plot_time_int = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(4) h.txt_row_pos(1) h.ui_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','center','String','-1 1','Callback',@update_study_cfg);
    %% Plot Freq Interval
    h.edit_plot_freq_int_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(4) h.txt_row_pos(2) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Plot Frequencies (Hz)');
    h.edit_plot_freq_int = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(4) h.txt_row_pos(2) h.ui_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','center','String','1 60','Callback',@update_study_cfg);
    %% Plot Color Axis
    h.edit_plot_caxis_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(4) h.txt_row_pos(3) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Plot Power Scale (%)');
    h.edit_plot_caxis = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(4) h.txt_row_pos(3) h.ui_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','center','String','-100 100','Callback',@update_study_cfg);
    %% Edit Study Name
    h.edit_study_name_txt = uicontrol(h.study_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.txt_clmn_pos(4) h.txt_row_pos(4) h.txt_length h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','Study Name Prefix:','Tooltip','Study Name to be used as a prefix for saving Datasets.');
    h.edit_study_name = uicontrol(h.study_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[h.ui_clmn_pos(4) h.txt_row_pos(4) h.ui_length*4 h.ui_height],...
        'FontSize',10,'HorizontalAlignment','left','String','SimMEEG','Tooltip','Study Name to be used as a prefix for saving Datasets.');
    %% Save Study
    h.btn_save_study = uicontrol(h.study_panel,'BackgroundColor',[1 .8 .8],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.91 .79 .075 h.ui_height],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','Save Study','Callback',@save_study);
    %% Load Study
    h.btn_load_study = uicontrol(h.study_panel,'BackgroundColor',[.8 1 .8],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.91 .57 .075 h.ui_height],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','Load Study','Callback',@load_study);
    %% New Study
    h.btn_new_study = uicontrol(h.study_panel,'BackgroundColor',[.8 .8 1],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.91 .35 .075 h.ui_height],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','New Study','Callback',@new_study);
    %% Study Directory
    h.btn_datadir = uicontrol(h.study_panel,'BackgroundColor',[1 1 1]*.9,'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.91 .13 .075 h.ui_height],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','Set Data Dir','Callback',@get_datadir);
    
    
    %% Initializing h.cfg.study data
    h.cfg.study.srate=str2num(h.edit_srate.String);
    h.cfg.study.dur=str2num(h.edit_dur.String); % (sec) start and end times for whole trial
    h.cfg.study.lat_sim=[h.cfg.study.dur(1):1/h.cfg.study.srate:h.cfg.study.dur(2)-(1/h.cfg.study.srate)]; % latency of each trial
    h.cfg.study.num_samps=length(h.cfg.study.lat_sim);
    h.cfg.study.num_trials = str2num(h.edit_num_trials.String); %{h.edit_num_trials.Value});
    h.cfg.study.max_perm_plv = str2num(h.edit_max_perm_plv.String{h.edit_max_perm_plv.Value}); %h.cfg.study.num_trials/10; % integer of cfg.study.num_trials/10 = maximum number of permutations to search for PLV_trials (must be <10 or memory will fail on most computers) .
    h.cfg.study.noise_flag=h.edit_noise_flag.Value;  % Whitening noise to be added to each sources
    h.cfg.study.noise_amp_perc=str2num(h.edit_noise_amp_perc.String); % percent of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
    h.cfg.study.noise_freqs=str2num(h.edit_noise_freqs.String); % [1 100] Frequency (Hz) of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
    h.cfg.study.plv_thresh=str2num(h.edit_plv_thresh.String);   % stoppping criterion when search for best PLV/PLI matched to sig_PLV_targets, sig_PLI_targets, etc. (e.g., 0.05).
    h.cfg.study.plot_sim_flag=0; %h.radio_plot_sim_flag.Value;  % plots evoked, trial, and PLV results. Note: This can take time because of filtering.
    h.cfg.study.plot_time_int=str2num(h.edit_plot_time_int.String);   % time interval to plot
    h.cfg.study.plot_freq_int=str2num(h.edit_plot_freq_int.String); % frequencies for calculating and plotting PLV/PLI
    h.cfg.study.base_int=str2num(h.edit_base_int.String);    % base line interval for plotting
    h.cfg.study.poststim_int=str2num(h.edit_poststim_int.String);    % base line interval for plotting
    h.cfg.study.act_samps = (round( (h.cfg.study.poststim_int(1)-h.cfg.study.lat_sim(1))*h.cfg.study.srate ):round( (h.cfg.study.poststim_int(2)-h.cfg.study.lat_sim(1))*h.cfg.study.srate )) +1;
    h.cfg.study.ctrl_samps = (round( (h.cfg.study.base_int(1)-h.cfg.study.lat_sim(1))*h.cfg.study.srate ):round( (h.cfg.study.base_int(2)-h.cfg.study.lat_sim(1))*h.cfg.study.srate )) +1;
    h.cfg.study.plot_time_vals=h.cfg.study.plot_time_int(1):1/h.cfg.study.srate:h.cfg.study.plot_time_int(2);
    h.cfg.study.plot_freq_vals=h.cfg.study.plot_freq_int(1):h.cfg.study.plot_freq_int(2);
    h.cfg.study.pink_noise_slope = normrnd(1.8,.1,1); % Arbitrarily set to be randomized with overall mean 1.8  % seemed to best match the slope from a couple of resting-state data from LetterAll study
    if h.cfg.study.pink_noise_slope<1; h.cfg.study.pink_noise_slope=1; elseif h.cfg.study.pink_noise_slope>2; h.cfg.study.pink_noise_slope=2; end
    h.edit_pink_noise_slope.String = sprintf('%.1f',h.cfg.study.pink_noise_slope);
    %% Update Study Information
    update_study_cfg;
    
    %% %%%%% Panel "Simulate Sources" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h.run_simulation_panel = uipanel(h.main_fig,'Title','Simulate Sources','FontSize',12,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[.6 0 1],...
        'Units','normalize','Position',[.87 .865 .125 .13],'Visible','on');
    %% Btn Update cfg
    h.btn_update_cfg = uicontrol(h.run_simulation_panel,'BackgroundColor',[1 .8 .8],'ForegroundColor',[0 0 0],'Style','pushbutton','Units','normalize',...
        'Position',[.05 .78 .9 .2],'Visible','on','UserData',1,...
        'FontSize',10,'HorizontalAlignment','center','String','Update Params','Callback',@update_cfg);
    %% Btn Simulate Source data
    h.btn_run_sim = uicontrol(h.run_simulation_panel,'BackgroundColor',[.8 1 .8],'ForegroundColor',[0 0 0],'Style','pushbutton','Units','normalize',...
        'Position',[.05 .78-.25 .9 .2],'Visible','on','UserData',1,...
        'FontSize',10,'HorizontalAlignment','center','String','Sim Source Data','Callback',@run_sim);
    %% Btn "Calc Source PLV" Data
    h.btn_plot_sim_data = uicontrol(h.run_simulation_panel,'BackgroundColor',[.8 .8 1],'ForegroundColor',[0 0 0],'Style','pushbutton','Units','normalize',...
        'Position',[.05 .78-.5 .9 .2],'Visible','on','UserData',1,...
        'FontSize',10,'HorizontalAlignment','center','String','Calc Source PLV','Callback',@plot_SimSignals);
    %% Edit: Wavelet Time-Bandwitdh parameter
    % delete(h.edit_wavelet_TB_txt); delete(h.edit_wavelet_TB);
    h.edit_wavelet_TB_txt = uicontrol(h.run_simulation_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.1 .02 .55 .2],...
        'FontSize',10,'HorizontalAlignment','left','String','Wavelet TB');
    h.edit_wavelet_TB = uicontrol(h.run_simulation_panel,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[.65 .02 .2 .2],'Tooltip',sprintf('The larger the time-bandwidth parameter, \nthe more spread out the wavelet is in time\nand narrower the wavelet is in frequency.'),...
        'FontSize',9,'HorizontalAlignment','center','String','30');
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tabs for Source, PLV, and PLI parameters and Simulate M/EEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Creating Tabs
    h.tabgrp = uitabgroup(h.main_fig,'Units','normalize','Position',[.005 .005 .995 .86],'Visible','on');
    % h.tab_anatomy = uitab(h.tabgrp,'Title','Anatomy','BackgroundColor',[1 1 1],'Foregroundcolor',[0 0 0]);
    h.tab_power = uitab(h.tabgrp,'Title','Event-Related Power','BackgroundColor',[1 1 1],'Foregroundcolor',[0 .6 0]);
    h.tab_PLV = uitab(h.tabgrp,'Title','Phase-Locking Value (PLV)','BackgroundColor',[1 1 1],'Foregroundcolor',h.plv_clr(1,:));
    h.tab_PLI = uitab(h.tabgrp,'Title','Phase-Lag Index (PLI)','BackgroundColor',[1 1 1],'Foregroundcolor',h.plv_clr(2,:));
    h.tab_PAC = uitab(h.tabgrp,'Title','Phase-Amplitude Coupling (PAC)','BackgroundColor',[1 1 1],'Foregroundcolor',[.2 .75 0]);
    h.tab_sim_meeg = uitab(h.tabgrp,'Title','Simulate M/EEG','BackgroundColor',[1 1 1],'Foregroundcolor',[0 .6 .6]);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tab "Source Modeling" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    create_tab_source_modeling;
    h.tab_monte_carlo = uitab(h.tabgrp,'Title','Monte Carlo','BackgroundColor',[1 1 1],'Foregroundcolor',[0 .2 .6]);
    % h.tab_preprocess = uitab(h.tabgrp,'Title','Data Preprocess','BackgroundColor',[1 1 1],'Foregroundcolor',[.3 .6 .5]);
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tab "Event-Related Power" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%% SOURCE TFR Paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% standard parameters for text and edit boxes
    h.tfr_xpos=[.05 .35 .65];
    h.tfr_ypos=.55; h.tfr_pos_size=[.265 .3];
    h.sig_ypos=.3; h.sig_pos_size=[.265 .2];
    h.prepost_ypos=.05;
    h.source_edit_length = .035; h.source_edit_height = .03; h.source_edit_ypos  = [.92 .88];
    %% Signal Window Type
    h.menu_sig_win_type_txt = uicontrol(h.tab_power,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 .965 .1 .03],'FontSize',10,'HorizontalAlignment','right','String','Windowing Type');
    h.menu_sig_win_type = uicontrol(h.tab_power,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[.11 .968 .08 .03],'FontSize',8,'HorizontalAlignment','center','String',{'Hanning' 'Gaussian' 'Triangular','Blackmann'},'Value',3,'Callback',@update_study_cfg);
    h.menu_sig_win_type.Enable='on'; % currently setting it so that it only does triangular windowing
    %% Button Select TFR ROI
    h.btn_ROI_source1 = uicontrol(h.tab_power,'BackgroundColor',[1 1 1]*.9,'ForegroundColor',h.src_clr(1,:),'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.menu_sig_win_type.Position([1 3]))+.01 .965 .075 .035],'Visible','on','UserData',1,...
        'FontSize',10,'HorizontalAlignment','center','String','Select ROI','Callback',@bl_rbbox_gca_v2);
    %% Menu: ROI
    h.menu_tfr_roi_idx_txt = uicontrol(h.tab_power,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.btn_ROI_source1.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .035 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','ROI #');
    h.menu_tfr_roi_idx = uicontrol(h.tab_power,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(h.menu_tfr_roi_idx_txt.Position([1 3]))+.005 .968 .065 .03],'FontSize',8,'HorizontalAlignment','center','String',{''},'Value',1,'Callback',{@set_current_ROI,0});    % menu listing TFR ROIs
    %% Menu: Source #
    h.menu_source_idx_txt = uicontrol(h.tab_power,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_tfr_roi_idx.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .04 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Source');
    h.menu_source_idx = uicontrol(h.tab_power,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(h.menu_source_idx_txt.Position([1 3]))+.005 .968 .035 .03],'FontSize',8,'HorizontalAlignment','center','String',{'1' '2' '3'},'Value',1,'Callback',{@set_current_source,0});    % menu listing TFR ROIs
    %% Edit: TFR ROI limits
    % Frequency
    h.edit_tfr_roi_ylim_txt = uicontrol(h.tab_power,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_source_idx.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .105 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Frequency Range');
    h.edit_tfr_roi_ylim = uicontrol(h.tab_power,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_ylim_txt.Position([1 3]))+.005 h.menu_sig_win_type_txt.Position(2) .05 .035],...
        'Tag','TFR   Frequency Range','FontSize',10,'HorizontalAlignment','center','String','0 0','Callback',@set_tfr_roi_cfg);
    % Latency
    h.edit_tfr_roi_xlim_txt = uicontrol(h.tab_power,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_ylim.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .095 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Latency Range');
    h.edit_tfr_roi_xlim = uicontrol(h.tab_power,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_xlim_txt.Position([1 3]))+.005 h.menu_sig_win_type_txt.Position(2) .075 .035],...
        'Tag','TFR   Latency Range','FontSize',10,'HorizontalAlignment','center','String','0 0','Callback',@set_tfr_roi_cfg);
    % Rise Time
    h.edit_tfr_roi_risetime_txt = uicontrol(h.tab_power,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_xlim.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .095 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Rise/Fall Time');
    h.edit_tfr_roi_risetime = uicontrol(h.tab_power,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_risetime_txt.Position([1 3]))+.005 h.menu_sig_win_type_txt.Position(2) .045 .035],...
        'Tag','TFR   Latency Range','FontSize',10,'HorizontalAlignment','center','String','0 0','Callback',@set_tfr_roi_cfg);
    %% Edit: Source Power & Evoked Percent
    sig_tag = {'Sig 1 ' 'Sig 2 ' 'Sig 3 '};
    prepost_tag = {'Pre 1 ' 'Pre 2 ' 'Pre 3 '};
    box_type = {'Power' 'PLV' 'PLI' 'Latency Range' 'Frequency Range' };
    
    for v=1:3
        %% Signal Power Percent
        h.edit_sig_power_perc_txt(v) = uicontrol(h.tab_power,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(v,:),'Units','normalize',...
            'Position',[h.tfr_xpos(v) h.source_edit_ypos(1) .105 .0275],...
            'FontSize',10,'HorizontalAlignment','right','String','Power (%) Signal ');
        h.edit_sig_power_perc(v) = uicontrol(h.tab_power,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(v,:),'Style','edit','Units','normalize',...
            'Position',[sum(h.edit_sig_power_perc_txt(v).Position([1 3]))+.005 h.source_edit_ypos(1) h.source_edit_length h.source_edit_height],...
            'Tag',[sig_tag{v} box_type{1}],'FontSize',10,'HorizontalAlignment','center','String','100','Callback',@set_tfr_roi_cfg);
        %% PrePost power Percent
        h.edit_prepost_power_perc_txt(v) = uicontrol(h.tab_power,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(v,:),'Units','normalize',...
            'Position',[sum(h.edit_sig_power_perc(v).Position([1 3]))+.005 h.source_edit_ypos(1) .05 .0275],...
            'FontSize',10,'HorizontalAlignment','right','String','Prepost');
        h.edit_prepost_power_perc(v) = uicontrol(h.tab_power,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(v,:),'Style','edit','Units','normalize',...
            'Position',[sum(h.edit_prepost_power_perc_txt(v).Position([1 3]))+.005 h.source_edit_ypos(1) h.source_edit_length h.source_edit_height],...
            'Tag',[prepost_tag{v} box_type{1}],'FontSize',10,'HorizontalAlignment','center','String','100','Callback',@set_tfr_roi_cfg);
        
        %% Signal Evoked Percent
        h.edit_sig_evoked_perc_txt(v) = uicontrol(h.tab_power,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(v,:),'Units','normalize',...
            'Position',[h.tfr_xpos(v) h.source_edit_ypos(2) .105 .0275],...
            'FontSize',10,'HorizontalAlignment','right','String','Evoked (%) Signal ');
        h.edit_sig_evoked_perc(v) = uicontrol(h.tab_power,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(v,:),'Style','edit','Units','normalize',...
            'Position',[sum(h.edit_sig_evoked_perc_txt(v).Position([1 3]))+.005 h.source_edit_ypos(2) h.source_edit_length h.source_edit_height],...
            'Tag',[sig_tag{v} box_type{1}],'FontSize',10,'HorizontalAlignment','center','String','50','Callback',@set_tfr_roi_cfg);
        %% PrePost Evoked Percent
        h.edit_prepost_evoked_perc_txt(v) = uicontrol(h.tab_power,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(v,:),'Units','normalize',...
            'Position',[sum(h.edit_sig_evoked_perc(v).Position([1 3]))+.005 h.source_edit_ypos(2) .05 .0275],...
            'FontSize',10,'HorizontalAlignment','right','String','Prepost');
        h.edit_prepost_evoked_perc(v) = uicontrol(h.tab_power,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(v,:),'Style','edit','Units','normalize',...
            'Position',[sum(h.edit_prepost_evoked_perc_txt(v).Position([1 3]))+.005 h.source_edit_ypos(2) h.source_edit_length h.source_edit_height],...
            'Tag',[prepost_tag{v} box_type{1}],'FontSize',10,'HorizontalAlignment','center','String','0','Callback',@set_tfr_roi_cfg);
    end
    %% Source 1 Power Axes
    h.ax_power(1)=axes(h.tab_power,'Position',[h.tfr_xpos(1) h.tfr_ypos h.tfr_pos_size]);
    h.tfr_data.s1=surf(h.ax_power(1),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
    view(h.ax_power(1),0,90); shading(h.ax_power(1),'flat'); axis(h.ax_power(1),'tight'); axis(h.ax_power(1),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
    hold on; h.tfr_data.tfr_zero_line(1) = plot([0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title('Source 1 Power','Color',h.src_clr(1,:)); colormap(jet(255)); caxis(h.caxis_power);
    % disableDefaultInteractivity(h.ax_power(1)); %h.tfr_data.s1.ButtonDownFcn=@bl_rbbox_gca_v2;
    h.ax_power(1).YLabel.String='Frequency (Hz)';
    %% Source 2 Power Axes
    h.ax_power(2)=axes(h.tab_power,'Position',[h.tfr_xpos(2) h.tfr_ypos h.tfr_pos_size]);
    h.tfr_data.s2=surf(h.ax_power(2),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
    view(h.ax_power(2),0,90); shading(h.ax_power(2),'flat'); axis(h.ax_power(2),'tight'); axis(h.ax_power(2),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
    hold on; h.tfr_data.tfr_zero_line(2) = plot([0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title('Source 2 Power','Color',h.src_clr(2,:)); colormap(jet(255)); caxis(h.caxis_power);
    % disableDefaultInteractivity(h.ax_power(2)); %h.tfr_data.s1.ButtonDownFcn=@bl_rbbox_gca_v2;
    %% Source 3 Power Axes
    h.ax_power(3)=axes(h.tab_power,'Position',[h.tfr_xpos(3) h.tfr_ypos h.tfr_pos_size]);
    h.tfr_data.s3=surf(h.ax_power(3),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
    view(h.ax_power(3),0,90); shading(h.ax_power(3),'flat'); axis(h.ax_power(3),'tight'); axis(h.ax_power(3),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
    hold on; h.tfr_data.tfr_zero_line(3) = plot([0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title('Source 3 Power','Color',h.src_clr(3,:)); colormap(jet(255)); caxis(h.caxis_power);
    % disableDefaultInteractivity(h.ax_power(3)); %h.tfr_data.s1.ButtonDownFcn=@bl_rbbox_gca_v2;
    pos=h.ax_power(3).Position; cb=colorbar;
    h.ax_power(3).Position=pos;
    %% Source Power Wave Axes
    h.ax_sig_waves(1)=axes(h.tab_power,'Position',[h.tfr_xpos(1) h.sig_ypos h.sig_pos_size]); box on; title('Source 1 Target Signal Waves','Color',h.src_clr(1,:));
    h.ax_sig_waves(2)=axes(h.tab_power,'Position',[h.tfr_xpos(2) h.sig_ypos h.sig_pos_size]); box on; title('Source 2 Target Signal Waves','Color',h.src_clr(2,:));
    h.ax_sig_waves(3)=axes(h.tab_power,'Position',[h.tfr_xpos(3) h.sig_ypos h.sig_pos_size]); box on; title('Source 3 Target Signal Waves','Color',h.src_clr(3,:));
    h.ax_prepost_waves(1)=axes(h.tab_power,'Position',[h.tfr_xpos(1) h.prepost_ypos h.sig_pos_size]); box on; title('Source 1 Target PrePost Waves','Color',h.src_clr(1,:));
    h.ax_prepost_waves(2)=axes(h.tab_power,'Position',[h.tfr_xpos(2) h.prepost_ypos h.sig_pos_size]); box on; title('Source 2 Target PrePost Waves','Color',h.src_clr(2,:));
    h.ax_prepost_waves(3)=axes(h.tab_power,'Position',[h.tfr_xpos(3) h.prepost_ypos h.sig_pos_size]); box on; title('Source 3 Target PrePost Waves','Color',h.src_clr(3,:));
    %% point position when dragging cursor point
    h.pt_txt = uicontrol(h.tab_power,'Style','text', 'BackgroundColor',h.tab_power.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[.35 0 .3 .025],'FontSize',10,'HorizontalAlignment','center','String','Lat = 0 ms   Amp = 0');
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tab "Phase-Locking Value (PLV)" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%% SOURCE PLV Paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Source 1 vs 2 PLV Axes
    h.ax_PLV(1)=axes(h.tab_PLV,'Position',[h.tfr_xpos(1) h.tfr_ypos h.tfr_pos_size]);
    h.plv_data.s1=surf(h.ax_PLV(1),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
    view(h.ax_PLV(1),0,90); shading(h.ax_PLV(1),'flat'); axis(h.ax_PLV(1),'tight'); axis(h.ax_PLV(1),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
    hold on; h.plv_data.tfr_zero_line(1) = plot([0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title('Source 1 vs 2 Target PLV','Color',h.plv_clr(1,:)); colormap(jet(255)); caxis(h.caxis_PLV);
    % disableDefaultInteractivity(h.ax_PLV(1)); %h.plv_data.s1.ButtonDownFcn=@bl_rbbox_gca_v2;
    h.ax_PLV(1).YLabel.String='Frequency (Hz)';
    %% Source 1 vs 3 PLV Axes
    h.ax_PLV(2)=axes(h.tab_PLV,'Position',[h.tfr_xpos(2) h.tfr_ypos h.tfr_pos_size]);
    h.plv_data.s2=surf(h.ax_PLV(2),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
    view(h.ax_PLV(2),0,90); shading(h.ax_PLV(2),'flat'); axis(h.ax_PLV(2),'tight'); axis(h.ax_PLV(2),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
    hold on; h.plv_data.tfr_zero_line(2) = plot([0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title('Source 1 vs 3 Target PLV','Color',h.plv_clr(2,:)); colormap(jet(255)); caxis(h.caxis_PLV);
    % disableDefaultInteractivity(h.ax_PLV(2)); %h.plv_data.s1.ButtonDownFcn=@bl_rbbox_gca_v2;
    %% Source 2 vs 3 PLV Axes
    h.ax_PLV(3)=axes(h.tab_PLV,'Position',[h.tfr_xpos(3) h.tfr_ypos h.tfr_pos_size]);
    h.plv_data.s3=surf(h.ax_PLV(3),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
    view(h.ax_PLV(3),0,90); shading(h.ax_PLV(3),'flat'); axis(h.ax_PLV(3),'tight'); axis(h.ax_PLV(3),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
    hold on; h.plv_data.tfr_zero_line(3) = plot([0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title('Source 2 vs 3 Target PLV','Color',h.plv_clr(3,:)); colormap(jet(255)); caxis(h.caxis_PLV);
    % disableDefaultInteractivity(h.ax_PLV(3)); %h.plv_data.s1.ButtonDownFcn=@bl_rbbox_gca_v2;
    pos=h.ax_PLV(3).Position; cb=colorbar;
    h.ax_PLV(3).Position=pos;
    %% Source PLV Wave Axes
    h.ax_sig_plv(1)=axes(h.tab_PLV,'Position',[h.tfr_xpos(1) h.sig_ypos h.sig_pos_size]); box on; title('Source 1 vs 2 Target Signal PLV Waves','Color',h.plv_clr(1,:));
    h.ax_sig_plv(2)=axes(h.tab_PLV,'Position',[h.tfr_xpos(2) h.sig_ypos h.sig_pos_size]); box on; title('Source 1 vs 3 Target Signal PLV Waves','Color',h.plv_clr(2,:));
    h.ax_sig_plv(3)=axes(h.tab_PLV,'Position',[h.tfr_xpos(3) h.sig_ypos h.sig_pos_size]); box on; title('Source 2 vs 3 Target Signal PLV Waves','Color',h.plv_clr(3,:));
    h.ax_prepost_plv(1)=axes(h.tab_PLV,'Position',[h.tfr_xpos(1) h.prepost_ypos h.sig_pos_size]); box on; title('Source 1 vs 2 Target PrePost PLV Waves','Color',h.plv_clr(1,:));
    h.ax_prepost_plv(2)=axes(h.tab_PLV,'Position',[h.tfr_xpos(2) h.prepost_ypos h.sig_pos_size]); box on; title('Source 1 vs 3 Target PrePost PLV Waves','Color',h.plv_clr(2,:));
    h.ax_prepost_plv(3)=axes(h.tab_PLV,'Position',[h.tfr_xpos(3) h.prepost_ypos h.sig_pos_size]); box on; title('Source 2 vs 3 Target PrePost PLV Waves','Color',h.plv_clr(3,:));
    
    for a=1:3;  h.ax_PLV(a).Toolbar.Visible = 'off'; h.ax_sig_plv(a).Toolbar.Visible = 'off'; h.ax_prepost_plv(a).Toolbar.Visible = 'off'; end
    %% Menu: ROI
    h.menu_tfr_roi_idx_txt_PLV = uicontrol(h.tab_PLV,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.btn_ROI_source1.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .035 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','ROI #');
    h.menu_tfr_roi_idx_PLV = uicontrol(h.tab_PLV,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(h.menu_tfr_roi_idx_txt_PLV.Position([1 3]))+.005 .968 .065 .03],'FontSize',8,'HorizontalAlignment','center','String',{''},'Value',1,'Callback',{@set_current_ROI,1});    % menu listing TFR ROIs
    %% Menu: Source #
    h.menu_source_idx_txt_PLV = uicontrol(h.tab_PLV,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_tfr_roi_idx_PLV.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .04 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Source');
    h.menu_source_idx_PLV = uicontrol(h.tab_PLV,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(h.menu_source_idx_txt_PLV.Position([1 3]))+.005 .968 .035 .03],'FontSize',8,'HorizontalAlignment','center','String',{'1' '2' '3'},'Value',1,'Callback',{@set_current_source,1});    % menu listing TFR ROIs
    %% ROI limits
    % Frequency
    h.edit_tfr_roi_ylim_txt_PLV = uicontrol(h.tab_PLV,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_source_idx_PLV.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .105 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Frequency Range');
    h.edit_tfr_roi_ylim_PLV = uicontrol(h.tab_PLV,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_ylim_txt_PLV.Position([1 3]))+.005 h.menu_sig_win_type_txt.Position(2) .05 .035],...
        'Tag','PLV   Frequency Range','FontSize',10,'HorizontalAlignment','center','String','0 0','Callback',@set_tfr_roi_cfg);
    % Latency
    h.edit_tfr_roi_xlim_txt_PLV = uicontrol(h.tab_PLV,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_ylim_PLV.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .095 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Latency Range');
    h.edit_tfr_roi_xlim_PLV = uicontrol(h.tab_PLV,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_xlim_txt_PLV.Position([1 3]))+.005 h.menu_sig_win_type_txt.Position(2) .075 .035],...
        'Tag','PLV   Latency Range','FontSize',10,'HorizontalAlignment','center','String','0 0','Callback',@set_tfr_roi_cfg);
    % Rise Time
    h.edit_tfr_roi_risetime_txt_PLV = uicontrol(h.tab_PLV,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_xlim_PLV.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .095 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Rise/Fall Time');
    h.edit_tfr_roi_risetime_PLV = uicontrol(h.tab_PLV,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_risetime_txt_PLV.Position([1 3]))+.005 h.menu_sig_win_type_txt.Position(2) .045 .035],...
        'Tag','PLV   Latency Range','FontSize',10,'HorizontalAlignment','center','String','0 0','Callback',@set_tfr_roi_cfg);
    %% Edit Boxes
    for v=1:3
        %% Signal Phase-Locking Value
        h.edit_sig_phase_locking_txt(v) = uicontrol(h.tab_PLV,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.plv_clr(v,:),'Units','normalize',...
            'Position',[h.tfr_xpos(v) h.source_edit_ypos(1) .12 .0275],...
            'FontSize',10,'HorizontalAlignment','left','String','Phase Locking: Signal ');
        h.edit_sig_phase_locking(v) = uicontrol(h.tab_PLV,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.plv_clr(v,:),'Style','edit','Units','normalize',...
            'Position',[sum(h.edit_sig_evoked_perc_txt(v).Position([1 3]))+.01 h.source_edit_ypos(1) h.source_edit_length h.source_edit_height],...
            'Tag',[sig_tag{v} box_type{2}],'FontSize',10,'HorizontalAlignment','center','String','0.4','Callback',@set_tfr_roi_cfg);
        %% PrePost Phase-Locking Value
        h.edit_prepost_phase_locking_txt(v) = uicontrol(h.tab_PLV,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.plv_clr(v,:),'Units','normalize',...
            'Position',[sum(h.edit_sig_evoked_perc(v).Position([1 3]))+.005 h.source_edit_ypos(1) .05 .0275],...
            'FontSize',10,'HorizontalAlignment','right','String','Prepost');
        h.edit_prepost_phase_locking(v) = uicontrol(h.tab_PLV,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.plv_clr(v,:),'Style','edit','Units','normalize',...
            'Position',[sum(h.edit_prepost_evoked_perc_txt(v).Position([1 3]))+.005 h.source_edit_ypos(1) h.source_edit_length h.source_edit_height],...
            'Tag',[prepost_tag{v} box_type{2}],'FontSize',10,'HorizontalAlignment','center','String','0.0','Callback',@set_tfr_roi_cfg);
        
    end
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tab "Phase-Lag Index (PLI)" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%% SOURCE PLI Paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Signal Phase Lag
    %% Menu: ROI
    h.menu_tfr_roi_idx_txt_PLI = uicontrol(h.tab_PLI,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.btn_ROI_source1.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .035 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','ROI #');
    h.menu_tfr_roi_idx_PLI = uicontrol(h.tab_PLI,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(h.menu_tfr_roi_idx_txt_PLI.Position([1 3]))+.005 .968 .065 .03],'FontSize',8,'HorizontalAlignment','center','String',{''},'Value',1,'Callback',{@set_current_ROI,2});    % menu listing TFR ROIs
    %% Menu: Source #
    h.menu_source_idx_txt_PLI = uicontrol(h.tab_PLI,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_tfr_roi_idx_PLI.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .04 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Source');
    h.menu_source_idx_PLI = uicontrol(h.tab_PLI,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(h.menu_source_idx_txt_PLI.Position([1 3]))+.005 .968 .035 .03],'FontSize',8,'HorizontalAlignment','center','String',{'1' '2' '3'},'Value',1,'Callback',{@set_current_source,2});    % menu listing TFR ROIs
    %% Edit: ROI limits
    % Frequency
    h.edit_tfr_roi_ylim_txt_PLI = uicontrol(h.tab_PLI,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_source_idx_PLI.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .105 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Frequency Range');
    h.edit_tfr_roi_ylim_PLI = uicontrol(h.tab_PLI,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_ylim_txt_PLI.Position([1 3]))+.005 h.menu_sig_win_type_txt.Position(2) .05 .035],...
        'Tag','PLI   Frequency Range','FontSize',10,'HorizontalAlignment','center','String','0 0','Callback',@set_tfr_roi_cfg);
    % Latency
    h.edit_tfr_roi_xlim_txt_PLI = uicontrol(h.tab_PLI,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_ylim_PLI.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .095 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Latency Range');
    h.edit_tfr_roi_xlim_PLI = uicontrol(h.tab_PLI,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_xlim_txt_PLI.Position([1 3]))+.005 h.menu_sig_win_type_txt.Position(2) .075 .035],...
        'Tag','PLV   Latency Range','FontSize',10,'HorizontalAlignment','center','String','0 0','Callback',@set_tfr_roi_cfg);
    % Rise Time
    h.edit_tfr_roi_risetime_txt_PLI = uicontrol(h.tab_PLI,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_xlim_PLI.Position([1 3]))+.01 h.menu_sig_win_type_txt.Position(2) .095 .0275],...
        'FontSize',10,'HorizontalAlignment','right','String','Rise/Fall Time');
    h.edit_tfr_roi_risetime_PLI = uicontrol(h.tab_PLI,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_tfr_roi_risetime_txt_PLI.Position([1 3]))+.005 h.menu_sig_win_type_txt.Position(2) .045 .035],...
        'Tag','PLV   Latency Range','FontSize',10,'HorizontalAlignment','center','String','0 0','Callback',@set_tfr_roi_cfg);
    %% Edit Boxes
    for v=1:3
        h.edit_sig_phase_lag_txt(v) = uicontrol(h.tab_PLI,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.plv_clr(v,:),'Units','normalize',...
            'Position',[h.tfr_xpos(v) h.source_edit_ypos(1)-.025 .11 .05],...
            'FontSize',10,'HorizontalAlignment','right','String','Phase Lag: Signal ');
        h.edit_sig_phase_lag(v) = uicontrol(h.tab_PLI,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.plv_clr(v,:),'Style','edit','Units','normalize',...
            'Position',[sum(h.edit_sig_evoked_perc_txt(v).Position([1 3]))+.01 h.source_edit_ypos(1) h.source_edit_length h.source_edit_height],...
            'Tag',[sig_tag{v} box_type{3}],'FontSize',10,'HorizontalAlignment','center','String','0.4','Callback',@set_tfr_roi_cfg);
        %% PrePost Evoked Percent
        h.edit_prepost_phase_lag_txt(v) = uicontrol(h.tab_PLI,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.plv_clr(v,:),'Units','normalize',...
            'Position',[sum(h.edit_sig_evoked_perc(v).Position([1 3]))+.005 h.source_edit_ypos(1) .05 .0275],...
            'FontSize',10,'HorizontalAlignment','right','String','Prepost');
        h.edit_prepost_phase_lag(v) = uicontrol(h.tab_PLI,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.plv_clr(v,:),'Style','edit','Units','normalize',...
            'Position',[sum(h.edit_prepost_evoked_perc_txt(v).Position([1 3]))+.005 h.source_edit_ypos(1) h.source_edit_length h.source_edit_height],...
            'Tag',[prepost_tag{v} box_type{3}],'FontSize',10,'HorizontalAlignment','center','String','0.0','Callback',@set_tfr_roi_cfg);
        
        %% Starting Phase relative to trial onset - determines the relative correlations among sources 1, 2, and 3
        %% Signal Starting Phase
        h.edit_sig_phase_start_txt(v) = uicontrol(h.tab_PLI,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.plv_clr(v,:),'Units','normalize',...
            'Position',[h.tfr_xpos(v) h.source_edit_ypos(2) .11 .0275],...
            'FontSize',10,'HorizontalAlignment','left','String','Starting Phase (deg): Signal ');
        h.edit_sig_phase_start(v) = uicontrol(h.tab_PLI,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.plv_clr(v,:),'Style','edit','Units','normalize',...
            'Position',[sum(h.edit_sig_phase_start_txt(v).Position([1 3]))+.005 h.source_edit_ypos(2) h.source_edit_length h.source_edit_height],'Tooltip',sprintf('inter-source correlation = cos( deg2rad (phase difference between sources) ) )'),...
            'Tag',[sig_tag{v} box_type{2}],'FontSize',10,'HorizontalAlignment','center','String','0','Callback',@update_cfg);
        %% PrePost Starting Phase
        h.edit_prepost_phase_start_txt(v) = uicontrol(h.tab_PLI,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.plv_clr(v,:),'Units','normalize',...
            'Position',[sum(h.edit_sig_phase_start(v).Position([1 3]))+.005 h.source_edit_ypos(2) .04 .0275],...
            'FontSize',10,'HorizontalAlignment','left','String','Prepost');
        h.edit_prepost_phase_start(v) = uicontrol(h.tab_PLI,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.plv_clr(v,:),'Style','edit','Units','normalize',...
            'Position',[sum( h.edit_prepost_phase_start_txt(v).Position([1 3]))+.005 h.source_edit_ypos(2) h.source_edit_length h.source_edit_height],'Tooltip',sprintf('inter-source correlation = cos( deg2rad (phase difference between sources) ) )'),...
            'Tag',[prepost_tag{v} box_type{2}],'FontSize',10,'HorizontalAlignment','center','String','0','Callback',@update_cfg);
        
    end
    %% Source 1 vs 2 PLI Axes
    h.ax_PLI(1)=axes(h.tab_PLI,'Position',[h.tfr_xpos(1) h.tfr_ypos h.tfr_pos_size]);
    h.pli_data.s1=surf(h.ax_PLI(1),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
    view(h.ax_PLI(1),0,90); shading(h.ax_PLI(1),'flat'); axis(h.ax_PLI(1),'tight'); axis(h.ax_PLI(1),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
    hold on; h.pli_data.tfr_zero_line(1) = plot([0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title('Source 1 vs 2 Target PLI','Color',h.plv_clr(1,:)); colormap(jet(255)); caxis(h.caxis_PLI);
    % disableDefaultInteractivity(h.ax_PLI(1)); %h.pli_data.s1.ButtonDownFcn=@bl_rbbox_gca_v2;
    h.ax_PLI(1).YLabel.String='Frequency (Hz)';
    %% Source 1 vs 3PLI Axes
    h.ax_PLI(2)=axes(h.tab_PLI,'Position',[h.tfr_xpos(2) h.tfr_ypos h.tfr_pos_size]);
    h.pli_data.s2=surf(h.ax_PLI(2),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
    view(h.ax_PLI(2),0,90); shading(h.ax_PLI(2),'flat'); axis(h.ax_PLI(2),'tight'); axis(h.ax_PLI(2),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
    hold on; h.pli_data.tfr_zero_line(2) = plot([0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title('Source 1 vs 3 Target PLI','Color',h.plv_clr(2,:)); colormap(jet(255)); caxis(h.caxis_PLI);
    % disableDefaultInteractivity(h.ax_PLI(2)); %h.pli_data.s1.ButtonDownFcn=@bl_rbbox_gca_v2;
    %% Source 2 vs 3 PLI Axes
    h.ax_PLI(3)=axes(h.tab_PLI,'Position',[h.tfr_xpos(3) h.tfr_ypos h.tfr_pos_size]);
    h.pli_data.s3=surf(h.ax_PLI(3),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
    view(h.ax_PLI(3),0,90); shading(h.ax_PLI(3),'flat'); axis(h.ax_PLI(3),'tight'); axis(h.ax_PLI(3),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
    hold on; h.pli_data.tfr_zero_line(3) = plot([0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title('Source 2 vs 3 PLI','Color',h.plv_clr(3,:)); colormap(jet(255)); caxis(h.caxis_PLI);
    % disableDefaultInteractivity(h.ax_PLI(3)); %h.pli_data.s1.ButtonDownFcn=@bl_rbbox_gca_v2;
    pos=h.ax_PLI(3).Position; cb=colorbar;
    h.ax_PLI(3).Position=pos;
    %% Source PLI Wave Axes
    h.ax_sig_pli(1)=axes(h.tab_PLI,'Position',[h.tfr_xpos(1) h.sig_ypos h.sig_pos_size]); box on; title('Source 1 vs 2 Target Signal PLI Waves','Color',h.plv_clr(1,:));
    h.ax_sig_pli(2)=axes(h.tab_PLI,'Position',[h.tfr_xpos(2) h.sig_ypos h.sig_pos_size]); box on; title('Source 1 vs 3 Target Signal PLI Waves','Color',h.plv_clr(2,:));
    h.ax_sig_pli(3)=axes(h.tab_PLI,'Position',[h.tfr_xpos(3) h.sig_ypos h.sig_pos_size]); box on; title('Source 2 vs 3 Target Signal PLI Waves','Color',h.plv_clr(3,:));
    h.ax_prepost_pli(1)=axes(h.tab_PLI,'Position',[h.tfr_xpos(1) h.prepost_ypos h.sig_pos_size]); box on; title('Source 1 vs 2 Target PrePost PLI Waves','Color',h.plv_clr(1,:));
    h.ax_prepost_pli(2)=axes(h.tab_PLI,'Position',[h.tfr_xpos(2) h.prepost_ypos h.sig_pos_size]); box on; title('Source 1 vs 3 Target PrePost PLI Waves','Color',h.plv_clr(2,:));
    h.ax_prepost_pli(3)=axes(h.tab_PLI,'Position',[h.tfr_xpos(3) h.prepost_ypos h.sig_pos_size]); box on; title('Source 2 vs 3 Target PrePost PLI Waves','Color',h.plv_clr(3,:));
    
    for a=1:3;  h.ax_PLI(a).Toolbar.Visible = 'off'; h.ax_sig_pli(a).Toolbar.Visible = 'off'; h.ax_prepost_pli(a).Toolbar.Visible = 'off'; end
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tab "Phase-Amplitude Coupling (PAC)" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%% Panel "Phase-Amplitude Coupling (PAC) Parameters" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delete(h.panel_PAC_params)
    h.panel_PAC_params = uipanel(h.tab_PAC,'Title','Phase-Amplitude Coupling (PAC) Parameters','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.01 .64 .98 .355],'Visible','on');
    % index (1) = Source 1 v Source 2; (2)=[1 3];  (2)=[2 1];  (2)=[2 3];  (2)=[3 ];  (2)=[3 2];
    %% Button Update PAC waves
    % delete(h.btn_update_PAC_waves)
    h.btn_update_PAC_waves = uicontrol(h.panel_PAC_params,'BackgroundColor',[1 1 1]*.9,'ForegroundColor',h.src_clr(1,:),'Style','pushbutton','Units','normalize',...
        'Position',[.87 .91 .12 .085],'Visible','on','UserData',1,...
        'FontSize',10,'HorizontalAlignment','center','String','Update PAC Waves','Callback',{@btn_update_PAC_waves,0});
    %% Menu & Text for Source contrasts for PAC parameters
    h.PAC_source_contrasts = [1 2; 1 3; 2 1; 2 3; 3 1; 3 2];
    xsize = .055; xpos = [];
    ysize = .075; ypos = .8-[ysize:ysize+.05:.95];
    xsize_edit = .075;
    for a=1:6
        %% Signal PAC
        % text for Carrier Source
        h.menu_PAC_sig_carrier_source_txt(a) = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,1),:),...
            'Units','normalize','Position',[.01 ypos(a) xsize ysize],...
            'FontSize',10,'HorizontalAlignment','left','String',sprintf('Source %.f',h.PAC_source_contrasts(a,1)));
        % Menu for Carrier Source Frequency
        h.menu_PAC_sig_carrier_source_freq(a) = uicontrol(h.panel_PAC_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
            'Position',[sum(h.menu_PAC_sig_carrier_source_txt(a).Position([1 3]))+.0025 h.menu_PAC_sig_carrier_source_txt(a).Position(2) xsize_edit ysize],'FontSize',8,...
            'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,1),:),'HorizontalAlignment','center','String',{''},'Value',1,'Callback',{@update_PAC,0});    % menu listing PAC carrier frequencies
        % text for Modulator Source
        h.menu_PAC_sig_modulator_source_txt(a) = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,2),:),...
            'Units','normalize','Position',[sum(h.menu_PAC_sig_carrier_source_freq(a).Position([1 3]))+.01 h.menu_PAC_sig_carrier_source_freq(a).Position(2) xsize ysize],...
            'FontSize',10,'HorizontalAlignment','left','String',sprintf('Source %.f',h.PAC_source_contrasts(a,2)));
        % Menu for Modulator Source Frequency
        h.menu_PAC_sig_modulator_source_freq(a) = uicontrol(h.panel_PAC_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
            'Position',[sum(h.menu_PAC_sig_modulator_source_txt(a).Position([1 3]))+.0025 h.menu_PAC_sig_modulator_source_txt(a).Position(2) xsize_edit ysize],'FontSize',8,...
            'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,2),:),'HorizontalAlignment','center','String',{''},'Value',1,'Callback',{@update_PAC,0});    % menu listing PAC carrier frequencies
        
        % Signal Depth of modulator
        h.edit_PAC_modulator_sig_depth(a) = uicontrol(h.panel_PAC_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
            'Position',[sum(h.menu_PAC_sig_modulator_source_freq(1).Position([1 3]))+.02 h.menu_PAC_sig_modulator_source_freq(a).Position(2) xsize_edit ysize],'FontSize',8,...
            'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,2),:),'HorizontalAlignment','center','String',{'0'},'Value',1,'Callback',{@update_PAC,0});    % menu listing PAC carrier frequencies
        h.edit_PAC_modulator_sig_depth_range(a) = uicontrol(h.panel_PAC_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
            'Position',[sum(h.edit_PAC_modulator_sig_depth(a).Position([1 3]))+.02  h.edit_PAC_modulator_sig_depth(a).Position(2) xsize_edit ysize],'FontSize',8,...
            'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,2),:),'HorizontalAlignment','center','String',{'0'},'Value',1,'Callback',{@update_PAC,0});    % menu listing PAC carrier frequencies
        
        %% Prepost PAC
        % text for Carrier Source
        h.menu_PAC_prepost_carrier_source_txt(a) = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,1),:),...
            'Units','normalize','Position',[sum(h.edit_PAC_modulator_sig_depth_range(a).Position([1 3]))+.05 ypos(a) xsize ysize],...
            'FontSize',10,'HorizontalAlignment','left','String',sprintf('Source %.f',h.PAC_source_contrasts(a,1)));
        % Menu for Carrier Source Frequency
        h.menu_PAC_prepost_carrier_source_freq(a) = uicontrol(h.panel_PAC_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
            'Position',[sum(h.menu_PAC_prepost_carrier_source_txt(a).Position([1 3]))+.0025 h.menu_PAC_prepost_carrier_source_txt(a).Position(2) xsize_edit ysize],'FontSize',8,...
            'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,1),:),'HorizontalAlignment','center','String',{''},'Value',1,'Callback',{@update_PAC,0});    % menu listing PAC carrier frequencies
        % text for Modulator Source
        h.menu_PAC_prepost_modulator_source_txt(a) = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,2),:),...
            'Units','normalize','Position',[sum(h.menu_PAC_prepost_carrier_source_freq(a).Position([1 3]))+.01 h.menu_PAC_prepost_carrier_source_freq(a).Position(2) xsize ysize],...
            'FontSize',10,'HorizontalAlignment','left','String',sprintf('Source %.f',h.PAC_source_contrasts(a,2)));
        % Menu for Modulator Source Frequency
        h.menu_PAC_prepost_modulator_source_freq(a) = uicontrol(h.panel_PAC_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
            'Position',[sum(h.menu_PAC_prepost_modulator_source_txt(a).Position([1 3]))+.0025 h.menu_PAC_prepost_modulator_source_txt(a).Position(2) xsize_edit ysize],'FontSize',8,...
            'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,2),:),'HorizontalAlignment','center','String',{''},'Value',1,'Callback',{@update_PAC,0});    % menu listing PAC carrier frequencies
        
        % Prepost Depth of modulator
        h.edit_PAC_modulator_prepost_depth(a) = uicontrol(h.panel_PAC_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
            'Position',[sum(h.menu_PAC_prepost_modulator_source_freq(a).Position([1 3]))+.02 h.menu_PAC_prepost_modulator_source_freq(a).Position(2) xsize_edit ysize],'FontSize',8,...
            'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,2),:),'HorizontalAlignment','center','String',{'0'},'Value',1,'Callback',{@update_PAC,0});    % menu listing PAC carrier frequencies
        h.edit_PAC_modulator_prepost_depth_range(a) = uicontrol(h.panel_PAC_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
            'Position',[sum(h.edit_PAC_modulator_prepost_depth(a).Position([1 3]))+.02  h.edit_PAC_modulator_prepost_depth(a).Position(2) xsize_edit ysize],'FontSize',8,...
            'Foregroundcolor',h.src_clr(h.PAC_source_contrasts(a,2),:),'HorizontalAlignment','center','String',{'0'},'Value',1,'Callback',{@update_PAC,0});    % menu listing PAC carrier frequencies
        
    end
    %% Text along top
    %% Signal heading texts
    % text for Signal Carrier and Modulator
    h.PAC_sig_carrier_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.menu_PAC_sig_carrier_source_txt(1).Position(1) .82 xsize ysize*2],...
        'FontSize',10,'HorizontalAlignment','left','String','Signal Carrier');
    h.PAC_sig_modulator_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.menu_PAC_sig_modulator_source_txt(1).Position(1) h.PAC_sig_carrier_txt.Position(2) xsize ysize*2],...
        'FontSize',10,'HorizontalAlignment','left','String','Signal Modulator');
    % text "Frequency" for Signal Carrier and Modulator
    h.PAC_sig_carrier_freq_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.menu_PAC_sig_carrier_source_freq(1).Position(1) h.PAC_sig_carrier_txt.Position(2) xsize ysize],...
        'FontSize',10,'HorizontalAlignment','left','String','Frequency');
    h.PAC_sig_modulator_freq_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.menu_PAC_sig_modulator_source_freq(1).Position(1) h.PAC_sig_carrier_txt.Position(2) xsize ysize],...
        'FontSize',10,'HorizontalAlignment','left','String','Frequency');
    
    % Text for Depth & Range
    h.PAC_sig_depth_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.edit_PAC_modulator_sig_depth(1).Position(1) h.PAC_sig_carrier_txt.Position(2) .1 ysize],...
        'FontSize',10,'HorizontalAlignment','left','String','Modulator Depth');
    h.PAC_sig_range_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.edit_PAC_modulator_sig_depth_range(1).Position(1) h.PAC_sig_carrier_txt.Position(2) h.PAC_sig_depth_txt.Position(3) ysize],...
        'FontSize',10,'HorizontalAlignment','left','String','Depth Range');
    %% Prepost heading texts
    % text for Signal Carrier and Modulator
    h.PAC_prepost_carrier_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.menu_PAC_prepost_carrier_source_txt(1).Position(1) h.PAC_sig_carrier_txt.Position(2) xsize ysize*2],...
        'FontSize',10,'HorizontalAlignment','left','String','PrePost Carrier');
    h.PAC_prepost_modulator_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.menu_PAC_prepost_modulator_source_txt(1).Position(1) h.PAC_sig_carrier_txt.Position(2) xsize ysize*2],...
        'FontSize',10,'HorizontalAlignment','left','String','PrePost Modulator');
    h.PAC_prepost_depth_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.edit_PAC_modulator_prepost_depth(1).Position(1) h.PAC_sig_carrier_txt.Position(2) .1 ysize],...
        'FontSize',10,'HorizontalAlignment','left','String','Prepost Depth');
    h.PAC_prepost_range_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.edit_PAC_modulator_prepost_depth_range(1).Position(1) h.PAC_sig_carrier_txt.Position(2) h.PAC_prepost_depth_txt.Position(3) ysize],...
        'FontSize',10,'HorizontalAlignment','left','String','Depth Range');
    % text "Frequency" for Signal Carrier and Modulator
    h.PAC_prepost_carrier_freq_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.menu_PAC_prepost_carrier_source_freq(1).Position(1) h.PAC_sig_carrier_txt.Position(2) xsize ysize],...
        'FontSize',10,'HorizontalAlignment','left','String','Frequency');
    h.PAC_prepost_modulator_freq_txt = uicontrol(h.panel_PAC_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.menu_PAC_prepost_modulator_source_freq(1).Position(1) h.PAC_prepost_carrier_txt.Position(2) xsize ysize],...
        'FontSize',10,'HorizontalAlignment','left','String','Frequency');
    
    %% %%%% Panel "Target PAC Waves" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delete(h.panel_PAC_waves)
    h.panel_PAC_waves = uipanel(h.tab_PAC,'Title','Target PAC waves','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.01 .01 .98 .63],'Visible','on');
    %% Axes PAC waves
    % delete(h.axes_PAC_waves_sig); delete(h.axes_PAC_waves_prepost)
    xsize = .22; ysize = .25;
    xpos1 = [.02 .27 .02 .27 .02 .27];
    xpos2 = [.02 .27 .02 .27 .02 .27]+.5;
    ypos =  [.7  .7 .375  .375 .05 .05];
    
    for a=1:6
        h.axes_PAC_waves_sig(a)=axes(h.panel_PAC_waves,'Position',[xpos1(a) ypos(a) xsize ysize]); box on; title(sprintf('Signal PAC: Source %.f mod by %.f ',h.PAC_source_contrasts(a,:)));
        h.axes_PAC_waves_prepost(a)=axes(h.panel_PAC_waves,'Position',[xpos2(a) ypos(a) xsize ysize]); box on; title(sprintf('PrePost PAC: Source %.f mod by %.f ',h.PAC_source_contrasts(a,:)));
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tab "Simulate M/EEG" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%% Panel "Forward Model: Anatomy & Source Locations" on Tab "Simulate M/EEG"  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h.panel_anatomy = uipanel(h.tab_sim_meeg,'Title','Forward Model: Anatomy & Source Locations','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.01 .01 .49 .98],'Visible','on');
    
    %% Load SimMEEG Anatomy.mat - must specific struct (see SimMEEG.m program for details)
    h.btn_load_default_anatomy = uicontrol(h.panel_anatomy,'BackgroundColor',[.2 .2 1]*1,'ForegroundColor',[1 1 1]*1,'Style','pushbutton','Units','normalize',...
        'Position',[.01 .96 .225 .035],...
        'FontSize',10,'HorizontalAlignment','center','String','SimMEEG Anatomy','Callback',@load_default_anatomy);
    %% Load Brainstorm Anatomy from subj anat directory
    h.btn_load_bst_anatomy = uicontrol(h.panel_anatomy,'BackgroundColor',[1 1 1]*.4,'ForegroundColor',[1 1 1]*1,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_load_default_anatomy.Position([1 3]))+.01 .96 .25 .035],...
        'FontSize',10,'HorizontalAlignment','center','String','BrainStorm Anatomy','Callback',@load_bst_anatomy);
    %% Load FieldTrip Anatomy.mat - must have mri, sens_eeg, sens_MEG, etc. saved inside the *.mat file
    h.btn_load_ft_anatomy = uicontrol(h.panel_anatomy,'BackgroundColor',[.1 .6 .1]*1,'ForegroundColor',[1 1 1]*1,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_load_bst_anatomy.Position([1 3]))+.01 .96 .25 .035],...
        'FontSize',10,'HorizontalAlignment','center','String','FieldTrip Anatomy','Callback',@sm_load_FieldTrip_Anatomy);
    %% New Anatomy by clearing h.anatomy
    h.btn_clear_SimMEEG_anatomy = uicontrol(h.panel_anatomy,'BackgroundColor',[1 .9 .8],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_load_ft_anatomy.Position([1 3]))+.05 .96 .185 .035],...
        'FontSize',10,'HorizontalAlignment','center','String','Clear Anatomy','Callback',@sm_clear_SimMEEG_Anatomy);
    %% Save FieldTrip Anatomy.mat - must have mri, sens_eeg, sens_MEG, etc. saved inside the *.mat file
    h.btn_save_SimMEEG_anatomy = uicontrol(h.panel_anatomy,'BackgroundColor',[.8 .9 1],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_load_ft_anatomy.Position([1 3]))+.05 h.btn_clear_SimMEEG_anatomy.Position(2)-.07 .185 .035],...
        'FontSize',10,'HorizontalAlignment','center','String','Save Anatomy','Callback',@sm_save_SimMEEG_Anatomy);
    %% Anatomy File Name
    h.mri_path_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 .925 .2 .03],'FontSize',10,'HorizontalAlignment','right','String','Anatomy Folder:');
    h.mri_path = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.22 h.mri_path_txt.Position(2) .8-sum(h.btn_load_default_anatomy.Position([1 3]))+.01 .03],'FontSize',10,'HorizontalAlignment','left','String','Anatomy Folder');
    h.mri_file_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 h.mri_path_txt.Position(2)-.0325 .2 .03],'FontSize',10,'HorizontalAlignment','right','String','Anatomy File:');
    h.anatomy_file_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.22 h.mri_file_txt.Position(2) .8-sum(h.btn_load_default_anatomy.Position([1 3]))+.01 .03],...
        'FontSize',10,'HorizontalAlignment','left','String','Default MRI','Callback',@set_anat_file);
    %% Load MRI - created in Field Trip format
    h.btn_load_anatomy = uicontrol(h.panel_anatomy,'BackgroundColor',[1 .7 .7],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.01 h.anatomy_file_txt.Position(2)-.035 .2 .03],...
        'FontSize',10,'HorizontalAlignment','center','String','Load MRI','Callback',@sm_load_mri);
    h.mri_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.btn_load_anatomy.Position([1 3]))+.01 h.btn_load_anatomy.Position(2) .8-sum(h.btn_load_anatomy.Position([1 3]))+.01 .03],...
        'FontSize',10,'HorizontalAlignment','left','String','MRI Info:____________________________________');
    %% Load Sensors - created in Field Trip format
    h.btn_load_sensors = uicontrol(h.panel_anatomy,'BackgroundColor',[1 .8 .4],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.01 h.btn_load_anatomy.Position(2)-.035 .2 .03],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','Load Sensors','Callback',@sm_load_sensors);
    h.sensors_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.btn_load_sensors.Position([1 3]))+.01 h.btn_load_sensors.Position(2) .8-sum(h.btn_load_sensors.Position([1 3]))+.01 .03],...
        'FontSize',10,'HorizontalAlignment','left','String','Sensor Info:____________________________________');
    h.menu_sens_type = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[1-.21 h.sensors_txt.Position(2)+.01 .18 .025],...
        'FontSize',8,'HorizontalAlignment','center','String',{'MEG' 'EEG'},'Value',1,'Callback',@menu_head_model_CallBack);
    %% Load HeadModel - created in Field Trip format
    h.btn_load_headmodel = uicontrol(h.panel_anatomy,'BackgroundColor',[.8 1 .8],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.01 h.btn_load_sensors.Position(2)-.035 .2 .03],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','Create HeadModel','Callback',@sm_create_headmodel);
    h.headmodel_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.btn_load_headmodel.Position([1 3]))+.01 h.btn_load_headmodel.Position(2) .8-sum(h.btn_load_headmodel.Position([1 3]))+.01 .03],...
        'FontSize',10,'HorizontalAlignment','left','String','Head Model Info:________________________________');
    %% Load Lead Field - created in Field Trip format
    h.btn_load_leadfield = uicontrol(h.panel_anatomy,'BackgroundColor',[.8 .8 1],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.01 h.btn_load_headmodel.Position(2)-.035 .2 .03],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','Create LeadField','Callback',@sm_create_leadfield);
    h.leadfield_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.btn_load_leadfield.Position([1 3]))+.01 h.btn_load_leadfield.Position(2) .8-sum(h.btn_load_leadfield.Position([1 3]))+.01 .03],...
        'FontSize',10,'HorizontalAlignment','left','String','Lead Field Info:_________________________________');
    %% Menu Head Model - Whole Head or Cortical Surface
    % h.menu_head_model_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
    %     'Position',[1-.38 h.btn_load_leadfield.Position(2)-.035 .16 .02],...
    %     'FontSize',10,'HorizontalAlignment','left','String','Head Model');
    h.menu_head_model = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[1-.21 h.btn_load_leadfield.Position(2)+.01 .18 .025],...
        'FontSize',8,'HorizontalAlignment','center','String',{'Volume' 'Cortical Surface'},'Value',1,'Callback',@menu_head_model_CallBack);
    %% Select Source Locations using bl_set_source_locs.m
    h.btn_select_source_locs = uicontrol(h.panel_anatomy,'BackgroundColor',[1 1 1]*.9,'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.01 h.btn_load_leadfield.Position(2)-.07 .2 .03],'Enable','inactive',...
        'FontSize',10,'HorizontalAlignment','center','String','Source Locations','Callback',@select_source_locs);
    %% Txt for locs ori and amp
    h.source_amp_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.btn_select_source_locs.Position([1 3]))+.01 h.btn_select_source_locs.Position(2) .25 .03],...
        'FontSize',10,'HorizontalAlignment','left','String','Amplitudes (nA):');
    h.source_loc_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.btn_select_source_locs.Position([1 3]))+.01 h.source_amp_txt.Position(2)-.035 .25 .03],...
        'FontSize',10,'HorizontalAlignment','left','String','Locations X,Y,Z (mm):');
    h.source_ori_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.btn_select_source_locs.Position([1 3]))+.01 h.source_loc_txt.Position(2)-.035 .265 .03],...
        'FontSize',10,'HorizontalAlignment','left','String','Orientations Az,El (deg):');
    txt_wdth=(.9-sum(h.source_ori_txt.Position([1 3])))/3;  % width for source text following "Orientations" text
    h.source1_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(1,:),'Units','normalize',...
        'Position',[sum(h.source_ori_txt.Position([1 3]))+.01 h.btn_load_leadfield.Position(2)-.035 txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','Source 1');
    h.source2_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(2,:),'Units','normalize',...
        'Position',[sum(h.source1_txt.Position([1 3]))+.01 h.source1_txt.Position(2) txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','Source 2');
    h.source3_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(3,:),'Units','normalize',...
        'Position',[sum(h.source2_txt.Position([1 3]))+.01 h.source1_txt.Position(2) txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','Source 3');
    %% Source Locations
    h.edit_source_locs(1) = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(1,:),'Style','edit','Units','normalize',...
        'Tag','edit_source1_ori','Position',[h.source1_txt.Position(1) h.source1_txt.Position(2)-.065 txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','0 0 0','Callback',@edit_source_CallBack);
    h.edit_source_locs(2) = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(2,:),'Style','edit','Units','normalize',...
        'Tag','edit_source2_ori','Position',[h.source2_txt.Position(1) h.source1_txt.Position(2)-.065 txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','0 0 0','Callback',@edit_source_CallBack);
    h.edit_source_locs(3) = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(3,:),'Style','edit','Units','normalize',...
        'Tag','edit_source3_ori','Position',[h.source3_txt.Position(1) h.source1_txt.Position(2)-.065 txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','0 0 0','Callback',@edit_source_CallBack);
    %% Source Amplitudes
    h.edit_source_amp(1) = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(1,:),'Style','edit','Units','normalize',...
        'Tag','edit_source1_amp','Position',[h.source1_txt.Position(1) h.source1_txt.Position(2)-.03 txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','30','Callback',@edit_source_CallBack);
    h.edit_source_amp(2) = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(2,:),'Style','edit','Units','normalize',...
        'Tag','edit_source2_amp','Position',[h.source2_txt.Position(1) h.source1_txt.Position(2)-.03 txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','30','Callback',@edit_source_CallBack);
    h.edit_source_amp(3) = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(3,:),'Style','edit','Units','normalize',...
        'Tag','edit_source3_amp','Position',[h.source3_txt.Position(1) h.source1_txt.Position(2)-.03 txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','30','Callback',@edit_source_CallBack);
    %% Source Orientations
    h.edit_source_ori(1) = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(1,:),'Style','edit','Units','normalize',...
        'Tag','edit_source1_ori','Position',[h.source1_txt.Position(1) h.source1_txt.Position(2)-.1 txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','180 0','Callback',@edit_source_CallBack);
    h.edit_source_ori(2) = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(2,:),'Style','edit','Units','normalize',...
        'Tag','edit_source2_ori','Position',[h.source2_txt.Position(1) h.source1_txt.Position(2)-.1 txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','245 305','Callback',@edit_source_CallBack);
    h.edit_source_ori(3) = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(3,:),'Style','edit','Units','normalize',...
        'Tag','edit_source3_ori','Position',[h.source3_txt.Position(1) h.source1_txt.Position(2)-.1 txt_wdth .03],...
        'FontSize',10,'HorizontalAlignment','center','String','295 305','Callback',@edit_source_CallBack);
    %% Menu Orientations normal to surface
    h.menu_ori_normal_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.edit_source_ori(1).Position(1) h.edit_source_ori(1).Position(2)-.04 .25 .025],...
        'FontSize',10,'HorizontalAlignment','left','String','Orientation Constraint');
    h.menu_ori_normal = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(h.menu_ori_normal_txt.Position([1 3]))+.01 h.menu_ori_normal_txt.Position(2)+.01 .18 .025],...
        'FontSize',8,'HorizontalAlignment','center','String',{'Random' 'Cortical Surface'},'Value',1,'Callback',@menu_ori_normal_CallBack);
    %% Plot 3D MRI
    h.btn_plot_3D_mri = uicontrol(h.panel_anatomy,'BackgroundColor',[.9 .9 1],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.01  h.menu_ori_normal_txt.Position(2)-.035 .2 .03],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','plot 3D','Callback',@plot_3D_mri);
    %% Plot HeadModel
    h.btn_headmodel = uicontrol(h.panel_anatomy,'BackgroundColor',[1 .9 1],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_plot_3D_mri.Position([1 3]))+.01 h.btn_plot_3D_mri.Position(2) .2 .03],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','plot HeadModel','Callback',@sm_plot_headmodel);
    %% Plot Lead Field Dipoles
    h.btn_leadfield_grids = uicontrol(h.panel_anatomy,'BackgroundColor',[.9 1 .9],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_headmodel.Position([1 3]))+.01 h.btn_headmodel.Position(2) .2 .03],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','plot LF dipoles','Callback',@sm_plot_leadfield_grid);
    %% Toggle channel names on/off
    h.toggle_chan_text_OnOff = uicontrol(h.panel_anatomy,'BackgroundColor',[.9 .9 1],'ForegroundColor',[0 .6 0],'Style','togglebutton','Units','normalize',...
        'Position',[sum(h.btn_leadfield_grids.Position([1 3]))+.02 h.btn_leadfield_grids.Position(2) .15 .03],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','Labels On','Value',1,'Callback',@chan_text_OnOff_Callback);
    %% Toggle channels on/off
    h.toggle_sens_OnOff = uicontrol(h.panel_anatomy,'BackgroundColor',[.9 .9 1],'ForegroundColor',[0 .6 0],'Style','togglebutton','Units','normalize',...
        'Position',[sum(h.toggle_chan_text_OnOff.Position([1 3]))+.01 h.toggle_chan_text_OnOff.Position(2) h.toggle_chan_text_OnOff.Position(3:4)],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','Sensors On','Value',1,'Callback',@chan_text_OnOff_Callback);
    %% Sliders: Transparencies SCALP
    h.slider_transparency_scalp_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 h.toggle_chan_text_OnOff.Position(2)-.03 .225 .025],...
        'FontSize',10,'HorizontalAlignment','center','String','Scalp Transparency');
    h.slider_transparency_scalp = uicontrol(h.panel_anatomy,'BackgroundColor',[0.8980 0.7549 0.5431],'ForegroundColor',[1 1 1]*0,'Style','slider','Units','normalize',...
        'Position',[h.slider_transparency_scalp_txt.Position(1) h.slider_transparency_scalp_txt.Position(2)-.0275 h.slider_transparency_scalp_txt.Position(3) .025],...
        'FontSize',8,'Value',0.25,'Max',1,'Min',0,'Callback',{@set_transparency,'scalp'});
    
    % Sliders: Transparencies for BRAIN
    h.slider_transparency_brain_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 h.slider_transparency_scalp.Position(2)-.03 h.slider_transparency_scalp.Position(3) h.slider_transparency_scalp_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','Brain Transparency');
    h.slider_transparency_brain = uicontrol(h.panel_anatomy,'BackgroundColor',[1 1 1]*.9,'ForegroundColor',[1 1 1]*0,'Style','slider','Units','normalize',...
        'Position',[h.slider_transparency_brain_txt.Position(1) h.slider_transparency_brain_txt.Position(2)-.0275 h.slider_transparency_brain_txt.Position(3) .025],...
        'FontSize',8,'Value',0.25,'Max',1,'Min',0,'Callback',{@set_transparency,'brain'});
    
    % Sliders: Transparencies for TOPO
    h.slider_transparency_topo_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 h.slider_transparency_brain.Position(2)-.03 h.slider_transparency_brain.Position(3) h.slider_transparency_brain_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','Topo Transparency');
    h.slider_transparency_topo = uicontrol(h.panel_anatomy,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','slider','Units','normalize',...
        'Position',[h.slider_transparency_topo_txt.Position(1) h.slider_transparency_topo_txt.Position(2)-.0275 h.slider_transparency_topo_txt.Position(3) .025],...
        'FontSize',8,'Value',0.25,'Max',1,'Min',0,'Callback',{@set_transparency,'topo'});
    
    % Sliders: Transparencies for HEADMODEL
    h.slider_transparency_hdm_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 h.slider_transparency_topo.Position(2)-.03 h.slider_transparency_topo.Position(3) h.slider_transparency_topo_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','HDM Transparency');
    h.slider_transparency_hdm = uicontrol(h.panel_anatomy,'BackgroundColor',[1 .8 1],'ForegroundColor',[1 1 1]*0,'Style','slider','Units','normalize',...
        'Position',[h.slider_transparency_hdm_txt.Position(1) h.slider_transparency_hdm_txt.Position(2)-.0275 h.slider_transparency_hdm_txt.Position(3) .025],...
        'FontSize',8,'Value',0.25,'Max',1,'Min',0,'Callback',{@set_transparency,'hdm'});
    
    % Sliders: Transparencies for LEADFIELD dipoles
    h.slider_transparency_lf_grids_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 h.slider_transparency_hdm.Position(2)-.03 h.slider_transparency_hdm.Position(3) h.slider_transparency_hdm_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','LF Transparency');
    h.slider_transparency_lf_grids = uicontrol(h.panel_anatomy,'BackgroundColor',[.8 1 .8],'ForegroundColor',[1 1 1]*0,'Style','slider','Units','normalize',...
        'Position',[h.slider_transparency_lf_grids_txt.Position(1) h.slider_transparency_lf_grids_txt.Position(2)-.0275 h.slider_transparency_lf_grids_txt.Position(3) .025],...
        'FontSize',8,'Value',0.25,'Max',1,'Min',0,'Callback',{@set_transparency,'lf'});
    %% Button:  Light on/off
    h.toggle_light_OnOff_anat = uicontrol(h.panel_anatomy,'BackgroundColor',[1 1 1]*.9,'ForegroundColor',[1 0 0],'Style','togglebutton','Units','normalize',...
        'Position',[h.slider_transparency_lf_grids.Position(1) h.slider_transparency_lf_grids.Position(2)-.05 .15 .035],...
        'FontSize',8,'Value',0,'String', 'Light Off','Callback',@toggle_light_OnOff_axes_anatomy);
    %% Axes Anatomy
    h.axes_anatomy = axes(h.panel_anatomy,'Position',[.3 .05 .5 .35]); axis off;
    % slider topo caxis scaling relative to yscale
    %% Sliders for Topo Scale
    h.slider_topo_scale = uicontrol(h.panel_anatomy,'BackgroundColor',[1 1 1]*.9,'ForegroundColor',[1 1 1]*0,'Style','slider','Units','normalize',...
        'Position',[.925 .2 .03 .2],...
        'FontSize',8,'Value',1,'Max',1,'Min',0.01,'Callback',@set_topo_caxis);
    h.slider_topo_scale_txt = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.slider_topo_scale.Position(1)-.02 h.slider_topo_scale.Position(2)-.06  .08 .05],...
        'FontSize',10,'HorizontalAlignment','center','String','Topo Scale');
    h.slider_topo_scale_text_max = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.slider_topo_scale.Position([1 3]))-.05 sum(h.slider_topo_scale.Position([2 4]))+.01  .08 .025],...
        'FontSize',10,'HorizontalAlignment','center','String',['100 ' char(181) 'V']);
    h.slider_topo_scale_text_val = uicontrol(h.panel_anatomy,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.slider_topo_scale.Position(1)-.08 sum(h.slider_topo_scale.Position([2]))  .08 .025],...
        'FontSize',10,'HorizontalAlignment','right','String','');
    
    %% %%%% Panel "Project Source Data --> Sensors" on Tab "Simulate M/EEG" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h.panel_sim_data = uipanel(h.tab_sim_meeg,'Title','Forward Model: Project Source Data --> Sensors','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.5 .01 .49 .98],'Visible','on');
    %% Subpanel Generate Sensor Noise
    h.panel_sens_noise = uipanel(h.panel_sim_data,'Title','Generate Sensor Noise','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.01 .825 .98 .175],'Visible','on');
    %% Menu: Noise Projection Type 'Synthetic', 'Brain' or 'Real'
    h.menu_noise_projection_txt = uicontrol(h.panel_sens_noise,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 .79 .18 .2],...
        'FontSize',10,'HorizontalAlignment','left','String','Noise Projection');
    h.menu_noise_projection = uicontrol(h.panel_sens_noise,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(h.menu_noise_projection_txt.Position([1 3]))+.01 h.menu_noise_projection_txt.Position(2) .25 .2],...
        'FontSize',8,'HorizontalAlignment','center','String',{'Synthetic Sensor Noise' 'Synthetic Brain Noise' 'Synthetic MVAR/GAN Noise (not implemented yet)' 'Real Sensor Noise'},...
        'Value',1,'Callback',@menu_noise_projection_CallBack);
    %     'FontSize',8,'HorizontalAlignment','center','String',{'Sim Sensor' 'Sim Brain' 'Real Sensor' 'GAN Sensor'},'Value',1,'Callback',@menu_noise_projection_CallBack);
    %% btn Sim Brain Noise
    h.btn_sim_sens_noise = uicontrol(h.panel_sens_noise,'BackgroundColor',[1 .9 .7],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.menu_noise_projection.Position([1 3]))+.02 h.menu_noise_projection_txt.Position(2) .175 .23],...
        'FontSize',8,'HorizontalAlignment','center','String','Simulate Noise','Value',1,'Callback',@sim_sens_noise);
    %% btn Get Real Sensor Data for Noise
    h.btn_get_real_noise_dataset = uicontrol(h.panel_sens_noise,'BackgroundColor',[.9 1 .9],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[h.btn_sim_sens_noise.Position(1) h.btn_sim_sens_noise.Position(2) .175 .23],...
        'Visible','off','FontSize',8,'HorizontalAlignment','center','String','Load Noise','Value',1,'Callback',{@sm_load_real_sensor_file,'noise'});
    %% btn FFT Phase Shuffle Sensor Data
    h.btn_fft_radnomize_phase = uicontrol(h.panel_sens_noise,'BackgroundColor',[1 .9 .9],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_get_real_noise_dataset.Position([1 3]))+.02 h.btn_get_real_noise_dataset.Position(2) .175 .23],...
        'Visible','on','FontSize',8,'HorizontalAlignment','center','String','Randomize Phase','Value',1,'Callback',@sm_fft_randomize_phase);
    
    %% Sub-Subpanel Brain Noise Parameters
    % delete(h.panel_brain_noise)
    h.panel_brain_noise = uipanel(h.panel_sens_noise,'Title','Brain Noise Sources','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[0 0 0],...
        'Units','normalized','Position',[.005 .01 .99 .77],'Visible','off');
    %% Edit Number of Noise sources
    h.edit_num_noise_sources_txt = uicontrol(h.panel_brain_noise,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 .65 .18 .35],...
        'FontSize',10,'HorizontalAlignment','left','String','Number per Trial');
    h.edit_num_noise_sources = uicontrol(h.panel_brain_noise,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_num_noise_sources_txt.Position([1 3]))+.01 h.edit_num_noise_sources_txt.Position(2) .07 .35],...
        'FontSize',8,'HorizontalAlignment','center','String','30','Value',1);
    %% Edit Noise Source Amplitudes Mean +/- Stdev
    h.edit_brain_noise_amp_txt = uicontrol(h.panel_brain_noise,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_num_noise_sources.Position([1 3]))+.01 h.edit_num_noise_sources_txt.Position(2) .275 .35],...
        'FontSize',10,'HorizontalAlignment','left','String','Noise Amp [mean stdev]');
    h.edit_brain_noise_amp = uicontrol(h.panel_brain_noise,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_brain_noise_amp_txt.Position([1 3]))+.01 h.edit_num_noise_sources_txt.Position(2) .07 .35],...
        'FontSize',8,'HorizontalAlignment','center','String','30 20','Value',1);
    %% Noise Source Locations
    h.btn_brain_locs_txt = uicontrol(h.panel_brain_noise,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 h.menu_noise_projection_txt.Position(2)-.55 .25 .35],...
        'FontSize',10,'HorizontalAlignment','left','String','Brain Noise Locations:');
    h.btn_rand_brain_noise_locs = uicontrol(h.panel_brain_noise,'BackgroundColor',[1 .9 .9],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_brain_locs_txt.Position([1 3]))+.01 h.btn_brain_locs_txt.Position(2) .15 .35],...
        'FontSize',8,'HorizontalAlignment','center','String','Random Locs','Value',1,'Callback',@btn_rand_brain_noise_locs_Callback);
    h.btn_load_brain_noise_locs = uicontrol(h.panel_brain_noise,'BackgroundColor',[1 1 1]*.9,'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_rand_brain_noise_locs.Position([1 3]))+.01 h.btn_brain_locs_txt.Position(2) .125 .35],...
        'FontSize',8,'HorizontalAlignment','center','String','Load Locs','Value',1,'Callback',@btn_load_brain_noise_locs_Callback);
    h.btn_load_default_mode_noise_locs = uicontrol(h.panel_brain_noise,'BackgroundColor',[.9 1 .9],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_load_brain_noise_locs.Position([1 3]))+.01 h.btn_brain_locs_txt.Position(2) .175 .35],...
        'FontSize',8,'HorizontalAlignment','center','String','Default Mode Locs','Value',1,'Callback',@btn_load_default_mode_noise_locs_Callback);
    h.btn_plot_noise_locs = uicontrol(h.panel_brain_noise,'BackgroundColor',[.9 .9 1],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_load_default_mode_noise_locs.Position([1 3]))+.01 h.btn_brain_locs_txt.Position(2) .175 .35],...
        'FontSize',8,'HorizontalAlignment','center','String','Plot Noise Locs','Value',1,'Callback',@btn_plot_noise_locs);
    
    %% Subpanel "Functions" buttons
    % delete(h.panel_project_btns)
    h.panel_project_btns = uipanel(h.panel_sim_data,'Title','Functions','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.01 .71 .98 .115],'Visible','on');
    %% Edit Leadfield gain factor
    h.edit_leadfield_gain_txt = uicontrol(h.panel_project_btns,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 .6 .175 .325],...
        'FontSize',10,'HorizontalAlignment','left','String','Lead Field Gain');
    h.edit_leadfield_gain = uicontrol(h.panel_project_btns,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_leadfield_gain_txt.Position([1 3]))+.01 h.edit_leadfield_gain_txt.Position(2) .07 .35],...
        'FontSize',8,'HorizontalAlignment','center','String','1e6','Value',1);
    %% Btn Simulate M/EEG by forward projecting the data
    h.btn_run_sim_meeg = uicontrol(h.panel_project_btns,'BackgroundColor',[.9 1 .9],'ForegroundColor',[0 0 0],'Style','pushbutton','Units','normalize',...
        'Position',[.05 .2 .255 .35],'Visible','on','UserData',1,...
        'FontSize',10,'HorizontalAlignment','center','String','Simulate M/EEG','Callback',@sim_meeg);
    %% Radio to open "Preprocessing" Panel
    h.radio_panel_preprocess = uicontrol(h.panel_project_btns,'Style','radiobutton', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[0 .5 .5],'Units','normalize',...
        'Position',[sum(h.btn_run_sim_meeg.Position([1 3]))+.02  h.edit_leadfield_gain_txt.Position(2) .25 .35],...
        'FontSize',10,'HorizontalAlignment','left','String','Preprocessing Panel','Value',0,'Callback',@radio_panel_preprocess_CallBack);
    %% Btn Filter data
    h.btn_filter_data = uicontrol(h.panel_project_btns,'BackgroundColor',[1 1 1]*.9,'ForegroundColor',[0 .5 .5],'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_run_sim_meeg.Position([1 3]))+.02  h.btn_run_sim_meeg.Position(2) .11 .35],'Visible','on','UserData',1,...
        'FontSize',8,'HorizontalAlignment','center','String','Filter Data','Callback',@sm_filter_data_Callback);
    %% Btn Restore Org data
    h.btn_restore_sens_org = uicontrol(h.panel_project_btns,'BackgroundColor',[.9 1 .9],'ForegroundColor',[0 .5 .5],'Style','pushbutton','Units','normalize',...
        'Position',[sum(h.btn_filter_data.Position([1 3]))+.01  h.btn_filter_data.Position(2) .135 .35],'Visible','on','UserData',1,...
        'FontSize',8,'HorizontalAlignment','center','String','Restore Data','Callback',@sm_restore_sens_data_Callback, ...
        'Tooltip','restores sensor data back to orginal data before preprocessing was applied');
    
    %% Btn Load Real M/EEG
    h.btn_load_real_meeg = uicontrol(h.panel_project_btns,'BackgroundColor',[.9 1 .9],'ForegroundColor',[0 0 0],'Style','pushbutton','Units','normalize',...
        'Position',[.805 .65 .185 .35],'Visible','on','UserData',1,...
        'FontSize',8,'HorizontalAlignment','center','String','Load Real M/EEG','Callback',{@sm_load_real_sensor_file,'final'});
    %% Btn Load Source waves from user-defined file
    h.btn_load_sig_waves = uicontrol(h.panel_project_btns,'BackgroundColor',[.8 .8 1],'ForegroundColor',[0 0 0],'Style','pushbutton','Units','normalize',...
        'Position',[.605 .65 .185 .35],'Visible','on','UserData',1,...
        'FontSize',8,'HorizontalAlignment','center','String','Load Source Waves','Callback',{@sm_load_source_waves,'sg'});
    %% Redundant buttons because clicking on the listbox items will plot the waveforms
    %     %% Btn Plot Source Waves
    %     h.btn_plot_source_data = uicontrol(h.panel_project_btns,'BackgroundColor',[.9 .8 1],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
    %         'Position',[.01 .05 .185 .4],...
    %         'FontSize',8,'HorizontalAlignment','center','String','Plot Source Waves','Value',1,'Callback',@plot_source_data);
    %% btn Save BRANE Lab dataset <-- Not Implemented
    % h.btn_save_bl_data = uicontrol(h.panel_project_btns,'BackgroundColor',[1 .8 .8],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
    %     'Position',[sum(h.btn_plot_sens_data.Position([1 3]))+.01 h.btn_run_sim_meeg.Position(2) .185 .4],...
    %     'FontSize',8,'HorizontalAlignment','center','String','Save BRANELab','Value',1,'Callback',@save_bl_data);
    % %% btn Save SimMEEG dataset
    % h.btn_save_bs_data = uicontrol(h.panel_project_btns,'BackgroundColor',[1 .8 .8],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
    %     'Position',[sum(h.btn_plot_sens_data.Position([1 3]))+.01 h.btn_run_sim_meeg.Position(2) .185 .4],...
    %     'FontSize',8,'HorizontalAlignment','center','String','Save SimMEEG','Value',1,'Callback',@save_bs_data);
    %% btn Save BESA dataset
    h.btn_save_BESA_data = uicontrol(h.panel_project_btns,'BackgroundColor',[1 .8 .8],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.605 .2 .185 .35],...
        'FontSize',8,'HorizontalAlignment','center','String','Write BESA','Value',1,'Callback',@save_BESA_data);
    %% btn Save Field Trip dataset
    h.btn_save_ft_data = uicontrol(h.panel_project_btns,'BackgroundColor',[1 .9 .8],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.805 .2 .185 .35],...
        'FontSize',8,'HorizontalAlignment','center','String','Save Field Trip','Value',1,'Callback',@save_ft_data);
    %% Subpanel "Waveform Data"
    % delete(h.panel_sens_data)
    h.panel_sens_data = uipanel(h.panel_sim_data,'Title','Waveform Data','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.01 .01 .98 .7],'Visible','on');
    %% ListBox Sources to Plot Data
    h.listbox_sources_txt = uicontrol(h.panel_sens_data,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[.6 .4 .9],'Units','normalize',...
        'Position',[.01 .95 .11 .035],...
        'FontSize',10,'HorizontalAlignment','left','String','Sources','FontWeight','bold');
    h.listbox_sources = uicontrol(h.panel_sens_data,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[.6 .4 .9],'Style','listbox','Units','normalize',...
        'Position',[.01 h.listbox_sources_txt.Position(2)-.11 .12 .11],...
        'FontSize',8,'HorizontalAlignment','center','String',{'Source 1' 'Source 2' 'Source 3'},'Value',1,'Callback',@plot_source_data);
    % changing color of list box items for sources 1, 2, and 3
    Data(1).name = 'Source 1'; Data(1).Color = h.src_clr(1,:)*255;
    Data(2).name = 'Source 2'; Data(2).Color = h.src_clr(2,:)*255;
    Data(3).name = 'Source 3'; Data(3).Color = h.src_clr(3,:)*255;
    pre = '<HTML><FONT color="'; post = '</FONT></HTML>';
    for v=1:3
        clr_str = reshape( dec2hex( Data(v).Color, 2 )',1, 6);
        h.listbox_sources.UserData.clr_str{v} = [pre clr_str '">' Data(v).name post];
    end
    h.listbox_sources.String = h.listbox_sources.UserData.clr_str;
    %% ListBox Sensors to Plot Data
    h.listbox_chans_txt = uicontrol(h.panel_sens_data,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.chan_clr*.45,'Units','normalize',...
        'Position',[.01 h.listbox_sources.Position(2)-.05 .11 .035],...
        'FontSize',10,'HorizontalAlignment','left','String','Sensors','FontWeight','bold');
    h.listbox_chans = uicontrol(h.panel_sens_data,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.chan_clr*.45,'Style','listbox','Units','normalize',...
        'Position',[.01 .1 .12 h.listbox_chans_txt.Position(2)-.1],...
        'FontSize',8,'HorizontalAlignment','center','String','','Value',1,'Max',2,'Callback',@plot_sens_data);
    %% Btn Plot Sensor Waves
    h.btn_plot_sens_data = uicontrol(h.panel_sens_data,'BackgroundColor',[1 1 .8],'ForegroundColor',[.8 1 1]*.2,'Style','pushbutton','Units','normalize',...
        'Position',[h.listbox_chans.Position(1) h.listbox_chans.Position(2)-.06 .15 .05],...
        'FontSize',8,'HorizontalAlignment','center','String','RePlot Sensors','Value',1,'Callback',@plot_sens_data);
    %% Axes plotting Noise, Signal, Final
    h.axes_sens_noise = axes(h.panel_sens_data,'Position',[.23 .675 .95-.23 .225]);
    h.axes_sens_signal = axes(h.panel_sens_data,'Position',[.23 .375 .95-.23 .225]);
    h.axes_sens_final = axes(h.panel_sens_data,'Position',[.23  .075 .95-.23 .225]); xlabel('Time (sec)'); ylabel('Amp (microV)');
    h.edit_yscale_txt = uicontrol(h.panel_sens_data,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.23 .95 .15 .035],...
        'FontSize',10,'HorizontalAlignment','left','String',['Y Scale (' char(181) 'V):'],'Value',1);
    h.edit_yscale = uicontrol(h.panel_sens_data,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_yscale_txt.Position([1 3]))+.01 h.edit_yscale_txt.Position(2) .125 .045],...
        'FontSize',10,'HorizontalAlignment','center','String','-100 100','Value',1,'Callback',@update_y_scale);
    %% Btn Topography
    h.btn_plot_topo = uicontrol(h.panel_sens_data,'BackgroundColor',[1 1 .8],'ForegroundColor',[.8 1 1]*.2,'Style','pushbutton','Units','normalize',...
        'Position',[.8 h.edit_yscale_txt.Position(2) .15 .05],...
        'FontSize',8,'HorizontalAlignment','center','String','Plot Topo','Value',1,'Callback',@plot_topo_data);
    
    %% %%%% Panel "Preprocessing & Analysis" on Tab "Simulate M/EEG" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delete(h.panel_preprocess);
    h.panel_preprocess = uipanel(h.tab_sim_meeg,'Title','Preprocessing & Analysis','FontSize',10,'BackgroundColor',[.96 1 1],'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.01 .01 .49 .98],'Visible','off');
    
    %% SubPanel Filtering
    % delete(h.panel_filtering)
    h.panel_filtering = uipanel(h.panel_preprocess,'Title','Filtering Parameters','FontSize',10,'BackgroundColor',h.panel_preprocess.BackgroundColor,'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.01 1-.12 .98 .12],'Visible','on');
    %% Menu: Filter Type
    h.menu_filt_type_txt = uicontrol(h.panel_filtering,'Style','text', 'BackgroundColor',h.panel_filtering.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 .65 .135 .28],...
        'FontSize',10,'HorizontalAlignment','left','String','Filter Type');
    h.menu_filt_type = uicontrol(h.panel_filtering,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(h.menu_filt_type_txt.Position([1 3]))+.01 h.menu_filt_type_txt.Position(2) .175 .28],...
        'FontSize',8,'HorizontalAlignment','center','String',{'bandpassfir' 'highpassfir' 'lowpassfir' 'bandstopfir' 'bandpassiir' 'highpassiir' 'lowpassiir' 'bandstopiir' },...
        'Value',1,'Callback',@menu_filt_type_Callback);
    %% Menu: Filter Method
    h.menu_filt_method_txt = uicontrol(h.panel_filtering,'Style','text', 'BackgroundColor',h.panel_filtering.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 h.menu_filt_type.Position(2)-.42 .155 .28],...
        'FontSize',10,'HorizontalAlignment','left','String','Filter Method');
    h.menu_filt_method = uicontrol(h.panel_filtering,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(h.menu_filt_method_txt.Position([1 3]))+.01 h.menu_filt_method_txt.Position(2) .155 .28],...
        'FontSize',8,'HorizontalAlignment','center','String',{'equiripple' 'kaiserwin'},'Value',1,'Callback',@menu_filt_type_Callback);
    %% Edit: Filter Order
    h.edit_filt_order_txt = uicontrol(h.panel_filtering,'Style','text', 'BackgroundColor',h.panel_filtering.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_filt_method.Position([1 3]))+.01 h.menu_filt_method.Position(2) .12 h.menu_filt_method_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','left','String','Filter Order');
    h.edit_filt_order = uicontrol(h.panel_filtering,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_filt_order_txt.Position([1 3]))+.01 h.edit_filt_order_txt.Position(2) .045 h.edit_filt_order_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Callback',@menu_filt_type_Callback,...
        'Tooltip','if set to 0 then default is to use the PassbandRipple and StopbandAttenuation to determine minimum filter order');
    %% Edit: Filter Bands - [lpStartFreq lpStopFreq hpStartFreq hpStopFreq]
    h.edit_filt_band_txt = uicontrol(h.panel_filtering,'Style','text', 'BackgroundColor',h.panel_filtering.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_filt_type.Position([1 3]))+.01 h.menu_filt_type.Position(2) .125 h.menu_filt_type_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','left','String','Filter Band');
    h.edit_filt_band = uicontrol(h.panel_filtering,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_filt_band_txt.Position([1 3]))+.005 h.edit_filt_band_txt.Position(2) .155 h.edit_filt_band_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','0 0 0 0',...
        'Tooltip',sprintf('freqs = frequencies for  [lpStartFreq lpStopFreq hpStartFreq hpStopFreq]\nlowpass filter at 20-Hz with 4-Hz stop-band transiton = [0 0 20 24].\nhighpass filter at 10-Hz with 2-Hz stop-band transiton = [8 10 0 0].\nbandpass filter at 10- to 20-Hz with 2- and 4-Hz stop-band transitons = [8 10 20 24].\nbandstop filter at 10- to 20-Hz with 2- and 4-Hz stop-band transitons = [10 12 18 20].\n') );
    %% Edit: StopbandAttenuation
    h.edit_filt_StopbandAttenuation_txt = uicontrol(h.panel_filtering,'Style','text', 'BackgroundColor',h.panel_filtering.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_filt_band.Position([1 3]))+.01 h.edit_filt_band.Position(2) .29 h.edit_filt_band_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','left','String','StopBand Attenuation (dB)');
    h.edit_filt_StopbandAttenuation = uicontrol(h.panel_filtering,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_filt_StopbandAttenuation_txt.Position([1 3]))+.005 h.edit_filt_StopbandAttenuation_txt.Position(2) .05 h.edit_filt_StopbandAttenuation_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','60');
    %% Radio: plot filter design
    h.radio_plot_filter_design = uicontrol(h.panel_filtering,'Style','radiobutton', 'BackgroundColor',h.panel_preprocess.BackgroundColor,'Foregroundcolor',[0 .5 .5],'Units','normalize',...
        'Position',[sum(h.edit_filt_order.Position([1 3]))+.05  .1 .2 .35],...
        'FontSize',8,'HorizontalAlignment','left','String','plot filter design','Value',1);
    %% Btn Create Filter Design
    h.btn_create_filter_design = uicontrol(h.panel_filtering,'BackgroundColor',[1 1 1]*.9,'ForegroundColor',[0 .5 .5],'Style','pushbutton','Units','normalize',...
        'Position',[.835 .1 .15 .35],'Visible','on','UserData',1,'FontWeight','normal',...
        'FontSize',10,'HorizontalAlignment','center','String','Create Filter','Callback',@sm_filter_data_Callback,...
        'Tooltip',sprintf('see help in "filter_data_new.m" for more details'));
    
    %% Sub-Panel: "FieldTrip's Time-Frequency Response (TFR) Analysis" before source modeling
    % delete(h.panel_inv_ft_tfr_params)
    h.panel_inv_ft_tfr_params = uipanel(h.panel_preprocess ,'Title','FieldTrip''s Time-Frequency Response (TFR) Analysis','FontSize',10,...
        'BackgroundColor',[1 1 .9],'Foregroundcolor',[.8 .25 .25],...
        'Units','normalize','Position',[.01 h.panel_filtering.Position(2)-.21 .98 .2],'Visible','on');
    %% Menu: FieldTrip TFR "Method"
    h.menu_preproc_ft_tfr_method_txt = uicontrol(h.panel_inv_ft_tfr_params,'Style','text', 'BackgroundColor',h.panel_inv_ft_tfr_params.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 .725 .125 .2],...
        'FontSize',10,'HorizontalAlignment','left','String','Method');
    h.menu_preproc_ft_tfr_method = uicontrol(h.panel_inv_ft_tfr_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[sum(h.menu_preproc_ft_tfr_method_txt.Position([1 3]))+.01 h.menu_preproc_ft_tfr_method_txt.Position(2)+.025 .12 h.menu_preproc_ft_tfr_method_txt.Position(4)],...
        'FontSize',8,'HorizontalAlignment','center','String',{'mtmfft' 'mtmconvol' 'wavelet' 'tfr' 'mvar'},'Value',1,'Callback',@sm_ft_tfr_menu_Callback,...
        'Tooltip',sprintf('see "ft_freqanalysis.m" for more information'));
    %% Menu: FieldTrip TFR "Taper Type"
    h.menu_preproc_ft_tfr_taper_txt = uicontrol(h.panel_inv_ft_tfr_params,'Style','text', 'BackgroundColor',h.panel_inv_ft_tfr_params.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 h.menu_preproc_ft_tfr_method_txt.Position(2)-.23 .1 .225],...
        'FontSize',10,'HorizontalAlignment','left','String','Taper Type');
    h.menu_preproc_ft_tfr_taper = uicontrol(h.panel_inv_ft_tfr_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[h.menu_preproc_ft_tfr_method.Position(1) h.menu_preproc_ft_tfr_taper_txt.Position(2)+.025 .12 h.menu_preproc_ft_tfr_taper_txt.Position(4)],...
        'FontSize',8,'HorizontalAlignment','center','String',{'dpss' 'hanning'},'Value',1,...
        'Tooltip',sprintf('"dpss", "hanning" for now'));
    %% Menu: FieldTrip TFR "Output Type"
    h.menu_preproc_ft_tfr_output_txt = uicontrol(h.panel_inv_ft_tfr_params,'Style','text', 'BackgroundColor',h.panel_inv_ft_tfr_params.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[.01 h.menu_preproc_ft_tfr_taper_txt.Position(2)-.23 .1 .225],...
        'FontSize',10,'HorizontalAlignment','left','String','Output Type');
    h.menu_preproc_ft_tfr_output = uicontrol(h.panel_inv_ft_tfr_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[h.menu_preproc_ft_tfr_method.Position(1) h.menu_preproc_ft_tfr_output_txt.Position(2)+.025 .12 h.menu_preproc_ft_tfr_output_txt.Position(4)],...
        'FontSize',8,'HorizontalAlignment','center','String',{'pow' 'powandcsd' 'fourier'},'Value',1,'Callback',@sm_ft_tfr_menu_Callback,...
        'Tooltip',sprintf('"pow" return the power-spectra\n"powandcsd" return the power and the cross-spectra\n"fourier" return the complex Fourier-spectra'));
    %% Edit: FieldTrip TFR "Frequency Bands"
    h.edit_preproc_ft_tfr_foilim_txt = uicontrol(h.panel_inv_ft_tfr_params,'Style','text', 'BackgroundColor',h.panel_inv_ft_tfr_params.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_preproc_ft_tfr_method.Position([1 3]))+.01 h.menu_preproc_ft_tfr_method_txt.Position(2) .185 h.menu_preproc_ft_tfr_method_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','left','String','Frequency Bands (Hz)');
    h.edit_preproc_ft_tfr_foilim = uicontrol(h.panel_inv_ft_tfr_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_preproc_ft_tfr_foilim_txt.Position([1 3]))+.005 h.menu_preproc_ft_tfr_method.Position(2) .075 h.menu_preproc_ft_tfr_method_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','1 20','Callback',@sm_ft_tfr_menu_Callback,...
        'Tooltip',sprintf('if 2 numbers then frequency band of interest [begin end] using "cfg.foilim"\nif 1 number or >2 numbers then individual frequencies using "cfg.foi"'));
    %% Edit: FieldTrip TFR "Taper Smoothing"
    h.edit_preproc_ft_tfr_tapsmofrq_txt = uicontrol(h.panel_inv_ft_tfr_params,'Style','text', 'BackgroundColor',h.panel_inv_ft_tfr_params.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_preproc_ft_tfr_taper.Position([1 3]))+.01 h.menu_preproc_ft_tfr_taper_txt.Position(2) .185 h.menu_preproc_ft_tfr_method_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','left','String','Taper Smoothing (Hz)','Visible','off');
    h.edit_preproc_ft_tfr_tapsmofrq = uicontrol(h.panel_inv_ft_tfr_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_preproc_ft_tfr_tapsmofrq_txt.Position([1 3]))+.005 h.menu_preproc_ft_tfr_taper.Position(2) .075 h.edit_preproc_ft_tfr_tapsmofrq_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','4','Callback',@sm_ft_tfr_menu_Callback,'Visible','off',...
        'Tooltip',sprintf('amount of spectral smoothing through\nmulti-tapering of individual frequencies in "Frequency Band".\nNote, 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.'));
    %% Edit: FieldTrip TFR "t_ftimwin"
    h.edit_preproc_ft_tfr_t_ftimwin_txt = uicontrol(h.panel_inv_ft_tfr_params,'Style','text', 'BackgroundColor',h.panel_inv_ft_tfr_params.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_preproc_ft_tfr_taper.Position([1 3]))+.01 h.menu_preproc_ft_tfr_output_txt.Position(2) .185 h.menu_preproc_ft_tfr_method_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','left','String','Time Window (sec)','Visible','off');
    h.edit_preproc_ft_tfr_t_ftimwin = uicontrol(h.panel_inv_ft_tfr_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_preproc_ft_tfr_t_ftimwin_txt.Position([1 3]))+.005 h.menu_preproc_ft_tfr_output.Position(2) .075 h.edit_preproc_ft_tfr_t_ftimwin_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','0.1','Callback',@sm_ft_tfr_menu_Callback,'Visible','off',...
        'Tooltip',sprintf('vector 1 x numfoi, length of time window (in seconds)'));
    %% Edit: FieldTrip TFR "toi" - center(s) of time-window above
    h.edit_preproc_ft_tfr_toi_txt = uicontrol(h.panel_inv_ft_tfr_params,'Style','text', 'BackgroundColor',h.panel_inv_ft_tfr_params.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_preproc_ft_tfr_taper.Position([1 3]))+.01 h.edit_preproc_ft_tfr_t_ftimwin_txt.Position(2)-.23 .185 h.menu_preproc_ft_tfr_method_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','left','String','Time Overlap (%)','Visible','off');
    h.edit_preproc_ft_tfr_toi = uicontrol(h.panel_inv_ft_tfr_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_preproc_ft_tfr_toi_txt.Position([1 3]))+.005 h.edit_preproc_ft_tfr_toi_txt.Position(2)+.025 .075 h.edit_preproc_ft_tfr_toi_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','50','Callback',@sm_ft_tfr_menu_Callback,'Visible','off',...
        'Tooltip',sprintf('vector 1 x numtoi, the times on which the analysis\nwindows should be centered (in seconds), or a string\nsuch as "50%" or "all" (default).  Both string options\nuse all timepoints available in the data, but "all"\ncenters a spectral estimate on each sample, whereas\nthe percentage specifies the degree of overlap between\nthe shortest time windows from cfg.t_ftimwin'));
    %% Edit: FieldTrip TFR "wavelet width"
    h.edit_preproc_ft_tfr_wt_width_txt = uicontrol(h.panel_inv_ft_tfr_params,'Style','text', 'BackgroundColor',h.panel_inv_ft_tfr_params.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_preproc_ft_tfr_toi.Position([1 3]))+.01 h.menu_preproc_ft_tfr_method_txt.Position(2) .195 h.menu_preproc_ft_tfr_method_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','left','String','Wavelet Width (cycles)','Visible','off');
    h.edit_preproc_ft_tfr_wt_width = uicontrol(h.panel_inv_ft_tfr_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_preproc_ft_tfr_wt_width_txt.Position([1 3]))+.005 h.edit_preproc_ft_tfr_wt_width_txt.Position(2)+.025 .075 h.edit_preproc_ft_tfr_wt_width_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','7','Callback',@sm_ft_tfr_menu_Callback,'Visible','off',...
        'Tooltip',sprintf('"cfg.width" = number of cycles, of the wavelet (default = 7)'));
    %% Edit: FieldTrip TFR "wavelet StdDev"
    h.edit_preproc_ft_tfr_gwidth_txt = uicontrol(h.panel_inv_ft_tfr_params,'Style','text', 'BackgroundColor',h.panel_inv_ft_tfr_params.BackgroundColor,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_preproc_ft_tfr_toi.Position([1 3]))+.01 h.menu_preproc_ft_tfr_taper_txt.Position(2) .195 h.menu_preproc_ft_tfr_method_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','left','String','Wavelet StdDev','Visible','off');
    h.edit_preproc_ft_tfr_gwidth = uicontrol(h.panel_inv_ft_tfr_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_preproc_ft_tfr_gwidth_txt.Position([1 3]))+.005 h.edit_preproc_ft_tfr_gwidth_txt.Position(2)+.025 .075 h.edit_preproc_ft_tfr_gwidth_txt.Position(4)],...
        'FontSize',10,'HorizontalAlignment','center','String','3','Callback',@sm_ft_tfr_menu_Callback,'Visible','off',...
        'Tooltip',sprintf('"cfg.gwidth" = determines the length of the used wavelets in standard\ndeviations of the implicit Gaussian kernel and should\nbe choosen >= 3; (default = 3)'));
    
    %% Btn: FieldTrip TFR "Calc FT TFR"
    h.btn_preproc_ft_ = uicontrol(h.panel_inv_ft_tfr_params,'BackgroundColor',[1 1 1]*.9,'ForegroundColor',h.panel_inv_ft_tfr_params.ForegroundColor,'Style','pushbutton','Units','normalize',...
        'Position',[.835 .025 .15 .2],'Visible','on','FontWeight','normal',...
        'FontSize',9,'HorizontalAlignment','center','String','Calc FT TFR','Callback',@sm_ft_freqanalysis);
    
    
    
    
    %%  %%%% "Loading and Saving" SimMEEG datasets
    %% Load SimMEEG dataset
    %      delete(h.btn_load_SimMEEG_dataset)
    h.btn_load_SimMEEG_dataset = uicontrol(h.tab_sim_meeg,'BackgroundColor',[.9 1 .9],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.78 .96 .1 .035],'Visible','on',...
        'FontSize',9,'HorizontalAlignment','center','String','Load SimMEEG Dataset','Callback',@sm_load_SimMEEG_dataset);
    %% Save SimMEEG dataset
    %      delete(h.btn_save_SimMEEG_dataset)
    h.btn_save_SimMEEG_dataset = uicontrol(h.tab_sim_meeg,'BackgroundColor',[1 .8 .8],'ForegroundColor',[1 1 1]*0,'Style','pushbutton','Units','normalize',...
        'Position',[.89 .96 .1 .035],'Visible','on',...
        'FontSize',9,'HorizontalAlignment','center','String','Save SimMEEG Dataset','Callback',@sm_save_SimMEEG_dataset);
    
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tab "Monte Carlo" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%% Panel "Monte Carlo Parameters"  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delete(h.tab_monte_carlo); h.tab_monte_carlo = uitab(h.tabgrp,'Title','Monte Carlo','BackgroundColor',[1 1 1],'Foregroundcolor',[0 .2 .6]);
    h.panel_monte_carlo = uipanel(h.tab_monte_carlo,'Title','Monte Carlo Parameters','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[0 0 0],...
        'Units','normalize','Position',[.005 .91 .99 .09],'Visible','on');
    %% Num Reruns --> conduct full Monte-Carlo again
    h.edit_monte_num_reruns_txt = uicontrol(h.panel_monte_carlo,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 .5 .15 .45],'FontSize',10,'HorizontalAlignment','left','String','Number of Monte Carlo Full ReRuns');
    h.edit_monte_num_reruns = uicontrol(h.panel_monte_carlo,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_monte_num_reruns_txt.Position([1 3]))+.001 h.edit_monte_num_reruns_txt.Position(2) .03 .45],...
        'FontSize',10,'HorizontalAlignment','center','String','1',...
        'Tooltip',sprintf('Number of Full ReRuns of the Monte-Carlo Simulations'));
    %% Number Simulations per Simulation Parameter
    h.edit_monte_num_sims_txt = uicontrol(h.panel_monte_carlo,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_monte_num_reruns.Position([1 3]))+.005 .5 .2 .45],'FontSize',10,'HorizontalAlignment','left','String','Number of Simulation Runs/Simulation Parameter');
    h.edit_monte_num_sims = uicontrol(h.panel_monte_carlo,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_monte_num_sims_txt.Position([1 3]))+.001 h.edit_monte_num_sims_txt.Position(2) .03 .45],...
        'FontSize',10,'HorizontalAlignment','center','String','100',...
        'Tooltip',sprintf('For each "Simulation Run" a new set of source waveforms with new PLV/PLI values will be simulated.\n One Simulation Run is performed for each of the Range Values for each of the "Simulation Parameters" below (e.g., Amplitudes, Locations, Orientations, and SNR).'));
    %% Run # to start Simulations
    h.edit_monte_run_num_txt = uicontrol(h.panel_monte_carlo,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 .05 .15 .45],'FontSize',10,'HorizontalAlignment','left','String','Start Number for Simulation Runs');
    h.edit_monte_run_num = uicontrol(h.panel_monte_carlo,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_monte_run_num_txt.Position([1 3]))+.001 h.edit_monte_run_num_txt.Position(2) .03 .45],...
        'FontSize',10,'HorizontalAlignment','center','String','1',...
        'Tooltip',sprintf('Starting number to begin at for naming the simulation runs'));
    %% Btn Update Monte Carlo
    h.btn_update_monte_carlo = uicontrol(h.panel_monte_carlo,'BackgroundColor',[1 .8 .8],'ForegroundColor','k','Style','pushbutton','Units','normalize',...
        'Position',[.875 .56 .12 .48],...
        'FontSize',10,'HorizontalAlignment','center','String','Update Monte Carlo','Value',1,'Callback',@update_monte_carlo);
    %% Btn Run Monte Carlo
    h.btn_run_monte_carlo = uicontrol(h.panel_monte_carlo,'BackgroundColor',[.8 1 .8],'ForegroundColor','k','Style','pushbutton','Units','normalize',...
        'Position',[.875 .05 .12 .48],...
        'FontSize',10,'HorizontalAlignment','center','String','Run Monte Carlo','Value',1,'Callback',@run_monte_carlo);
    %% Radio to run Source Modeling
    h.radio_monte_source_modeling = uicontrol(h.panel_monte_carlo,'Style','radiobutton', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.btn_update_monte_carlo.Position(1)-.165 h.btn_update_monte_carlo.Position(2) .15 .48],...
        'FontSize',10,'HorizontalAlignment','left','String','Perform Inverse Modeing','Value',1,'Callback',@radio_monte_source_modeling_CallBack);
    %% Radio convert BS to FieldTrip format and save with Monte Carlo run
    h.radio_monte_save_FT_data = uicontrol(h.panel_monte_carlo,'Style','radiobutton', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.btn_run_monte_carlo.Position(1)-.165 h.btn_run_monte_carlo.Position(2) .15 .48],...
        'FontSize',10,'HorizontalAlignment','left','String','Save with Field Trip format','Value',1);
    
    %% %%%% Panel "Study Parameters"  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delete(h.panel_monte_study_params);
    h.panel_monte_study_params = uipanel(h.tab_monte_carlo,'Title','Study & Source PLV/PLI Parameters','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[.6 0 1],...
        'Units','normalize','Position',[.005 h.panel_monte_carlo.Position(2)-.15 .99 .15],'Visible','on');
    %% txt and ui sizes
    txt_xsize = .2; txt_ysize = .2;
    ui_xsize = .08; ui_ysize = .2;
    %% SNR Range
    h.edit_monte_SNR_range_txt = uicontrol(h.panel_monte_study_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 .76 txt_xsize txt_ysize],...
        'FontSize',10,'HorizontalAlignment','left','String','Signal-to-Noise Range (min:step:max)');
    h.edit_monte_SNR_range = uicontrol(h.panel_monte_study_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[sum(h.edit_monte_SNR_range_txt.Position([1 3]))+.005 h.edit_monte_SNR_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','-10 : 3 : -1',...
        'Tooltip',sprintf('Mean +/- StDev (set in last column) for SNR values for each Simulation Run'));
    h.edit_monte_SNR_StdDev = uicontrol(h.panel_monte_study_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[sum(h.edit_monte_SNR_range.Position([1 3]))+.005 h.edit_monte_SNR_range.Position(2) ui_xsize/2 ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','1','Tooltip',sprintf('Standard Deviation for the Means in the Range Values to the left.\n Set to 0 for fixed Means.\n Set to 0 for fixed Means.'));
    %% Number of Trials
    h.edit_monte_num_trials_txt = uicontrol(h.panel_monte_study_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 h.edit_monte_SNR_range_txt.Position(2)-txt_ysize-.02 txt_xsize txt_ysize],...
        'FontSize',10,'HorizontalAlignment','left','String','Number Trials Range (min:step:max)');
    h.edit_monte_num_trials = uicontrol(h.panel_monte_study_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[h.edit_monte_SNR_range.Position(1) h.edit_monte_num_trials_txt.Position(2) h.edit_monte_SNR_range.Position(3:4)],...
        'FontSize',10,'HorizontalAlignment','center','String','90 : 9 : 108',...
        'Tooltip',sprintf('Mean +/- StDev (set in last column) Change in PLV & PLI relative to initial simulation PLV & PLI for each Simulation Run\nSet to 0 for fixed PLV & PLI per Simulation Run\n'));
    h.edit_monte_num_trials_StdDev = uicontrol(h.panel_monte_study_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[h.edit_monte_SNR_StdDev.Position(1) h.edit_monte_num_trials.Position(2) h.edit_monte_SNR_StdDev.Position(3:4)],...
        'FontSize',10,'HorizontalAlignment','center','String','9','Tooltip',sprintf('Standard Deviation for the Means in the Range Values to the left.\n Set to 0 for fixed Means.'));
    %% PLV/PLI Range Change
    h.edit_monte_plv_range_txt = uicontrol(h.panel_monte_study_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 h.edit_monte_num_trials_txt.Position(2)-txt_ysize-.02 txt_xsize txt_ysize],...
        'FontSize',10,'HorizontalAlignment','left','String','PLV/PLI Change Range (min:step:max)');
    h.edit_monte_plv_range = uicontrol(h.panel_monte_study_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[h.edit_monte_num_trials.Position(1) h.edit_monte_plv_range_txt.Position(2) h.edit_monte_num_trials.Position(3:4)],...
        'FontSize',10,'HorizontalAlignment','center','String','-0.1 : 0.1 : 0.1',...
        'Tooltip',sprintf('Mean +/- StDev (set in last column) Change in PLV & PLI relative to initial simulation PLV & PLI for each Simulation Run\nSet to 0 for fixed PLV & PLI per Simulation Run\n'));
    h.edit_monte_plv_StdDev = uicontrol(h.panel_monte_study_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[h.edit_monte_num_trials_StdDev.Position(1) h.edit_monte_plv_range.Position(2) h.edit_monte_num_trials_StdDev.Position(3:4)],...
        'FontSize',10,'HorizontalAlignment','center','String','0.05','Tooltip',sprintf('Standard Deviation for the Means in the Range Values to the left.\n Set to 0 for fixed Means.'));
    
    %% %%%% Panel "Simulate M/EEG Parameters"  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delete(h.panel_monte_random_params);
    h.panel_monte_random_params = uipanel(h.tab_monte_carlo,'Title','Simulate M/EEG Parameters','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[0 .6 .6],...
        'Units','normalize','Position',[.005 h.panel_monte_study_params.Position(2)-.25 .99 .25],'Visible','on');
    %% Source 1 Amp Range
    txt_xsize = .21; txt_ysize = .125;
    ui_xsize = .06; ui_ysize = .125;
    %% Text Titles
    h.edit_monte_amp_range_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 1-(2*txt_ysize)-.01 txt_xsize txt_ysize],'FontSize',10,'HorizontalAlignment','left','String','Amplitude Range (min:step:max)');
    h.edit_monte_loc_range_X_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 h.edit_monte_amp_range_txt.Position(2)-txt_ysize-.01 txt_xsize txt_ysize],'FontSize',10,'HorizontalAlignment','left','String','Location Change Range X (min:step:max)');
    h.edit_monte_loc_range_Y_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 h.edit_monte_loc_range_X_txt.Position(2)-txt_ysize-.01 txt_xsize txt_ysize],'FontSize',10,'HorizontalAlignment','left','String','Location Change Range Y (min:step:max)');
    h.edit_monte_loc_range_Z_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 h.edit_monte_loc_range_Y_txt.Position(2)-txt_ysize-.01 txt_xsize txt_ysize],'FontSize',10,'HorizontalAlignment','left','String','Location Change Range Z (min:step:max)');
    h.edit_monte_ori_Az_range_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 h.edit_monte_loc_range_Z_txt.Position(2)-txt_ysize-.01 txt_xsize txt_ysize],'FontSize',10,'HorizontalAlignment','left','String','Orientation Range Az (min:step:max)');
    h.edit_monte_ori_El_range_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.005 h.edit_monte_ori_Az_range_txt.Position(2)-txt_ysize-.01 txt_xsize txt_ysize],'FontSize',10,'HorizontalAlignment','left','String','Orientation Range El (min:step:max)');
    %% Source titles
    h.edit_monte_source1_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(1,:),'Units','normalize',...
        'Position',[.005+txt_xsize 1-(1*txt_ysize)-.01 ui_xsize txt_ysize],'FontSize',10,...
        'HorizontalAlignment','center','String','Source 1');
    h.edit_monte_source2_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(2,:),'Units','normalize',...
        'Position',[sum(h.edit_monte_source1_txt.Position([1 3]))+.005 h.edit_monte_source1_txt.Position(2) ui_xsize txt_ysize],'FontSize',10,...
        'HorizontalAlignment','center','String','Source 2');
    h.edit_monte_source3_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.src_clr(3,:),'Units','normalize',...
        'Position',[sum(h.edit_monte_source2_txt.Position([1 3]))+.005 h.edit_monte_source1_txt.Position(2) ui_xsize txt_ysize],'FontSize',10,...
        'HorizontalAlignment','center','String','Source 3');
    h.edit_monte_range_StdDev_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.edit_monte_source3_txt.Position([1 3]))+.005 h.edit_monte_source1_txt.Position(2) ui_xsize txt_ysize],'FontSize',10,...
        'HorizontalAlignment','center','String','Std Dev','Tooltip',sprintf('Standard Deviation for the Means in the Range Values to the left.\n Set to 0 for fixed Means.'));
    %% Source Amp Range
    h.edit_monte_source_amp_range(1) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(1,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source1_txt.Position(1) h.edit_monte_amp_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','30',...
        'Tooltip',sprintf('Mean +/- StDev (set in last column) for Source Amplitude Distribution to be sampled from for each Simulation Run\n Setting a range 10:20:50 will generate 100 Simulation Runs for the 3 Amplitude values sampled from a distribution with a mean at these values +/- StdDev.'));
    h.edit_monte_source_amp_range(2) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(2,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source2_txt.Position(1) h.edit_monte_amp_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','30',...
        'Tooltip',sprintf('Mean +/- StDev (set in last column) for Source Amplitude Distribution to be sampled from for each Simulation Run\n Setting a range 10:20:50 will generate 100 Simulation Runs for the 3 Amplitude values sampled from a distribution with a mean at these values +/- StdDev.'));
    h.edit_monte_source_amp_range(3) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(3,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source3_txt.Position(1) h.edit_monte_amp_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','30',...
        'Tooltip',sprintf('Mean +/- StDev (set in last column) for Source Amplitude Distribution to be sampled from for each Simulation Run\n Setting a range 10:20:50 will generate 100 Simulation Runs for the 3 Amplitude values sampled from a distribution with a mean at these values +/- StdDev.'));
    h.edit_monte_source_amp_StdDev = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[h.edit_monte_range_StdDev_txt.Position(1) h.edit_monte_amp_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','5','Tooltip',sprintf('Standard Deviation for the Means in the Range Values to the left.\n Set to 0 for fixed Means.'));
    %% Source Location Change Range - relative to True source locations
    % X Range
    h.edit_monte_source_loc_range_X(1) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(1,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source1_txt.Position(1) h.edit_monte_loc_range_X_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Location relative to initial simulation Locations for each Simulation Run\nSet to 0 for fixed locations per Simulation Run\n'));
    h.edit_monte_source_loc_range_X(2) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(2,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source2_txt.Position(1) h.edit_monte_loc_range_X_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Location relative to initial simulation Locations for each Simulation Run\nSet to 0 for fixed locations per Simulation Run\n'));
    h.edit_monte_source_loc_range_X(3) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(3,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source3_txt.Position(1) h.edit_monte_loc_range_X_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Location relative to initial simulation Locations for each Simulation Run\nSet to 0 for fixed locations per Simulation Run\n'));
    h.edit_monte_source_loc_StdDev_X = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[h.edit_monte_range_StdDev_txt.Position(1) h.edit_monte_loc_range_X_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Standard Deviation for the Means in the Range Values to the left.\n Set to 0 for fixed Means.'));
    % Y Range
    h.edit_monte_source_loc_range_Y(1) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(1,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source1_txt.Position(1) h.edit_monte_loc_range_Y_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Location relative to initial simulation Locations for each Simulation Run\nSet to 0 for fixed locations per Simulation Run\n'));
    h.edit_monte_source_loc_range_Y(2) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(2,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source2_txt.Position(1) h.edit_monte_loc_range_Y_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Location relative to initial simulation Locations for each Simulation Run\nSet to 0 for fixed locations per Simulation Run\n'));
    h.edit_monte_source_loc_range_Y(3) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(3,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source3_txt.Position(1) h.edit_monte_loc_range_Y_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Location relative to initial simulation Locations for each Simulation Run\nSet to 0 for fixed locations per Simulation Run\n'));
    h.edit_monte_source_loc_StdDev_Y = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[h.edit_monte_range_StdDev_txt.Position(1) h.edit_monte_loc_range_Y_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Standard Deviation for the Means in the Range Values to the left.\n Set to 0 for fixed Means.'));
    % Z Range
    h.edit_monte_source_loc_range_Z(1) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(1,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source1_txt.Position(1) h.edit_monte_loc_range_Z_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Location relative to initial simulation Locations for each Simulation Run\nSet to 0 for fixed locations per Simulation Run\n'));
    h.edit_monte_source_loc_range_Z(2) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(2,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source2_txt.Position(1) h.edit_monte_loc_range_Z_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Location relative to initial simulation Locations for each Simulation Run\nSet to 0 for fixed locations per Simulation Run\n'));
    h.edit_monte_source_loc_range_Z(3) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(3,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source3_txt.Position(1) h.edit_monte_loc_range_Z_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Location relative to initial simulation Locations for each Simulation Run\nSet to 0 for fixed locations per Simulation Run\n'));
    h.edit_monte_source_loc_StdDev_Z = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[h.edit_monte_range_StdDev_txt.Position(1) h.edit_monte_loc_range_Z_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Standard Deviation for the Means in the Range Values to the left.\n Set to 0 for fixed Means.'));
    %% Source Orientation Range - relative to True source orientation
    % Az Range
    h.edit_monte_source_ori_range_Az(1) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(1,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source1_txt.Position(1) h.edit_monte_ori_Az_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Orientations relative to initial simulation Orientations for each Simulation Run\nSet to 0 for fixed Orientations per Simulation Run\n'));
    h.edit_monte_source_ori_range_Az(2) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(2,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source2_txt.Position(1) h.edit_monte_ori_Az_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Orientations relative to initial simulation Orientations for each Simulation Run\nSet to 0 for fixed Orientations per Simulation Run\n'));
    h.edit_monte_source_ori_range_Az(3) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(3,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source3_txt.Position(1) h.edit_monte_ori_Az_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Orientations relative to initial simulation Orientations for each Simulation Run\nSet to 0 for fixed Orientations per Simulation Run\n'));
    h.edit_monte_source_ori_StdDev_Az = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[h.edit_monte_range_StdDev_txt.Position(1) h.edit_monte_ori_Az_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Standard Deviation for the Means in the Range Values to the left.\n Set to 0 for fixed Means.'));
    % El Range
    h.edit_monte_source_ori_range_El(1) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(1,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source1_txt.Position(1) h.edit_monte_ori_El_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Orientations relative to initial simulation Orientations for each Simulation Run\nSet to 0 for fixed Orientations per Simulation Run\n'));
    h.edit_monte_source_ori_range_El(2) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(2,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source2_txt.Position(1) h.edit_monte_ori_El_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Orientations relative to initial simulation Orientations for each Simulation Run\nSet to 0 for fixed Orientations per Simulation Run\n'));
    h.edit_monte_source_ori_range_El(3) = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.src_clr(3,:),'Style','edit','Units','normalize',...
        'Position',[h.edit_monte_source3_txt.Position(1) h.edit_monte_ori_El_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Mean +/- StDev (set in last column) Change in Source Orientations relative to initial simulation Orientations for each Simulation Run\nSet to 0 for fixed Orientations per Simulation Run\n'));
    h.edit_monte_source_ori_StdDev_El = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor','k','Style','edit','Units','normalize',...
        'Position',[h.edit_monte_range_StdDev_txt.Position(1) h.edit_monte_ori_El_range_txt.Position(2) ui_xsize ui_ysize],...
        'FontSize',10,'HorizontalAlignment','center','String','0','Tooltip',sprintf('Standard Deviation for the Means in the Range Values to the left.\n Set to 0 for fixed Means.'));
    %% Listbox to select MEG/EEG to calculate during Mont Carlo
    h.listbox_monte_sens_type_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.chan_clr*.45,'Units','normalize',...
        'Position',[sum(h.edit_monte_source_amp_StdDev.Position([1 3]))+.02 1-.1 .05 .08],...
        'FontSize',10,'HorizontalAlignment','left','String',sprintf('Sensors'),'FontWeight','bold');
    h.listbox_monte_sens_type = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.chan_clr*.45,'Style','listbox','Units','normalize',...
        'Position',[h.listbox_monte_sens_type_txt.Position(1) h.listbox_monte_sens_type_txt.Position(2)-.23 h.listbox_monte_sens_type_txt.Position(3) .225],...
        'FontSize',8,'HorizontalAlignment','center','String',{'MEG' 'EEG'},'Value',1,'Max',3);
    %% Radio Filtering after simulation of all data
    h.radio_monte_filter_flag = uicontrol(h.panel_monte_random_params,'Style','radiobutton', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[h.listbox_monte_sens_type.Position(1) h.listbox_monte_sens_type.Position(2)-.13 .16 .08],...
        'FontSize',10,'HorizontalAlignment','left','String','Filter Sensor & Source Data','Value',0);
    %% Menu: Noise Projection Type 'Synthetic', 'Brain' or 'Real'
    h.menu_monte_synthetic_real_data_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.listbox_monte_sens_type_txt.Position([1 3]))+.01 h.listbox_monte_sens_type_txt.Position(2) .125 .08],...
        'FontSize',10,'HorizontalAlignment','left','String','Sensor Noise Type','FontWeight','bold',...
        'Tooltip',sprintf('Sets Sensor Noise Type\n"Synthetic" parameters based on those set in "Simulate M/EEG" Panel\n"Real Sensor Noise" requires *.mat files to have saved variables: "data [channels x samples x trials]", "srate", "lat"'));
    h.menu_monte_synthetic_real_data = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[h.menu_monte_synthetic_real_data_txt.Position(1) h.menu_monte_synthetic_real_data_txt.Position(2)-.09 h.menu_monte_synthetic_real_data_txt.Position(3) .08],...
        'FontSize',8,'HorizontalAlignment','left','String',{'Synthetic Sensor Noise' 'Synthetic Brain Noise' 'Synthetic MVAR/GAN Noise (not implemented yet)' 'Real Sensor Noise'},...
        'Value',1,'Max',3,'Callback',@menu_monte_synthetic_real_data_Callback,...
        'Tooltip',sprintf('Sets Sensor Noise Type\n"Synthetic" parameters based on those set in "Simulate M/EEG" Panel\n"Real Sensor Noise" requires *.mat files to have saved variables: "data [channels x samples x trials]", "srate", "lat"'));
    %% Menu: "synthetic" or "Real" Source waveforms and Sensor noise
    h.menu_monte_synthetic_real_source_txt = uicontrol(h.panel_monte_random_params,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.menu_monte_synthetic_real_data_txt.Position([1 3]))+.02 h.listbox_monte_sens_type_txt.Position(2) .125 .08],...
        'FontSize',10,'HorizontalAlignment','left','String','Source Type','FontWeight','bold',...
        'Tooltip',sprintf('Sets Source waveforms:\n"Synthetic" simulates source from parameters based on those set in "Simulate M/EEG" Panel\n"Real Sources" requires *.mat files to have variables: "source_data [sources x samples x trials]", "srate", "lat"'));
    h.menu_monte_synthetic_real_source = uicontrol(h.panel_monte_random_params,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','popupmenu','Units','normalize',...
        'Position',[h.menu_monte_synthetic_real_source_txt.Position(1) h.menu_monte_synthetic_real_source_txt.Position(2)-.09 .125 .08],...
        'FontSize',8,'HorizontalAlignment','left','String',{'Synthetic Sources' 'Real Sources'},'Value',1,'Max',3,'Callback',@menu_monte_synthetic_real_data_Callback,...
        'Tooltip',sprintf('Sets Source waveforms:\n"Synthetic" simulates source from parameters based on those set in "Simulate M/EEG" Panel\n"Real Sources" requires *.mat files to have variables: "source_data [sources x samples x trials]", "srate", "lat"'));
    %% Btn Set Real Data Dir
    h.btn_set_real_datadir = uicontrol(h.panel_monte_random_params,'BackgroundColor',[.9 .9 .9],'ForegroundColor','k','Style','pushbutton','Units','normalize',...
        'Position',[h.menu_monte_synthetic_real_data.Position(1) h.menu_monte_synthetic_real_data.Position(2)-.15 .080 .1],...
        'FontSize',10,'HorizontalAlignment','center','String','Set Real Dir','Value',1,'Visible','off','Callback',@btn_set_real_datadir);
    %% Text Real Data Dir
    h.real_datadir_txt = uicontrol(h.panel_monte_random_params,'Style','text','BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.btn_set_real_datadir.Position([1 3]))+.01 h.btn_set_real_datadir.Position(2) .35 .09],...
        'FontSize',10,'HorizontalAlignment','left','String','Real Data Dir: ____________________________','Visible','off');
    
    %% %%%% Panel "Source Modeling Parameters"  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % delete(h.panel_monte_inverse_model);
    h.panel_monte_inverse_model = uipanel(h.tab_monte_carlo,'Title','Source Modeling Parameters','FontSize',10,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[1 .3 .5],...
        'Units','normalize','Position',[.005 h.panel_monte_random_params.Position(2)-.3 .99 .3],'Visible','on');
    %% Listbox to select Inv Soln to calculate during Mont Carlo
    h.listbox_monte_inv_soln_txt = uicontrol(h.panel_monte_inverse_model,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',h.chan_clr*.45,'Units','normalize',...
        'Position',[.01 1-.085 .105 .075],...
        'FontSize',10,'HorizontalAlignment','center','String',sprintf('Inverse Models'),'FontWeight','bold');
    h.listbox_monte_inv_soln = uicontrol(h.panel_monte_inverse_model,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',h.chan_clr*.45,'Style','listbox','Units','normalize',...
        'Position',[h.listbox_monte_inv_soln_txt.Position(1) h.listbox_monte_inv_soln_txt.Position(2)-.75 h.listbox_monte_inv_soln_txt.Position(3) .725],...
        'FontSize',8,'HorizontalAlignment','center','String',{'SPA' 'SIA' 'MIA' 'LCMV' 'eLORETA' 'sLORETA' 'MNE' 'dics' 'pcc' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'},'Value',1,'Max',11);
    %% Listbox: Analyses Type for Source Modeling
    h.listbox_monte_inv_analyses_txt = uicontrol(h.panel_monte_inverse_model,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[sum(h.listbox_monte_inv_soln.Position([1 3]))+.005 h.listbox_monte_inv_soln_txt.Position(2:4)],...
        'FontSize',10,'HorizontalAlignment','center','String',sprintf('Analysis Methods'),'FontWeight','bold');
    h.listbox_monte_inv_analyses = uicontrol(h.panel_monte_inverse_model,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','listbox','Units','normalize',...
        'Position',[h.listbox_monte_inv_analyses_txt.Position(1) h.listbox_monte_inv_analyses_txt.Position(2)-.385 h.listbox_monte_inv_analyses_txt.Position(3) .365],'Enable','on',...
        'FontSize',8,'HorizontalAlignment','center','String',{'none','Seeded Functional Connectivity' 'more to come'},'Max',2,'Value',1,'Callback',@sm_menu_inv_analyses_CallBack);
    
    %     %% Radio to run Functional Connectivity (FC) analyses after source modeling
    %     h.radio_monte_FC_analysis = uicontrol(h.panel_monte_inverse_model,'Style','radiobutton', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
    %         'Position',[h.listbox_monte_inv_soln.Position(1) h.listbox_monte_inv_soln.Position(2)-.085 .175 .075],...
    %         'FontSize',10,'HorizontalAlignment','left','String','Perform Connectivity Analyses','Value',1);
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%%% Button "Export to BST" <-- Not implemented yet because closing SimMEEG call exporting to Brainstorm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    %% Btn Save to Brainstorm if Brainstorm is running
    %     h.btn_output_bst_data = uicontrol(h.main_fig,'BackgroundColor',[1 1 1]*.4,'ForegroundColor',[1 1 1],'Style','pushbutton','Units','normalize',...
    %         'Position',[.88 .84 .1 .025],'Visible','on','UserData',1,'Enable','inactive','visible','off',...
    %         'FontSize',10,'HorizontalAlignment','center','String','Export to BST','Callback',@sm_export_brainstorm_data);
    %     if isappdata(0, 'BrainstormRunning')
    %         h.btn_output_bst_data.Enable = 'on';
    %     else
    %         warning('Brainstorm is not running. You will not be able to export directly into Brainstorm.');
    %     end
    
    %% %%%%% Edit "FontSize Gain" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Edit FontSize Gain - for gaining all fontsizes within the main_fig
    h.edit_fontsize_gain_txt = uicontrol(h.main_fig,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.6 .97 .125 .025],...
        'FontSize',10,'HorizontalAlignment','right','String','FontSize Gain:');
    h.edit_fontsize_gain = uicontrol(h.main_fig,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
        'Position',[sum(h.edit_fontsize_gain_txt.Position([1 3]))+.01 h.edit_fontsize_gain_txt.Position(2) .035 h.edit_fontsize_gain_txt.Position(4)],...
        'FontSize',8,'HorizontalAlignment','center','String','1.0','Callback',@sm_fig_fontsize_gain);
    
    
    %% %%%%% Panel "Wait For" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h.waitfor_panel = uipanel(h.main_fig,'Title','Waiting For Execution','FontSize',12,'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor',[.6 0 1],...
        'Units','normalize','Position',[.005 .005 .99 .99],'FontWeight','bold','Visible','off');
    h.waitfor_txt = uicontrol(h.waitfor_panel,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
        'Position',[.25 .25 .5 .5],'FontSize',14,'HorizontalAlignment','left','String','Wait for Message');
    h.btn_waitfor_close = uicontrol(h.waitfor_panel,'BackgroundColor',[1 .8 .8],'ForegroundColor',[1 0 0],'Style','pushbutton','Units','normalize',...
        'Position',[.975 .975 .025 .025],'Visible','on',...
        'FontSize',10,'HorizontalAlignment','center','String','X','Callback',@hide_waitfor_panel);
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Misc Figure setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Running extra commands
    set_current_ROI([],[],[]); %kl added an extra parameter for curr_tab to be passed
    set_current_source([],[],[]); %kl added an extra parameter for curr_tab to be passed
    update_tfr_roi_cfg([],[]);
    update_monte_carlo;
    h.load_study_flag = 0;  % 1=for loading a saved study; 0=not loading saved study
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZING Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h.sim_data = [];
    
    if exist('bst','var')           % Brainstorm called in SimMEEG directly and thus loading anaotmy directly from Brainstorm
        % bs_anatomy =
        sm_bst2ft_anatomy_from_bst_files(bst);  % converting
    else
        load_default_anatomy();     % loading default anatomy
    end
    initialize_cfg_source_data();
    menu_head_model_CallBack();
    plot_3D_mri();
    %% Turning off plot toolbar
    for a=1:3
        h.ax_power(a).Toolbar.Visible = 'off'; h.ax_PLV(a).Toolbar.Visible = 'off'; h.ax_PLI(a).Toolbar.Visible = 'off';
        h.ax_sig_waves(a).Toolbar.Visible = 'off'; h.ax_prepost_waves(a).Toolbar.Visible = 'off';
        h.ax_sig_plv(a).Toolbar.Visible = 'off'; h.ax_prepost_plv(a).Toolbar.Visible = 'off';
        h.ax_sig_pli(a).Toolbar.Visible = 'off'; h.ax_prepost_pli(a).Toolbar.Visible = 'off';
        disableDefaultInteractivity(h.ax_power(a)); disableDefaultInteractivity(h.ax_PLV(a));  disableDefaultInteractivity(h.ax_PLI(a));
        disableDefaultInteractivity(h.ax_sig_waves(a)); disableDefaultInteractivity(h.ax_prepost_waves(a));
        disableDefaultInteractivity(h.ax_sig_plv(a));   disableDefaultInteractivity(h.ax_prepost_plv(a));
        disableDefaultInteractivity(h.ax_sig_pli(a));   disableDefaultInteractivity(h.ax_prepost_pli(a));
    end
    % update_graphs()
    h.panel_PAC_params.Visible = 'off';
    
    %% Main Figure created
    h.start_flag = 0;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%% Initialize TFR source data %%%%%%%%%%%%%%%
function initialize_cfg_source_data(varargin) %% Initializing h.cfg.source data
global h
h.num_sig_freqs=1;
for v=1:3
    h.cfg.source.sig_freqs(v,h.num_sig_freqs,:)             = [0 0];
    h.cfg.source.sig_amp_perc(v,h.num_sig_freqs)            = 0;
    h.cfg.source.prepost_amp_perc(v,h.num_sig_freqs)        = 0;
    h.cfg.source.sig_amp_perc_std(v,h.num_sig_freqs)        = 0;
    h.cfg.source.prepost_amp_perc_std(v,h.num_sig_freqs)    = 0;
    h.cfg.source.sig_evoked_perc(v,h.num_sig_freqs)         = 0;
    h.cfg.source.prepost_evoked_perc(v,h.num_sig_freqs)     = 0;
    h.cfg.source.sig_durs(v,h.num_sig_freqs)                = 0;
    h.cfg.source.sig_start(v,h.num_sig_freqs)               = 0;
    h.cfg.source.sig_win_type(v,h.num_sig_freqs)            = 0;
    h.cfg.source.sig_win_rise_time(v,h.num_sig_freqs)       = 0;
    
    h.cfg.source.sig_PLV_targets(v,h.num_sig_freqs)         = 0;  %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
    h.cfg.source.prepost_PLV_targets(v,h.num_sig_freqs)     = 0;  %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
    h.cfg.source.sig_PLI_targets(v,h.num_sig_freqs)         = 0;  %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
    h.cfg.source.prepost_PLI_targets(v,h.num_sig_freqs)     = 0;  %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
    h.cfg.source.sig_phase_lag(v,h.num_sig_freqs)           = 0;  % phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval relative to first sample --> (phase_lag/360)*2*pi);    cos(phase_lag) = correlation of signal relative to zero-phase onset
    h.cfg.source.prepost_phase_lag(v,h.num_sig_freqs)       = 0;  % phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval relative to first sample --> (phase_lag/360)*2*pi);    cos(phase_lag) = correlation of signal relative to zero-phase onset
    h.cfg.source.vx_ori(v,:)                                = [0 0 1];  % source orientations (X, Y, Z)
    h.cfg.source.vx_idx(v)                                  = 1;    % source's voxel index from leadfield positions
    h.cfg.source.vx_amp(v)                                  = 30;   % nAmps
    h.cfg.source.vx_locs(v,:)                               = [0 0 0];   % source locations (X, Y, Z)
end
h.cfg.study.source_locs_mm = [-5   -95     5; -55   -10     5; 60   -15     5]; %ones(3,3); %[];

% h.cfg.source.phase_amp_contrasts = []; %[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; % sig_contrasts = cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
% h.cfg.source.sig_phase_amp_freq_idx = []; % [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0]; %(sig_contrasts x Nfreqs);  % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
% h.cfg.source.prepost_phase_amp_freq_idx = []; % [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ];   %(PLV_contrasts x Nfreqs);     % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
% h.cfg.source.sig_phase_amp_depth_perc = []; % [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; % amplitude-modulation depth as a percentage of signal's amplitude (sig_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
% h.cfg.source.prepost_phase_amp_depth_perc = []; % [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; % depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
% h.cfg.source.sig_phase_amp_depth_perc_range = []; % [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; % +/- range of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: sig_phase_amp_depth_perc +/- sig_phase_amp_depth_perc_range must be within [0 100]
% h.cfg.source.prepost_phase_amp_depth_perc_range = []; % [0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; % +/- range devitaion of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: prepost_phase_amp_depth_perc +/- prepost_phase_amp_depth_perc_range must be within [0 100]

%% Initializing PAC matrices
h.cfg.source.phase_amp_contrasts                        = [1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; % These are fixed --> sig_contrasts = cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
h.cfg.source.sig_phase_amp_freq_idx                     = zeros(size(h.cfg.source.phase_amp_contrasts,1), 1); %[0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; %(sig_contrasts x Nfreqs);  % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
h.cfg.source.prepost_phase_amp_freq_idx                 = zeros(size(h.cfg.source.phase_amp_contrasts,1), 1);   %(PLV_contrasts x Nfreqs);     % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
h.cfg.source.sig_phase_amp_depth_perc                   = h.cfg.source.sig_phase_amp_freq_idx; % amplitude-modulation depth as a percentage of signal's amplitude (sig_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
h.cfg.source.prepost_phase_amp_depth_perc               = h.cfg.source.sig_phase_amp_freq_idx; % depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
h.cfg.source.sig_phase_amp_depth_perc_range             = h.cfg.source.sig_phase_amp_freq_idx; % +/- range of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: sig_phase_amp_depth_perc +/- sig_phase_amp_depth_perc_range must be within [0 100]
h.cfg.source.prepost_phase_amp_depth_perc_range         = h.cfg.source.sig_phase_amp_freq_idx; % +/- range devitaion of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: prepost_phase_amp_depth_perc +/- prepost_phase_amp_depth_perc_range must be within [0 100]

%% %%%%% Run Simulation %%%%%%%%%%%%%%%
function run_sim(src,hobj)
global h

h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Running PLV/PLI Simulation.\n\nPlease wait ...'); drawnow;

update_study_cfg(src,hobj);
update_source_cfg(src,hobj);
[h.sim_data.sig_final,h.sim_data.sig_wav,h.sim_data.prepost_wav,h.sim_data.noise_wav,h.sim_data.cfg,h.sim_data.prepost_win,h.sim_data.sig_win] = SimSignals(h.cfg);
% SimSignals_v2 - sets phases for all freqs within an ROI to be same starting phase such that signals line up, BUT not using becuase I haven't fully checked if PLV/PLI are correct after doing this.
% [h.sim_data.sig_final,h.sim_data.sig_wav,h.sim_data.prepost_wav,h.sim_data.noise_wav,h.sim_data.cfg,h.sim_data.prepost_win,h.sim_data.sig_win] = SimSignals_v2(h.cfg);

% display intersource correlations
for t=1:90; r1(:,:,t) = corr(h.sim_data.sig_final(h.cfg.study.act_samps,:,t)); end
r2 = corr(nanmean(h.sim_data.sig_final(h.cfg.study.act_samps,:,:),3));
tb = table(mean(r1,3),r2,'RowNames',{'Source 1' 'Source 2' 'Source 3' },'VariableNames',{'Inter-Source Correlations' 'Inter-Source Correlations Evoked'});
h.sim_data.intersource_correlations = tb;
display(tb);

h.sim_data.source_waveform_type = 'Synthetic Sources';
try h.sim_data = rmfield(h.sim_data,'sig_final_org'); catch; end

h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');

%% %%%%% Save & Load cfg parameters %%%%%%%%%%%%%%%
function save_study(varargin)
global h
[fname,fpath]=uiputfile( '*.mat','Save SimMEEG Study',sprintf('%s_study_parameters',h.edit_study_name.String) );
h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Saving \n\n%s\n\nPlease wait...',fname); drawnow;
cfg = h.cfg;
if isfield (h,'monte_params'); monte_params = h.monte_params; else; monte_params=[]; end

%% Not saving "sim_data" or "inv_data" anymore  -- files are large and inneficient to load back in
% if isfield (h,'sim_data')
%     if isfield(h.sim_data,'sens_final_org')
%         % restoring original data and removing org to save space
%         h.sim_data.sens_final = h.sim_data.sens_final_org;
%         h.sim_data.sens_noise_final = h.sim_data.sens_noise_final_org;
%         h.sim_data.sens_sig_data = h.sim_data.sens_sig_data_org;
%         h.sim_data.sig_final = h.sim_data.sig_final_org;
%
%         %% removing original sens data
%         if isfield(h.sim_data,'sens_final_org')
%             h.sim_data = rmfield(h.sim_data,'sens_final_org');
%             h.sim_data = rmfield(h.sim_data,'sens_noise_final_org');
%             h.sim_data = rmfield(h.sim_data,'sens_sig_data_org');
%             h.sim_data = rmfield(h.sim_data,'sig_final_org');
%         end
%     end
%
%     sim_data = h.sim_data;  sim_data = double2single(sim_data);
% else
%     sim_data=[];
% end
% if isfield (h,'inv_soln'); inv_soln = h.inv_soln; else; inv_soln=[]; end
study_name = h.edit_study_name.String;
% save(fullfile(fpath,fname),'cfg','sim_data','inv_soln','study_name','monte_params');  % saves all figure data to be opened as a new study
% not saving "sim_data" or "inv_data" anymore  -- files are large and inneficient to load back in
save(fullfile(fpath,fname),'cfg','study_name','monte_params');  % saves all figure data to be opened as a new study
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');

% hm.Children(2).Children(1).String(1) = {'File Saved'}; hm.Children(2).Children(1).String(2) = {'Continue'};
h.main_fig.Name = fname;
function load_study(varargin)
global h
answ = questdlg('Are you sure you want to close this study and Load a saved study?','Load SimMEEG Study?','Yes','No','No');
switch answ
    case 'Yes'
        [fname,fpath]=uigetfile('.mat','Load SimMEEG Study');
        %         try
        h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Loading \n\n%s\n\nPlease wait...',fname); drawnow;
        
        if fname~=0
            load(fullfile(fpath,fname)); addToolbarExplorationButtons(gcf); % Adds zoom, rotate, ... buttons to figure toolbar
            
            h.load_study_flag = 1; % loading study
            
            %% clearing plots
            clear_TFR_plots;
            %% resetting cfg and data
            h.cfg = cfg;
            
            %% update study parameters
            h.edit_srate.String = sprintf('%.f',cfg.study.srate);
            h.edit_dur.String = sprintf('%.3f %.3f',cfg.study.dur);
            h.edit_num_trials.String = sprintf('%.f',cfg.study.num_trials);
            h.edit_sens_SNR.String = sprintf('%.1f',cfg.study.SNR);
            h.edit_noise_freqs.String = sprintf('%.1f %.1f',cfg.study.noise_freqs);
            h.edit_noise_amp_perc.String = sprintf('%.f',cfg.study.noise_amp_perc);
            h.edit_noise_flag.Value = cfg.study.noise_flag;
            h.edit_pink_noise_slope.String = sprintf('%.1f',cfg.study.pink_noise_slope);
            h.edit_base_int.String = sprintf('%.3f %.3f',cfg.study.base_int);
            h.edit_poststim_int.String = sprintf('%.3f %.3f',cfg.study.poststim_int);
            h.edit_plv_thresh.String = sprintf('%.2f',cfg.study.plv_thresh);
            h.edit_max_perm_plv.Value = find(strcmp(num2str(cfg.study.max_perm_plv),h.edit_max_perm_plv.String)==1);
            h.edit_plot_time_int.String = sprintf('%.3f %.3f',cfg.study.plot_time_int);
            h.edit_plot_freq_int.String = sprintf('%.1f %.1f',cfg.study.plot_freq_int);
            h.edit_plot_caxis.String = sprintf('%.f %.f',cfg.study.plot_caxis);
            h.edit_study_name.String = cfg.study.study_name;
            
            %% Not loading "sim_data" or "inv_data" anymore  -- files are large and inneficient to load back in
            %             sim_data = single2double(sim_data);
            %             h.inv_soln = inv_soln; h.sim_data = sim_data;
            %             %% update Listbox in Source Modeling Panel
            %             xnames='';
            %             for s=1:length(h.inv_soln); xnames{s} = h.inv_soln(s).ListBox_name; end
            %             h.listbox_inv_solns.String = xnames; h.listbox_inv_solns.Value=1;
            %
            %% clearing "sim_data" and "inv_soln"
            h.sim_data = []; h.inv_soln=[];
            h.listbox_inv_solns.String = ' ';  h.listbox_inv_solns.Value=1; h.current_inv_soln=[];
            
            %% update source locs, amps, and ori
            for v = 1:3
                h.edit_source_locs(v).String = sprintf('%.f %.f %.f',h.cfg.study.source_locs_mm(v,:));
                h.edit_source_amp(v).String = sprintf('%.f',h.cfg.source.vx_amp(v));
                [az,el] = cart2sph(h.cfg.source.vx_ori(v,1),h.cfg.source.vx_ori(v,2),h.cfg.source.vx_ori(v,3));
                h.edit_source_ori(v).String = sprintf('%.f %.f',rad2deg(az),rad2deg(el));
            end
            
            replot_study;
            h.load_study_flag = 0; % not loading study
            h.current_tfr_idx = 1; set_current_ROI('','',1);
            update_source_cfg();
            update_cfg; plot_3D_mri;
            btn_update_PAC_waves();
            
            
            %% update Monte Carlo
            if exist('monte_params','var'); h.monte_params = monte_params; load_monte_carlo_params(); end
            
        end
    case 'No'
end
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
function replot_study(varargin)
global h
replot_tfr_data_from_cfg();
src.Tag = ''; update_cfg(src,[]); % updating cfg based on targets set on screen tab's for ERP, PLV, PLI

if ~isempty(h.tfr_ROI); h.panel_PAC_params.Visible='on'; end
function replot_tfr_data_from_cfg(varargin)  % selects tfr_ROI
% Select a rectangular region of interest using rbbox.m with coordinates of the current axis and draw a box if plot_flag==1;
global h
rotate3d off
ha=h.ax_power(1);   % always Source 1 Power axes
plot_flag=1;

% creating new TFR ROI data for all 3 sources
for tfr_idx = 1:size(h.cfg.source.sig_freqs,2)
    for v = 1:3
        
        p1(1) = h.cfg.source.sig_start(v,tfr_idx);  % signal start
        offset(1) = h.cfg.source.sig_durs(v,tfr_idx);
        p1(2) = h.cfg.source.sig_freqs(v,tfr_idx,1);
        offset(2) = h.cfg.source.sig_freqs(v,tfr_idx,2);
        
        
        x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
        y = [p1(2) p1(2) offset(2) offset(2) p1(2)];
        
        % h.tfr_ROI(v).h(tfr_idx)=plot(x,y,'color',h.src_clr((v),:),'linewidth',2);
        h.tfr_ROI(v).h(tfr_idx)=plot(h.ax_power(v),x,y,'--','color',[1 1 1]*.7,'linewidth',1); %h.tfr_ROI(v).h(tfr_idx).Visible='off';
        h.tfr_ROI(v).h(tfr_idx).UserData.source_num=v;
        h.tfr_ROI(v).h(tfr_idx).UserData.tfr_idx=tfr_idx;
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.xpos=x;
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.ypos=y;
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_amp=[0 0 1 1 0 0]*h.cfg.source.sig_amp_perc(v,tfr_idx)/100;
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_amp=[1 1 0 0 1 1]*h.cfg.source.prepost_amp_perc(v,tfr_idx)/100;
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_plv=[0 0 1 1 0 0]*h.cfg.source.sig_PLV_targets(v,tfr_idx);
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_plv=[1 1 0 0 1 1]*h.cfg.source.prepost_PLV_targets(v,tfr_idx);
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_pli=[0 0 1 1 0 0]*h.cfg.source.sig_PLI_targets(v,tfr_idx);
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_pli=[1 1 0 0 1 1]*h.cfg.source.prepost_PLI_targets(v,tfr_idx);
        %         h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi= [h.cfg.study.lat_sim(1) x(1) x(1)+range(x)/3 x(2)-range(x)/3 x(2) h.cfg.study.lat_sim(end)]; %[h.cfg.study.plot_time_int(1) x(1) x(1)+range(x)/3 x(2)-range(x)/3 x(2) h.cfg.study.plot_time_int(2)];
        
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi = [h.cfg.study.lat_sim(1) x(1) x(1)+h.cfg.source.sig_win_rise_time(v,tfr_idx) x(2)-h.cfg.source.sig_win_rise_time(v,tfr_idx) x(2) h.cfg.study.lat_sim(end)];
        h.edit_tfr_roi_risetime.String = sprintf('%.3f',h.cfg.source.sig_win_rise_time(v,tfr_idx));
        h.edit_tfr_roi_risetime_PLV.String = sprintf('%.3f',h.cfg.source.sig_win_rise_time(v,tfr_idx));
        h.edit_tfr_roi_risetime_PLI.String = sprintf('%.3f',h.cfg.source.sig_win_rise_time(v,tfr_idx));
        
        h.menu_sig_win_type.Value = h.cfg.source.sig_win_type(v,tfr_idx);
        
        h.tfr_ROI(v).h(tfr_idx).ButtonDownFcn=@selected_ROI;
        
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_evk_perc = h.cfg.source.sig_evoked_perc(v,tfr_idx);
        h.edit_sig_evoked_perc(v).String = num2str(h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_evk_perc);
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_evk_perc = h.cfg.source.prepost_evoked_perc(v,tfr_idx);
        h.edit_prepost_evoked_perc(v).String = num2str(h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_evk_perc);
        
        % PLV info
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_phase_start = rad2deg(h.cfg.source.sig_phase_lag(v,tfr_idx));
        h.edit_sig_phase_start(v).String = num2str(h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_phase_start);
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_phase_start = rad2deg(h.cfg.source.prepost_phase_lag(v,tfr_idx));
        h.edit_prepost_phase_start(v).String = num2str(h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_phase_start);
        
        freqs=round(y(1)):round(y(3));
        
        for f=1:length(freqs)
            h.tfr_ROI(v).h(tfr_idx).UserData.waves(f,:)=sin(2*pi*h.cfg.study.lat_sim*freqs(f)+(rand(1)*2*pi));
            h.tfr_ROI(v).h(tfr_idx).UserData.freq_idx(f)=freqs(f);
        end
        h.tfr_ROI(v).h(tfr_idx).UserData.deleted=0;
    end
    
    h.menu_tfr_roi_idx.String = {1:length(h.tfr_ROI(1).h)};
    h.menu_tfr_roi_idx.Value=length(h.tfr_ROI(1).h);
    
    h.menu_tfr_roi_idx_PLV.String = {1:length(h.tfr_ROI(1).h)}; %kl
    h.menu_tfr_roi_idx_PLV.Value=length(h.tfr_ROI(1).h); %kl
    
    h.menu_tfr_roi_idx_PLI.String = {1:length(h.tfr_ROI(1).h)}; %kl
    h.menu_tfr_roi_idx_PLI.Value=length(h.tfr_ROI(1).h); %kl
    h.current_tfr_idx = tfr_idx;
    
    for v=1:3
        h.current_source_num=v; update_tfr_roi_cfg([],[]);
    end
    btn.Button=1; selected_ROI(h.tfr_ROI(v).h(tfr_idx),btn); % automatically plotting window function
    
end

update_source_cfg(); update_PAC();
function new_study(varargin)
global h
answ = questdlg('Are you sure you want to clear this study and start a new SimMEEG study?','New SimMEEG Study?','Yes','No','No');
switch answ
    case 'Yes'
        h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('\nCreating New study ...\n'); drawnow;
        h.new_study_flag = 1;
        clear_TFR_plots;
        initialize_cfg_source_data;
        h.new_study_flag = 0;
        h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
    case 'No'
end
function clear_TFR_plots(varargin)
global h
% clearing TFR and Wave plots
for a=1:3
    h.ax_power(a).clo;
    h.ax_sig_waves(a).clo; h.ax_prepost_waves(a).clo;
    h.ax_sig_plv(a).clo; h.ax_prepost_plv(a).clo;
    h.ax_sig_pli(a).clo; h.ax_prepost_pli(a).clo;
end

%% replotting zeros for Source data
h.tfr_data.s1=surf(h.ax_power(1),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
view(h.ax_power(1),0,90); shading(h.ax_power(1),'flat'); axis(h.ax_power(1),'tight'); axis(h.ax_power(1),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
hold on; h.tfr_data.tfr_zero_line(1) = plot(h.ax_power(1),[0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title(h.ax_power(1),'Source 1 Power','Color',h.src_clr(1,:)); colormap(h.ax_power(1),jet(255)); h.ax_power(1).CLim=h.caxis_power;
h.tfr_data.s2=surf(h.ax_power(2),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
view(h.ax_power(2),0,90); shading(h.ax_power(2),'flat'); axis(h.ax_power(2),'tight'); axis(h.ax_power(2),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
hold on; h.tfr_data.tfr_zero_line(2) = plot(h.ax_power(2),[0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title(h.ax_power(2),'Source 2 Power','Color',h.src_clr(2,:)); colormap(h.ax_power(2),jet(255)); h.ax_power(2).CLim=h.caxis_power;
h.tfr_data.s3=surf(h.ax_power(3),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
view(h.ax_power(3),0,90); shading(h.ax_power(3),'flat'); axis(h.ax_power(3),'tight'); axis(h.ax_power(3),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
hold on; h.tfr_data.tfr_zero_line(3) = plot(h.ax_power(3),[0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title(h.ax_power(3),'Source 3 Power','Color',h.src_clr(3,:)); colormap(h.ax_power(3),jet(255)); h.ax_power(3).CLim=h.caxis_power;
%% replotting zeros for PLV data
h.plv_data.s1=surf(h.ax_PLV(1),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
view(h.ax_PLV(1),0,90); shading(h.ax_PLV(1),'flat'); axis(h.ax_PLV(1),'tight'); axis(h.ax_PLV(1),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
hold on; h.plv_data.tfr_zero_line(1) = plot(h.ax_PLV(1),[0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title(h.ax_PLV(1),'Source 1 vs 2 PLV','Color',h.plv_clr(1,:)); colormap(h.ax_PLV(1),jet(255)); h.ax_PLV(1).CLim=h.caxis_power;
h.plv_data.s2=surf(h.ax_PLV(2),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
view(h.ax_PLV(2),0,90); shading(h.ax_PLV(2),'flat'); axis(h.ax_PLV(2),'tight'); axis(h.ax_PLV(2),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
hold on; h.plv_data.tfr_zero_line(2) = plot(h.ax_PLV(2),[0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title(h.ax_PLV(2),'Source 1 vs 3 PLV','Color',h.plv_clr(2,:)); colormap(h.ax_PLV(2),jet(255)); h.ax_PLV(2).CLim=h.caxis_power;
h.plv_data.s3=surf(h.ax_PLV(3),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
view(h.ax_PLV(3),0,90); shading(h.ax_PLV(3),'flat'); axis(h.ax_PLV(3),'tight'); axis(h.ax_PLV(3),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
hold on; h.plv_data.tfr_zero_line(3) = plot(h.ax_PLV(3),[0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title(h.ax_PLV(3),'Source 2 vs 3 PLV','Color',h.plv_clr(3,:)); colormap(h.ax_PLV(3),jet(255)); h.ax_PLV(3).CLim=h.caxis_power;
%% replotting zeros for PLI data
h.pli_data.s1=surf(h.ax_PLI(1),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
view(h.ax_PLI(1),0,90); shading(h.ax_PLI(1),'flat'); axis(h.ax_PLI(1),'tight'); axis(h.ax_PLI(1),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
hold on; h.pli_data.tfr_zero_line(1) = plot(h.ax_PLI(1),[0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title(h.ax_PLI(1),'Source 1 vs 2 PLI/dPLI','Color',h.plv_clr(1,:)); colormap(h.ax_PLI(1),jet(255)); h.ax_PLI(1).CLim=h.caxis_power;
h.pli_data.s2=surf(h.ax_PLI(2),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
view(h.ax_PLI(2),0,90); shading(h.ax_PLI(2),'flat'); axis(h.ax_PLI(2),'tight'); axis(h.ax_PLI(2),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
hold on; h.pli_data.tfr_zero_line(2) = plot(h.ax_PLI(2),[0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title(h.ax_PLI(2),'Source 1 vs 3 PLI/dPLI','Color',h.plv_clr(2,:)); colormap(h.ax_PLI(2),jet(255)); h.ax_PLI(2).CLim=h.caxis_power;
h.pli_data.s3=surf(h.ax_PLI(3),h.cfg.study.lat_sim,h.freq_vals,zeros(length(h.freq_vals),length(h.cfg.study.lat_sim)));
view(h.ax_PLI(3),0,90); shading(h.ax_PLI(3),'flat'); axis(h.ax_PLI(3),'tight'); axis(h.ax_PLI(3),[h.cfg.study.plot_time_vals([1 end]) h.cfg.study.plot_freq_vals([1 end])]);
hold on; h.pli_data.tfr_zero_line(3) = plot(h.ax_PLI(3),[0 0],[h.cfg.study.plot_freq_int],'k--'); box on; title(h.ax_PLI(3),'Source 2 vs 3 PLI/dPLI','Color',h.plv_clr(3,:)); colormap(h.ax_PLI(3),jet(255)); h.ax_PLI(3).CLim=h.caxis_power;

%% clearing tfr_ROI data
h.tfr_ROI=[];
%% clearing sim_data
h.sim_data = [];
h.axes_sens_noise.clo;
h.axes_sens_signal.clo;
h.axes_sens_final.clo;

%% updating edit boxes
for a=1:3
    % Signal
    h.edit_sig_power_perc(a).String     = '100';
    h.edit_sig_evoked_perc(a).String    = '50';
    h.edit_sig_phase_locking(a).String  = '0.4';
    h.edit_sig_phase_lag(a).String      = '0.2';
    h.edit_sig_phase_start(a).String    = '0';
    % Prepost
    h.edit_prepost_power_perc(a).String     = '100';
    h.edit_prepost_evoked_perc(a).String    = '0';
    h.edit_prepost_phase_locking(a).String  = '0';
    h.edit_prepost_phase_lag(a).String      = '0';
    h.edit_prepost_phase_start(a).String    = '0';
end
%% clearing menus
h.menu_tfr_roi_idx.Value=1; h.menu_tfr_roi_idx.String='none'; h.menu_tfr_roi_idx.Visible='off';
h.menu_tfr_roi_idx_PLV.Value=1; h.menu_tfr_roi_idx_PLV.String='none'; h.menu_tfr_roi_idx_PLV.Visible='off';
h.menu_tfr_roi_idx_PLI.Value=1; h.menu_tfr_roi_idx_PLI.String='none'; h.menu_tfr_roi_idx_PLI.Visible='off';
h.menu_source_idx.Value=1; h.menu_source_idx.Visible='off';
h.menu_source_idx_PLV.Value=1; h.menu_source_idx_PLV.Visible='off';
h.menu_source_idx_PLI.Value=1; h.menu_source_idx_PLI.Visible='off';
h.edit_tfr_roi_risetime.Visible = 'off'; h.edit_tfr_roi_risetime_PLV.Visible = 'off'; h.edit_tfr_roi_risetime_PLI.Visible = 'off';
h.edit_tfr_roi_xlim.Visible='off'; h.edit_tfr_roi_ylim.Visible='off';
h.edit_tfr_roi_xlim_PLV.Visible='off'; h.edit_tfr_roi_ylim_PLV.Visible='off';
h.edit_tfr_roi_xlim_PLI.Visible='off'; h.edit_tfr_roi_ylim_PLI.Visible='off';
%% resetting PAC
for a=1:6
    % Signal
    h.menu_PAC_sig_carrier_source_freq(a).Value=1;
    h.menu_PAC_sig_carrier_source_freq(a).String = {'none'};
    h.menu_PAC_sig_modulator_source_freq(a).Value=1;
    h.menu_PAC_sig_modulator_source_freq(a).String = {'none'};
    h.edit_PAC_modulator_sig_depth(a).String = {'0'};
    h.edit_PAC_modulator_sig_depth_range(a).String = {'0'};
    % Prepost
    h.menu_PAC_prepost_carrier_source_freq(a).Value=1;
    h.menu_PAC_prepost_carrier_source_freq(a).String = {'none'};
    h.menu_PAC_prepost_modulator_source_freq(a).Value=1;
    h.menu_PAC_prepost_modulator_source_freq(a).String = {'none'};
    h.edit_PAC_modulator_prepost_depth(a).String = {'0'};
    h.edit_PAC_modulator_prepost_depth_range(a).String = {'0'};
    h.axes_PAC_waves_sig(a).clo;
    h.axes_PAC_waves_prepost(a).clo;
end
h.panel_PAC_params.Visible = 'on';
%% Clearing Source Modeling
h.listbox_inv_solns.String = {''}; h.listbox_inv_solns.Value=1;
h.listbox_peaks_found.String={''}; h.listbox_peaks_found.Value=1;
h.axes_3D_images.clo;
h.axes_source_fft.clo;
h.axes_source_waves.clo;
h.axes_invSoln_errors_locs.clo;
h.axes_invSoln_errors_ori.clo;
h.axes_invSoln_errors_waves.clo;
%% replot default sources
h.cfg.study.source_locs_mm = [-5   -95     5; -55   -10     5; 60   -15     5]; %ones(3,3); %[];
for v = 1:3
    h.edit_source_locs(v).String = sprintf('%.f %.f %.f',h.cfg.study.source_locs_mm(v,:));
    h.edit_source_amp(v).String = sprintf('%.f',30);
    h.edit_source_ori(v).String = sprintf('%.f %.f',0,90);
end
update_source_data;
plot_3D_mri;
%% update cfg
update_cfg();
function load_monte_carlo_params(varargin)
global h

h.edit_monte_num_sims.String = h.monte_params.num_sims;
for v=1:3
    h.edit_monte_source_amp_range(v).String= num2str(h.monte_params.source_amp_range(v,:));
    h.edit_monte_source_loc_range_X(v).String = num2str(h.monte_params.source_loc_range_X(v,:));
    h.edit_monte_source_loc_range_Y(v).String = num2str(h.monte_params.source_loc_range_Y(v,:));
    h.edit_monte_source_loc_range_Z(v).String = num2str(h.monte_params.source_loc_range_Z(v,:));
    h.edit_monte_source_ori_range_Az(v).String = num2str(h.monte_params.source_ori_range_Az(v,:));
    h.edit_monte_source_ori_range_El(v).String = num2str(h.monte_params.source_ori_range_El(v,:));
end
h.edit_monte_source_amp_StdDev.String = h.monte_params.source_amp_StdDev;
h.edit_monte_source_loc_StdDev_X.String = h.monte_params.source_loc_StdDev_X;
h.edit_monte_source_loc_StdDev_Y.String = h.monte_params.source_loc_StdDev_Y;
h.edit_monte_source_loc_StdDev_Z.String = h.monte_params.source_loc_StdDev_Z;
h.edit_monte_source_ori_StdDev_Az.String = num2str(h.monte_params.source_ori_StdDev_Az);
h.edit_monte_source_ori_StdDev_El.String = num2str(h.monte_params.source_ori_StdDev_El);
h.edit_monte_SNR_range.String = num2str(h.monte_params.SNR_range); h.edit_monte_SNR_StdDev.String = num2str(h.monte_params.SNR_StdDev);
h.edit_monte_plv_range.String = num2str(h.monte_params.plv_range); h.edit_monte_plv_StdDev.String = num2str(h.monte_params.plv_StdDev);
h.listbox_monte_inv_soln.Value = h.monte_params.inv_soln;
h.radio_monte_FC_analysis.Value = h.monte_params.FC_analysis;
h.radio_monte_save_FT_data.Value = h.monte_params.save_FT_data;
h.listbox_monte_sens_type.Value = h.monte_params.sens_type;
function get_datadir(varargin)
global h
data_dir = h.data_dir;
h.data_dir = uigetdir(h.data_dir,'Set Data Directory');
if h.data_dir == 0  % user cancelled setting data directory
    h.data_dir = data_dir;
end

%% %%%%% Save & Load SimMEEG datasets
function sm_load_SimMEEG_dataset(varargin)
% load SimMEEG dataset that was stored using
global h
answ = questdlg(sprintf('This will close the current study\nand load a saved study and dataset.\n\nWould you like to CONTINUE?'),'Load SimMEEG Dataset?','Yes','No','No');
switch answ
    case 'Yes'
        [fname,fpath]=uigetfile('.mat','Load SimMEEG Dataset');
        %         try
        h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Loading \n\n%s\n\nPlease wait...',fname); drawnow;
        
        if any(fname~=0)
            load(fullfile(fpath,fname),'cfg','monte_params','sim_data','inv_soln');
            
            h.load_study_flag = 1; % loading study
            
            %% clearing plots
            clear_TFR_plots;
            %% resetting cfg and data
            cfg=sim_data.cfg;
            h.cfg = cfg;
             
            %% update study parameters
            h.edit_srate.String = sprintf('%.f',cfg.study.srate);
            h.edit_dur.String = sprintf('%.3f %.3f',cfg.study.dur);
            h.edit_num_trials.String = sprintf('%.f',cfg.study.num_trials);
            h.edit_sens_SNR.String = sprintf('%.1f',cfg.study.SNR);
            h.edit_noise_freqs.String = sprintf('%.1f %.1f',cfg.study.noise_freqs);
            h.edit_noise_amp_perc.String = sprintf('%.f',cfg.study.noise_amp_perc);
            h.edit_noise_flag.Value = cfg.study.noise_flag;
            h.edit_pink_noise_slope.String = sprintf('%.1f',cfg.study.pink_noise_slope);
            h.edit_base_int.String = sprintf('%.3f %.3f',cfg.study.base_int);
            h.edit_poststim_int.String = sprintf('%.3f %.3f',cfg.study.poststim_int);
            h.edit_plv_thresh.String = sprintf('%.2f',cfg.study.plv_thresh);
            h.edit_max_perm_plv.Value = find(strcmp(num2str(cfg.study.max_perm_plv),h.edit_max_perm_plv.String)==1);
            h.edit_plot_time_int.String = sprintf('%.3f %.3f',cfg.study.plot_time_int);
            h.edit_plot_freq_int.String = sprintf('%.1f %.1f',cfg.study.plot_freq_int);
            h.edit_plot_caxis.String = sprintf('%.f %.f',cfg.study.plot_caxis);
            h.edit_study_name.String = cfg.study.study_name;
            
            %% clearing "sim_data" and "inv_soln"
            h.sim_data = []; h.inv_soln=[];
            h.listbox_inv_solns.String = ' ';  h.listbox_inv_solns.Value=1; h.current_inv_soln=[];
            
            %% loading sim_data, inv_soln, ... if they exist
            sim_data = single2double(sim_data); h.sim_data = sim_data;
            if exist('inv_soln','var')
                h.inv_soln = inv_soln;
                % update Listbox in Source Modeling Panel
                xnames='';
                for s=1:length(h.inv_soln); xnames{s} = h.inv_soln(s).ListBox_name; end
                h.listbox_inv_solns.String = xnames; h.listbox_inv_solns.Value=1;
            end
            
            
            %% update source locs, amps, and ori
            for v = 1:3
                h.edit_source_locs(v).String = sprintf('%.f %.f %.f',h.cfg.study.source_locs_mm(v,:));
                h.edit_source_amp(v).String = sprintf('%.f',h.cfg.source.vx_amp(v));
                [az,el] = cart2sph(h.cfg.source.vx_ori(v,1),h.cfg.source.vx_ori(v,2),h.cfg.source.vx_ori(v,3));
                h.edit_source_ori(v).String = sprintf('%.f %.f',rad2deg(az),rad2deg(el));
            end
            
            %% replotting study
            replot_study;
            h.load_study_flag = 0; % not loading study
            h.current_tfr_idx = 1; set_current_ROI('','',1);
            update_source_cfg();
            update_cfg; plot_3D_mri;
            btn_update_PAC_waves();
            
            
            %% update Monte Carlo
            if exist('monte_params','var'); h.monte_params = monte_params; load_monte_carlo_params(); end
            h.monte_carlo_flag=0;
        end
    case 'No'
end
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');

function sm_save_SimMEEG_dataset(varargin)
global h

% get list of SimMEEG dataset files within the datadir
[fname,fpath]=uiputfile( '*.mat','Save SimMEEG DataSet',sprintf('%s_SimMEEG_Dataset.mat',h.edit_study_name.String) );
if fname~=0
    h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Saving \n\n%s\n\nPlease wait...',fname); drawnow;
    if exist(fullfile(fpath,fname),'file')
        answ = questdlg(sprintf('File Exists\n\nWould you like to overwrite it?\n'),'Save SimMEEG Dataset?','Yes','No','No');
    else
        answ = 'Yes';
    end
    switch answ
        case 'Yes'
            cfg = h.cfg;
            if isfield (h,'monte_params'); monte_params = h.monte_params; else; monte_params=[]; end
            
            % Not saving "sim_data" or "inv_data" anymore  -- files are large and inneficient to load back in
            if isfield(h,'sim_data')
                if isfield(h.sim_data,'sens_final_org')
                    % restoring original data and removing org to save space
                    h.sim_data.sens_final = h.sim_data.sens_final_org;
                    h.sim_data.sens_noise_final = h.sim_data.sens_noise_final_org;
                    h.sim_data.sens_sig_data = h.sim_data.sens_sig_data_org;
                    h.sim_data.sig_final = h.sim_data.sig_final_org;
                    
                    %% removing original sens data
                    if isfield(h.sim_data,'sens_final_org')
                        h.sim_data = rmfield(h.sim_data,'sens_final_org');
                        h.sim_data = rmfield(h.sim_data,'sens_noise_final_org');
                        h.sim_data = rmfield(h.sim_data,'sens_sig_data_org');
                        h.sim_data = rmfield(h.sim_data,'sig_final_org');
                    end
                end
                
                sim_data = h.sim_data;  sim_data = double2single(sim_data);
            else
                sim_data=[];
            end
            
            if isfield (h,'inv_soln'); inv_soln = h.inv_soln; else; inv_soln=[]; end
            study_name = h.edit_study_name.String;
            save(fullfile(fpath,fname),'cfg','sim_data','inv_soln','study_name','monte_params');  % saves all figure data to be opened as a new study
            h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
            
            % hm.Children(2).Children(1).String(1) = {'File Saved'}; hm.Children(2).Children(1).String(2) = {'Continue'};
            h.main_fig.Name = fname;
        case 'No'
    end
end



%% %%%% Updating cfg before simulate data
function update_cfg(src,hobj)
global h
try
    h.cfg.study.SNR = str2num(h.edit_sens_SNR.String);
    update_study_cfg(src,hobj);
    %     h.cfg.source=[];
    current_tfr_idx = h.current_tfr_idx;
    for r=1:length(h.tfr_ROI(1).h)
        h.current_tfr_idx = r;
        update_source_cfg(src,hobj);
    end
    h.current_tfr_idx = current_tfr_idx;
    update_plv_pli();
    for v=1:3
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_evk_perc = str2num(h.edit_sig_evoked_perc(v).String);
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_evk_perc = str2num(h.edit_prepost_evoked_perc(v).String);
        
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_phase_start = str2num(h.edit_sig_phase_start(v).String);
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_phase_start = str2num(h.edit_prepost_phase_start(v).String);
        
    end
    
    plot_source_tfr(); plot_source_waves();
    update_monte_carlo();
catch me
    if h.new_study_flag==0 && h.start_flag==0 && h.load_study_flag==0
        fprintf('Error in "update_cfg"\n%s\n',me.message);
        %         errordlg('Please make sure a ROI is present');
    end
end
function update_study_cfg(src,~)
global h

h.cfg.study.study_name = h.edit_study_name.String;

if exist('src','var')
    if strcmp(src.Tag,'edit_srate') || strcmp(src.Tag,'edit_dur')   % these will reset study back to beginning clear all plots
        
        ButtonName = questdlg('All ROIs will be deleted. Would you like to continue?', ...
            'Warning!', ...
            'Yes', 'No','No');
        switch ButtonName
            case 'Yes'
                if isfield(h,'tfr_ROI')
                    try
                        while length(h.tfr_ROI(1).h)>=1
                            delete_tfr_ROI;
                        end
                    catch
                    end
                end
                
                h.cfg.study.srate=str2num(h.edit_srate.String);
                h.cfg.study.nyq = h.cfg.study.srate/2;
                h.cfg.study.dur=str2num(h.edit_dur.String); % (sec) start and end times for whole trial
                h.cfg.study.lat_sim=[h.cfg.study.dur(1):1/h.cfg.study.srate:h.cfg.study.dur(2)-(1/h.cfg.study.srate)]; % latency of each trial
                h.cfg.study.num_samps=length(h.cfg.study.lat_sim);
                h.cfg.study.pink_noise_slope = normrnd(1.8,.1,1); % Arbitrarily set to be randomized with overall mean 1.8  % seemed to best match the slope from a couple of resting-state data from LetterAll study
                h.edit_pink_noise_slope.String = sprintf('%.1f',h.cfg.study.pink_noise_slope);
                if h.cfg.study.pink_noise_slope<1; h.cfg.study.pink_noise_slope=1; elseif h.cfg.study.pink_noise_slope>2; h.cfg.study.pink_noise_slope=2; end
                
                h.cfg.study.num_trials = str2num(h.edit_num_trials.String); %{h.edit_num_trials.Value});
                xr=h.cfg.study.num_trials/h.cfg.study.max_perm_plv;
                if mod(xr,1)~=0
                    h.cfg.study.num_trials = h.cfg.study.max_perm_plv*round(xr);
                    h.edit_num_trials.String = num2str(h.cfg.study.num_trials);
                    msgbox(sprintf('"Number Trials" must be integer multiple of "PLV Permutations"\n\nNumber Trials = %.f',h.cfg.study.num_trials));
                end
                % chekcing plot frequencies are witihin study parameter limits
                pfreqs = str2num(h.edit_plot_freq_int.String);
                if length(pfreqs)~=2; pfreqs = [1 h.cfg.study.nyq]; end
                if pfreqs(1)<=0; pfreqs(1) = 1; end
                if pfreqs(2) > h.cfg.study.nyq
                    w=warndlg(sprintf('\nPlot Frequency %.f Hz exceeds Nyquist Frequency %.f\n\n Automatically adjusted to Nyquist',pfreqs(2),h.cfg.study.nyq));
                    w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
                    pfreqs(2) = str2num(h.edit_srate.String)/2;
                end
                h.freq_vals = pfreqs(1):pfreqs(2);
                h.edit_plot_freq_int.String = sprintf('%.f %.f',pfreqs);
                
                h.cfg.study.noise_flag=h.edit_noise_flag.Value;  % Whitening noise to be added to each sources
                h.cfg.study.noise_amp_perc=str2num(h.edit_noise_amp_perc.String); % percent of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
                h.cfg.study.noise_freqs=str2num(h.edit_noise_freqs.String); % [1 100] Frequency (Hz) of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
                h.cfg.study.plv_thresh=str2num(h.edit_plv_thresh.String);   % stoppping criterion when search for best PLV/PLI matched to sig_PLV_targets, sig_PLI_targets, etc. (e.g., 0.05).
                h.cfg.study.plot_sim_flag=0; %h.radio_plot_sim_flag.Value;  % plots evoked, trial, and PLV results. Note: This can take time because of filtering.
                h.cfg.study.plot_time_int=str2num(h.edit_plot_time_int.String);   % time interval to plot
                h.cfg.study.plot_freq_int=str2num(h.edit_plot_freq_int.String); % frequencies for calculating and plotting PLV/PLI
                h.cfg.study.base_int=str2num(h.edit_base_int.String);    % base line interval for plotting
                ss = find(h.cfg.study.lat_sim<=h.cfg.study.base_int(1)); bs1=ss(end);
                ss = find(h.cfg.study.lat_sim<=h.cfg.study.base_int(2)); bs2=ss(end);
                h.cfg.study.base_samps = bs1:bs2;
                
                h.cfg.study.poststim_int=str2num(h.edit_poststim_int.String);    % base line interval for plotting
                
                h.cfg.study.plot_time_vals=h.cfg.study.plot_time_int(1):1/h.cfg.study.srate:h.cfg.study.plot_time_int(2);
                h.cfg.study.plot_freq_vals=h.cfg.study.plot_freq_int(1):h.cfg.study.plot_freq_int(2);
                update_graphs();
                
                % clear all tfr_ROI and plots
                clear_TFR_plots;
                
            case 'No'
                
        end % switch
    else
        
        h.cfg.study.srate=str2num(h.edit_srate.String);
        h.cfg.study.nyq = h.cfg.study.srate/2;
        h.cfg.study.dur=str2num(h.edit_dur.String); % (sec) start and end times for whole trial
        h.cfg.study.lat_sim=[h.cfg.study.dur(1):1/h.cfg.study.srate:h.cfg.study.dur(2)-(1/h.cfg.study.srate)]; % latency of each trial
        h.cfg.study.num_samps=length(h.cfg.study.lat_sim);
        h.cfg.study.pink_noise_slope = normrnd(1.8,.1,1); % Arbitrarily set to be randomized with overall mean 1.8  % seemed to best match the slope from a couple of resting-state data from LetterAll study
        if h.cfg.study.pink_noise_slope<1; h.cfg.study.pink_noise_slope=1; elseif h.cfg.study.pink_noise_slope>2; h.cfg.study.pink_noise_slope=2; end
        h.edit_pink_noise_slope.String = sprintf('%.1f',h.cfg.study.pink_noise_slope);
        h.cfg.study.max_perm_plv = str2num(h.edit_max_perm_plv.String{h.edit_max_perm_plv.Value}); %h.cfg.study.num_trials/10; % integer of cfg.study.num_trials/10 = maximum number of permutations to search for PLV_trials (must be <10 or memory will fail on most computers) .
        h.cfg.study.num_trials = str2num(h.edit_num_trials.String); %{h.edit_num_trials.Value});
        xr=h.cfg.study.num_trials/h.cfg.study.max_perm_plv;
        if mod(xr,1)~=0
            h.cfg.study.num_trials = h.cfg.study.max_perm_plv*round(xr);
            h.edit_num_trials.String = num2str(h.cfg.study.num_trials);
            msgbox(sprintf('"Number Trials" must be integer multiple of "PLV Permutations"\n\nNumber Trials = %.f',h.cfg.study.num_trials));
        end
        % chekcing plot frequencies are witihin study parameter limits
        pfreqs = str2num(h.edit_plot_freq_int.String);
        if length(pfreqs)~=2; pfreqs = [1 h.cfg.study.nyq]; end
        if pfreqs(1)<=0; pfreqs(1) = 1; end
        if pfreqs(2) > h.cfg.study.nyq
            w=warndlg(sprintf('\nPlot Frequency %.f Hz exceeds Nyquist Frequency %.f\n\n Automatically adjusted to Nyquist',pfreqs(2),h.cfg.study.nyq));
            w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
            pfreqs(2) = str2num(h.edit_srate.String)/2;
        end
        h.freq_vals = pfreqs(1):pfreqs(2);
        h.edit_plot_freq_int.String = sprintf('%.f %.f',pfreqs);
        
        h.cfg.study.noise_flag=h.edit_noise_flag.Value;  % Whitening noise to be added to each sources
        h.cfg.study.noise_amp_perc=str2num(h.edit_noise_amp_perc.String); % percent of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
        h.cfg.study.noise_freqs=str2num(h.edit_noise_freqs.String); % [1 100] Frequency (Hz) of noise to add to overall signal throughout the num_samps to whiten the data for time-freq analyses.
        h.cfg.study.plv_thresh=str2num(h.edit_plv_thresh.String);   % stoppping criterion when search for best PLV/PLI matched to sig_PLV_targets, sig_PLI_targets, etc. (e.g., 0.05).
        h.cfg.study.plot_sim_flag=0; %h.radio_plot_sim_flag.Value;  % plots evoked, trial, and PLV results. Note: This can take time because of filtering.
        h.cfg.study.plot_time_int=str2num(h.edit_plot_time_int.String);   % time interval to plot
        h.cfg.study.plot_freq_int=str2num(h.edit_plot_freq_int.String); % frequencies for calculating and plotting PLV/PLI
        h.cfg.study.base_int=str2num(h.edit_base_int.String);    % base line interval for plotting
        ss = find(h.cfg.study.lat_sim<=h.cfg.study.base_int(1)); bs1=ss(end);
        ss = find(h.cfg.study.lat_sim<=h.cfg.study.base_int(2)); bs2=ss(end);
        h.cfg.study.base_samps = bs1:bs2;
        h.cfg.study.poststim_int=str2num(h.edit_poststim_int.String);    % base line interval for plotting
        h.cfg.study.act_samps = (round( (h.cfg.study.poststim_int(1)-h.cfg.study.lat_sim(1))*h.cfg.study.srate ):round( (h.cfg.study.poststim_int(2)-h.cfg.study.lat_sim(1))*h.cfg.study.srate )) +1;
        h.cfg.study.plot_caxis = str2num(h.edit_plot_caxis.String);    % color scale for tfr plots
        
        h.cfg.study.plot_time_vals=h.cfg.study.plot_time_int(1):1/h.cfg.study.srate:h.cfg.study.plot_time_int(2);
        h.cfg.study.plot_freq_vals=h.cfg.study.plot_freq_int(1):h.cfg.study.plot_freq_int(2);
        update_graphs();
    end
end
function update_source_cfg(varargin)
global h

%% Extracting Source Info from selected tfr_ROI for each source
% h.num_sig_freqs=0;
try
    if isfield(h,'tfr_ROI')
        for v=1:length(h.tfr_ROI)   % sources
            for r = h.current_tfr_idx  %1:length(h.tfr_ROI(v).h)  % ROIs per source
                if isvalid(h.tfr_ROI(v).h(r))
                    %                 for f=1:length(h.tfr_ROI(v).h(r).UserData.freq_idx)
                    % TFR data from windowing function
                    %                     fidx(f)=find(h.tfr_data.s1.YData==h.tfr_ROI(v).h(r).UserData.freq_idx(f));
                    %                     tfr_data1(fidx(f),:,r,1,v)=h.tfr_ROI(v).h(r).UserData.roi.sig_win_final;     % sig_win
                    %                     tfr_data1(fidx(f),:,r,2,v)=h.tfr_ROI(v).h(r).UserData.roi.prepost_win_final;     % prepost_win
                    
                    % ROI freqs
                    h.num_sig_freqs = r; % h.num_sig_freqs+1;
                    h.cfg.source.sig_freqs(v,h.num_sig_freqs,:)             = [min(h.tfr_ROI(v).h(r).UserData.freq_idx) max(h.tfr_ROI(v).h(r).UserData.freq_idx)];
                    h.cfg.source.sig_amp_perc(v,h.num_sig_freqs)            = h.tfr_ROI(v).h(r).UserData.roi.sig_amp(3)*100; % max(h.tfr_ROI(v).h(r).UserData.sig_waves(f).YData)*100;
                    h.cfg.source.prepost_amp_perc(v,h.num_sig_freqs)        = h.tfr_ROI(v).h(r).UserData.roi.prepost_amp(2)*100; % max(h.tfr_ROI(v).h(r).UserData.prepost_waves(f).YData)*100;
                    h.cfg.source.sig_amp_perc_std(v,h.num_sig_freqs)        = 0;
                    h.cfg.source.prepost_amp_perc_std(v,h.num_sig_freqs)    = 0;
                    h.cfg.source.sig_evoked_perc(v,h.num_sig_freqs)         = h.tfr_ROI(v).h(r).UserData.roi.sig_evk_perc;
                    h.cfg.source.prepost_evoked_perc(v,h.num_sig_freqs)     = h.tfr_ROI(v).h(r).UserData.roi.prepost_evk_perc;
                    h.cfg.source.sig_durs(v,h.num_sig_freqs)                = h.tfr_ROI(v).h(r).XData(2)-h.tfr_ROI(v).h(r).XData(1);
                    h.cfg.source.sig_start(v,h.num_sig_freqs)               = h.tfr_ROI(v).h(r).XData(1);
                    h.cfg.source.sig_win_type(v,h.num_sig_freqs)            = h.menu_sig_win_type.Value;  %(3signals x Nfreqs) % type of windowing function (1)='Hann' (2)='Gauss' (3)='Triang'  (4)='Blackman' ;
                    h.cfg.source.sig_win_rise_time(v,h.num_sig_freqs)       = h.tfr_ROI(v).h(r).UserData.roi.x_roi(3)- h.tfr_ROI(v).h(r).UserData.roi.x_roi(2) ; %(3signals x Nfreqs) % rise_time of windowing function. If less than 1/2 sig_duration then window function will have a plateau of the difference in duration between the rise/fall time and the signal duration;
                    
                    h.cfg.source.sig_PLV_targets(v,h.num_sig_freqs)         = h.tfr_ROI(v).h(r).UserData.roi.sig_plv(3);  %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
                    h.cfg.source.prepost_PLV_targets(v,h.num_sig_freqs)     = h.tfr_ROI(v).h(r).UserData.roi.prepost_plv(1);  %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
                    h.cfg.source.sig_PLI_targets(v,h.num_sig_freqs)         = h.tfr_ROI(v).h(r).UserData.roi.sig_pli(3);  %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
                    h.cfg.source.prepost_PLI_targets(v,h.num_sig_freqs)     = h.tfr_ROI(v).h(r).UserData.roi.prepost_pli(1); %(PLV_contrasts x Nfreqs);     % for each signal contrast (1-2, 1-3, 2-3) x Nfreqs.
                    
                    h.cfg.source.sig_phase_lag(v,h.num_sig_freqs)           = deg2rad(str2num(h.edit_sig_phase_start(v).String)); %([0 0; 0 0; 0 0]/360)*2*pi;  % phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval relative to first sample --> (phase_lag/360)*2*pi);    cos(phase_lag) = correlation of signal relative to zero-phase onset
                    h.cfg.source.prepost_phase_lag(v,h.num_sig_freqs)       = deg2rad(str2num(h.edit_prepost_phase_start(v).String)); %([0 0; 0 0; 0 0]/360)*2*pi;  % phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval relative to first sample --> (phase_lag/360)*2*pi);    cos(phase_lag) = correlation of signal relative to zero-phase onset
                    
                    update_PAC();
                    
                    % Source info
                    h.cfg.source.vx_locs(v,:)                               = str2num(h.edit_source_locs(v).String);  % source locations X, Y, Z (mm)
                    az_el = deg2rad(str2num(h.edit_source_ori(v).String)); [x,y,z] = sph2cart(az_el(1),az_el(2),1);
                    h.cfg.source.vx_ori(v,:)                                = [x y z];  % source orientations (X, Y, Z)
                    h.cfg.source.vx_idx(v)                                  = find_nearest_voxel(h.cfg.source.vx_locs(v,:),h.anatomy.leadfield.voxel_pos);    % source's voxel index from leadfield positions
                    h.cfg.source.vx_amp(v)                                  = str2num(h.edit_source_amp(v).String);   % nAmps
                    h.cfg.source.vx_locs(v,:)                               = h.anatomy.leadfield.voxel_pos(h.cfg.source.vx_idx(v),:);
                    %                 end
                    %                 % ROI freqs
                    %                 h.cfg.source.sig_freqs(v,:,r)=[min(h.tfr_ROI(v).h(r).UserData.freq_idx) max(h.tfr_ROI(v).h(r).UserData.freq_idx)];
                end
            end
        end
    end
catch me
    h.menu_tfr_roi_idx.Value=1; h.menu_tfr_roi_idx_PLV.Value=1; h.menu_tfr_roi_idx_PLI.Value=1;
    fprintf('Error in "update_source_cfg"\n%s\n',me.message);
end

%% create tfr ROIs
function bl_rbbox_gca_v2(varargin)  % selects tfr_ROI
% Select a rectangular region of interest using rbbox.m with coordinates of the current axis and draw a box if plot_flag==1;
global h
h.panel_PAC_params.Visible = 'on';
rotate3d off
ha=h.ax_power(1);   % always Source 1 Power axes
plot_flag=1;
k = waitforbuttonpress;
point1 = get(ha,'CurrentPoint');    % button down detected
%point1 = hobj.IntersectionPoint;
h.finalRect = rbbox;                   % return figure units
point2 = get(ha,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions
x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
hold on
axis manual
if (h.ax_power(1).XLim(1)>x(1) || x(1)>h.ax_power(1).XLim(2)) || (h.ax_power(1).XLim(1)>x(2) || x(2)>h.ax_power(1).XLim(2)) || ...  % checking if box is outside of graph's x-limits
        (h.ax_power(1).YLim(1)>y(1) || y(1)>h.ax_power(1).YLim(2)) || (h.ax_power(1).YLim(1)>y(3) || y(3)>h.ax_power(1).YLim(2)) || ...        % checking if box is outside of graph's y-limits
        (length(unique(x))==1 || length(unique(y))==1) % bounding box not drawn correctly or too quickly
    
    w=warndlg(sprintf('\nSelected ROI is outside of "Source 1 Power" graph\n\nPlease select ROI within "Source 1 Power" graph'));
    w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
    return;
else % create ROI bounding rectangular box
    % creating new TFR ROI data for all 3 sources
    for v=1:3
        try
            tfr_idx=length(h.tfr_ROI(v).h)+1;    % next roi box to plot
        catch
            tfr_idx=1;
        end
        
        % h.tfr_ROI(v).h(tfr_idx)=plot(x,y,'color',h.src_clr((v),:),'linewidth',2);
        h.tfr_ROI(v).h(tfr_idx)=plot(h.ax_power(v),x,y,'--','color',[1 1 1]*.7,'linewidth',1); %h.tfr_ROI(v).h(tfr_idx).Visible='off';
        h.tfr_ROI(v).h(tfr_idx).UserData.source_num=v;
        h.tfr_ROI(v).h(tfr_idx).UserData.tfr_idx=tfr_idx;
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.xpos=x;
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.ypos=y;
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_amp=[0 0 1 1 0 0];
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_amp=[1 1 0 0 1 1];
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_plv=[0 0 1 1 0 0];
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_plv=[1 1 0 0 1 1];
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_pli=[0 0 1 1 0 0];
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_pli=[1 1 0 0 1 1];
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi= [h.cfg.study.lat_sim(1) x(1) x(1)+range(x)/3 x(2)-range(x)/3 x(2) h.cfg.study.lat_sim(end)]; %[h.cfg.study.plot_time_int(1) x(1) x(1)+range(x)/3 x(2)-range(x)/3 x(2) h.cfg.study.plot_time_int(2)];
        h.tfr_ROI(v).h(tfr_idx).ButtonDownFcn=@selected_ROI;
        
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_evk_perc = str2num(h.edit_sig_evoked_perc(v).String);
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_evk_perc = str2num(h.edit_prepost_evoked_perc(v).String);
        
        
        % PLV info
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_phase_start = str2num(h.edit_sig_phase_start(v).String);
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_phase_start = str2num(h.edit_prepost_phase_start(v).String);
        
        freqs=round(y(1)):round(y(3));
        
        for f=1:length(freqs)
            h.tfr_ROI(v).h(tfr_idx).UserData.waves(f,:)=sin(2*pi*h.cfg.study.lat_sim*freqs(f)+(rand(1)*2*pi));
            h.tfr_ROI(v).h(tfr_idx).UserData.freq_idx(f)=freqs(f);
        end
        h.tfr_ROI(v).h(tfr_idx).UserData.deleted=0;
    end
    
    h.menu_tfr_roi_idx.String = {1:length(h.tfr_ROI(1).h)};
    h.menu_tfr_roi_idx.Value=length(h.tfr_ROI(1).h);
    
    h.menu_tfr_roi_idx_PLV.String = {1:length(h.tfr_ROI(1).h)}; %kl
    h.menu_tfr_roi_idx_PLV.Value=length(h.tfr_ROI(1).h); %kl
    
    h.menu_tfr_roi_idx_PLI.String = {1:length(h.tfr_ROI(1).h)}; %kl
    h.menu_tfr_roi_idx_PLI.Value=length(h.tfr_ROI(1).h); %kl
    h.current_tfr_idx = tfr_idx;
    
    % for v=1:3
    btn.Button=1; selected_ROI(h.tfr_ROI(v).h(tfr_idx),btn); % automatically plotting window function
    % end
    h.current_source_num=1; update_tfr_roi_cfg([],[]);
end
update_source_cfg(); update_PAC();
function selected_ROI(src,hobj)
global h

if hobj.Button==1   % left-click
    h.current_tfr_roi=src;
    h.current_tfr_idx=src.UserData.tfr_idx;
    h.current_tfr_roi_idx = src.UserData.tfr_idx;
    h.menu_tfr_roi_idx.Value = h.current_tfr_idx;
    %     h.menu_tfr_roi_idx_PLV.Value = h.current_tfr_idx; %kl
    %     h.menu_tfr_roi_idx_PLI.Value = h.current_tfr_idx; %kl
    h.current_source_num=src.UserData.source_num;
    plot_source_waves();
    plot_source_tfr(); %(h.current_tfr_ROI);
    update_tfr_roi_cfg([],[]);
    set_current_ROI(src,[],0);
    
elseif hobj.Button==3   % right-click
    h.current_tfr_roi=src;
    h.current_tfr_idx=src.UserData.tfr_idx;
    h.current_tfr_roi_idx = src.UserData.tfr_idx;
    h.menu_tfr_roi_idx.Value = h.current_tfr_idx;
    %     h.menu_tfr_roi_idx_PLV.Value = h.current_tfr_idx; %kl
    %     h.menu_tfr_roi_idx_PLI.Value = h.current_tfr_idx; %kl
    h.current_source_num=src.UserData.source_num;
    c = uicontextmenu;
    src.UIContextMenu = c;
    m1 = uimenu(c,'Label','Delete ROI Box?','Callback',{@delete_tfr_ROI,src});
    %     waitfor(m1);
    set_current_ROI(src,[],0);
end
function delete_tfr_ROI(~,~,~)
global h
cp_idx = setxor(1:length(h.tfr_ROI(1).h),h.current_tfr_idx);    % keeping TFR ROIs
new_idx = 1:length(cp_idx);
tfr_ROI = h.tfr_ROI;

for v=1:3
    tfr_ROI(v).h=[];
    x(new_idx)=h.tfr_ROI(v).h(cp_idx);
    tfr_ROI(v).h = handle(x);
    for t=1:length(tfr_ROI(v).h); tfr_ROI(v).h(t).UserData.source_num=v; tfr_ROI(v).h(t).UserData.tfr_idx=t; end
    delete(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_win)
    delete(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_win)
    delete(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.sig_waves)
    delete(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_win)
    delete(h.tfr_ROI(v).h(h.current_tfr_idx));
end

h.tfr_ROI=tfr_ROI;
h.menu_tfr_roi_idx.String = {1:length(h.tfr_ROI(1).h)};
h.menu_tfr_roi_idx_PLV.String = {1:length(h.tfr_ROI(1).h)};
h.menu_tfr_roi_idx_PLI.String = {1:length(h.tfr_ROI(1).h)};

h.menu_tfr_roi_idx.Value=length(h.tfr_ROI(1).h);
h.menu_tfr_roi_idx_PLV.Value = h.menu_tfr_roi_idx.Value;
h.menu_tfr_roi_idx_PLI.Value = h.menu_tfr_roi_idx.Value;

h.current_tfr_idx = h.menu_tfr_roi_idx.Value;
update_tfr_roi_cfg([],[]);
if h.current_tfr_idx>0
    plot_source_tfr(); plot_source_waves(); %(ev);
end

% clear source info then reload it with remaining tfr_roi info
h.cfg.source = struct('sig_freqs',[],'sig_amp_perc',[],'prepost_amp_perc',[],'sig_amp_perc_std',[],'prepost_amp_perc_std',[],'sig_evoked_perc',[],...
    'prepost_evoked_perc',[],'sig_durs',[],'sig_start',[],'sig_win_type',[],'sig_win_rise_time',[],'sig_PLV_targets',[],...
    'prepost_PLV_targets',[],'sig_PLI_targets',[],'prepost_PLI_targets',[],'sig_phase_lag',[],'prepost_phase_lag',[],'vx_ori',[],...
    'vx_idx',[],'vx_amp',[],'vx_locs',[],'phase_amp_contrasts',[],'sig_phase_amp_freq_idx',[],'prepost_phase_amp_freq_idx',[],...
    'sig_phase_amp_depth_perc',[],'prepost_phase_amp_depth_perc',[],'sig_phase_amp_depth_perc_range',[],'prepost_phase_amp_depth_perc_range',[]);

%% Update sources and PAC
update_source_cfg(); update_PAC();


%% update fcns
function update_graphs()
global h
% Set X limits for TFR and sig waves
if isfield(h,'ax_power') % can change Xlim once axes exists
    for a=1:3
        h.ax_power(a).XLim          = str2num(h.edit_plot_time_int.String);
        h.ax_sig_waves(a).XLim      = str2num(h.edit_plot_time_int.String);
        h.ax_prepost_waves(a).XLim  = str2num(h.edit_plot_time_int.String);
        h.ax_PLV(a).XLim            = str2num(h.edit_plot_time_int.String);
        h.ax_sig_plv(a).XLim        = str2num(h.edit_plot_time_int.String);
        h.ax_prepost_plv(a).XLim    = str2num(h.edit_plot_time_int.String);
        h.ax_PLI(a).XLim            = str2num(h.edit_plot_time_int.String);
        h.ax_sig_plv(a).XLim        = str2num(h.edit_plot_time_int.String);
        h.ax_prepost_pli(a).XLim    = str2num(h.edit_plot_time_int.String);
        
        h.ax_power(a).Colormap    = jet(255);
        h.ax_PLV(a).Colormap      = jet(255);
        h.ax_PLI(a).Colormap      = jet(255);
    end
    
    % Set Y limits for TFR and sig waves
    for a=1:3
        h.ax_power(a).YLim  = str2num(h.edit_plot_freq_int.String); h.ax_power(a).CLim = [-100 100];
        h.ax_PLV(a).YLim    = str2num(h.edit_plot_freq_int.String); h.ax_PLV(a).CLim = [-1 1];
        h.ax_PLI(a).YLim    = str2num(h.edit_plot_freq_int.String); h.ax_PLI(a).CLim = [-1 1];
        h.ax_sig_plv(a).YLim        = [0 1];
        h.ax_prepost_plv(a).YLim    = [0 1];
        h.ax_sig_pli(a).YLim        = [-0.5 0.5];
        h.ax_prepost_pli(a).YLim    = [-0.5 0.5];
        
        h.tfr_data.tfr_zero_line(a).YData = [h.cfg.study.plot_freq_int];
        h.tfr_data.plv_zero_line(a).YData = [h.cfg.study.plot_freq_int];
        h.tfr_data.pli_zero_line(a).YData = [h.cfg.study.plot_freq_int];
    end
end
function update_tfr_roi_cfg(varargin)
global h
if ~isempty(h.menu_tfr_roi_idx.String{1}) && h.current_tfr_idx<=length(h.tfr_ROI(h.current_source_num).h)
    if  h.tfr_ROI(h.current_source_num).h(h.current_tfr_idx).UserData.deleted==0
        update_source_edit_boxes();
        h.menu_tfr_roi_idx.Visible            = 'on';     % Set ERP Panel: no TFR ROIs exist
        h.menu_tfr_roi_idx_txt.Visible        = 'on';
        h.menu_source_idx_txt.Visible         = 'on';
        h.menu_source_idx.Visible             = 'on';
        h.edit_tfr_roi_xlim.Visible           = 'on';
        h.edit_tfr_roi_xlim_txt.Visible       = 'on';
        h.edit_tfr_roi_ylim.Visible           = 'on';
        h.edit_tfr_roi_ylim_txt.Visible       = 'on';
        h.edit_tfr_roi_risetime_txt.Visible   = 'on';
        h.edit_tfr_roi_risetime.Visible       = 'on';
        
        h.menu_tfr_roi_idx_PLV.Visible            = 'on';     % Set PLV Panel: no TFR ROIs exist
        h.menu_tfr_roi_idx_txt_PLV.Visible        = 'on';
        h.menu_source_idx_txt_PLV.Visible         = 'on';
        h.menu_source_idx_PLV.Visible             = 'on';
        h.edit_tfr_roi_xlim_PLV.Visible           = 'on';
        h.edit_tfr_roi_xlim_txt_PLV.Visible       = 'on';
        h.edit_tfr_roi_ylim_PLV.Visible           = 'on';
        h.edit_tfr_roi_ylim_txt_PLV.Visible       = 'on';
        h.edit_tfr_roi_risetime_txt_PLV.Visible   = 'on';
        h.edit_tfr_roi_risetime_PLV.Visible       = 'on';
        
        h.menu_tfr_roi_idx_PLI.Visible            = 'on';     % Set PLI Panel: no TFR ROIs exist
        h.menu_tfr_roi_idx_txt_PLI.Visible        = 'on';
        h.menu_source_idx_txt_PLI.Visible         = 'on';
        h.menu_source_idx_PLI.Visible             = 'on';
        h.edit_tfr_roi_xlim_PLI.Visible           = 'on';
        h.edit_tfr_roi_xlim_txt_PLI.Visible       = 'on';
        h.edit_tfr_roi_ylim_PLI.Visible           = 'on';
        h.edit_tfr_roi_ylim_txt_PLI.Visible       = 'on';
        h.edit_tfr_roi_risetime_txt_PLI.Visible   = 'on';
        h.edit_tfr_roi_risetime_PLI.Visible       = 'on';
    end
else
    h.menu_tfr_roi_idx.Visible            = 'on';     %  Set ERP Panel: no TFR ROIs exist
    h.menu_tfr_roi_idx_txt.Visible        = 'on';
    h.menu_source_idx_txt.Visible         = 'off';
    h.menu_source_idx.Visible             = 'off';
    h.edit_tfr_roi_xlim.Visible           = 'off';
    h.edit_tfr_roi_xlim_txt.Visible       = 'off';
    h.edit_tfr_roi_ylim.Visible           = 'off';
    h.edit_tfr_roi_ylim_txt.Visible       = 'off';
    h.edit_tfr_roi_risetime_txt.Visible   = 'off';
    h.edit_tfr_roi_risetime.Visible       = 'off';
    
    h.menu_tfr_roi_idx_PLV.Visible            = 'off';     %  Set PLV Panel: no TFR ROIs exist
    h.menu_tfr_roi_idx_txt_PLV.Visible        = 'off';
    h.menu_source_idx_txt_PLV.Visible         = 'off';
    h.menu_source_idx_PLV.Visible             = 'off';
    h.edit_tfr_roi_xlim_PLV.Visible           = 'off';
    h.edit_tfr_roi_xlim_txt_PLV.Visible       = 'off';
    h.edit_tfr_roi_ylim_PLV.Visible           = 'off';
    h.edit_tfr_roi_ylim_txt_PLV.Visible       = 'off';
    h.edit_tfr_roi_risetime_txt_PLV.Visible   = 'off';
    h.edit_tfr_roi_risetime_PLV.Visible       = 'off';
    
    h.menu_tfr_roi_idx_PLI.Visible            = 'off';     % Set PLI Panel: no TFR ROIs exist
    h.menu_tfr_roi_idx_txt_PLI.Visible        = 'off';
    h.menu_source_idx_txt_PLI.Visible         = 'off';
    h.menu_source_idx_PLI.Visible             = 'off';
    h.edit_tfr_roi_xlim_PLI.Visible           = 'off';
    h.edit_tfr_roi_xlim_txt_PLI.Visible       = 'off';
    h.edit_tfr_roi_ylim_PLI.Visible           = 'off';
    h.edit_tfr_roi_ylim_txt_PLI.Visible       = 'off';
    h.edit_tfr_roi_risetime_txt_PLI.Visible   = 'off';
    h.edit_tfr_roi_risetime_PLI.Visible       = 'off';
    
end
h.menu_source_idx.Value = h.current_source_num;
function update_source_edit_boxes()
global h
try
    
    if isfield(h,'tfr_ROI')
        v=h.current_source_num;
        %     for v=1:3
        h.edit_tfr_roi_xlim.String              = sprintf('%.3f %.3f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.xpos(1:2)); %kl: tfr_roi to tfr_idk
        h.edit_tfr_roi_xlim.ForegroundColor     = h.src_clr(v,:);
        h.edit_tfr_roi_ylim.String              = sprintf('%.f %.f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.ypos([1 3]));
        h.edit_tfr_roi_ylim.ForegroundColor     = h.src_clr(v,:);
        h.edit_tfr_roi_risetime.String          = sprintf('%.3f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi(3)-h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi(2));
        h.edit_tfr_roi_risetime.ForegroundColor = h.src_clr(v,:);
        %         h.edit_sig_power_perc(v).String = sprintf('%.f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_amp(3)*100);
        %         h.edit_prepost_power_perc(v).String = sprintf('%.f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_amp(1)*100);
        
        h.edit_tfr_roi_xlim_PLV.String              = sprintf('%.3f %.3f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.xpos(1:2));
        h.edit_tfr_roi_xlim_PLV.ForegroundColor     = h.src_clr(v,:);
        h.edit_tfr_roi_ylim_PLV.String              = sprintf('%.f %.f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.ypos([1 3]));
        h.edit_tfr_roi_ylim_PLV.ForegroundColor     = h.src_clr(v,:);
        h.edit_tfr_roi_risetime_PLV.String          = sprintf('%.3f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi(3)-h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi(2));
        h.edit_tfr_roi_risetime_PLV.ForegroundColor = h.src_clr(v,:);
        %h.edit_sig_power_perc_PLV(v).String = sprintf('%.f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_amp(3)*100);
        %h.edit_prepost_power_perc_PLV(v).String = sprintf('%.f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_amp(1)*100);
        
        h.edit_tfr_roi_xlim_PLI.String              = sprintf('%.3f %.3f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.xpos(1:2));
        h.edit_tfr_roi_xlim_PLI.ForegroundColor     = h.src_clr(v,:);
        h.edit_tfr_roi_ylim_PLI.String              = sprintf('%.f %.f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.ypos([1 3]));
        h.edit_tfr_roi_ylim_PLI.ForegroundColor     = h.src_clr(v,:);
        h.edit_tfr_roi_risetime_PLI.String          = sprintf('%.3f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi(3)-h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi(2));
        h.edit_tfr_roi_risetime_PLI.ForegroundColor = h.src_clr(v,:);
        %h.edit_sig_power_perc_PLV(v).String = sprintf('%.f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_amp(3)*100);
        %h.edit_prepost_power_perc_PLV(v).String = sprintf('%.f',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_amp(1)*100);
        %     end
        
        % Power
        h.edit_sig_power_perc(v).String = num2str(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_amp(3)*100);
        h.edit_prepost_power_perc(v).String = num2str(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_amp(1)*100);
        % Evoked
        h.edit_sig_evoked_perc(v).String = num2str(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_evk_perc);
        h.edit_prepost_evoked_perc(v).String = num2str(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_evk_perc);
        % PLV
        h.edit_sig_phase_locking(v).String = num2str(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_plv(3));
        h.edit_prepost_phase_locking(v).String = num2str(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_plv(1));
        h.edit_sig_phase_start(v).String = num2str(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_phase_start);
        h.edit_prepost_phase_start(v).String = num2str(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_phase_start);
        %PLI
        h.edit_sig_phase_lag(v).String = num2str(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_pli(3));
        h.edit_prepost_phase_lag(v).String = num2str(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_pli(1));
        
    end
catch
end

function update_plv_pli
global h
%% Update cfg for PLV tab
for v = 1:3
    %% Setting PLV waveform amplitudes
    h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_plv(3:4) = [str2num(h.edit_sig_phase_locking(v).String) str2num(h.edit_sig_phase_locking(v).String)];
    h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_plv([1 2 5 6]) = repmat(str2num(h.edit_prepost_phase_locking(v).String),1,4);
    h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_pli(3:4) = [str2num(h.edit_sig_phase_lag(v).String) str2num(h.edit_sig_phase_lag(v).String)];
    h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_pli([1 2 5 6]) = repmat(str2num(h.edit_prepost_phase_lag(v).String),1,4);
    
    % S1 vs. S2
    ft=fit(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi',[0 0 1 1 0 0]','linearinterp'); xwin=ft(h.cfg.study.lat_sim);     % windowing freq waves
    sig_plv(:,v) = xwin;
    h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi_plv = h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi;
    h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi_plv([2 3 4 5]) = [0 0 0 0];
end
% multiply sig_plv to find where intervals overlap --> if 0=no overlap >0=overlap and PLV/PLI can exist  1=flat top
% S1 vs S2, S1 vs S3, and S2 vs S3
cx = [1 2; 1 3; 2 3]; % PLV contrasts
for v=1:3
    x12 = sig_plv(:,cx(v,1)).*sig_plv(:,cx(v,2)); xs = find(x12>0);
    if ~isempty(xs)    % overlap exists
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi_plv([2 5]) = h.cfg.study.lat_sim(xs([1 end]));
        max_val = max(x12); xs = find(x12==max_val);
        if length(xs)==1
            xs = [xs xs+1];
            h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi_plv(5) = h.cfg.study.lat_sim(xs(end)+1);
        end
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi_plv([3 4]) = h.cfg.study.lat_sim(xs([1 end]));
        % PLV
        ft=fit(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi_plv',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_plv','linearinterp'); xwin=ft(h.cfg.study.lat_sim);     % windowing freq waves
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_plv_final = xwin; % x12;
        ft=fit(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi_plv',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_plv','linearinterp'); xwin=ft(h.cfg.study.lat_sim);     % windowing freq waves
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_plv_final = xwin; % x12;
        % PLI
        ft=fit(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi_plv',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_pli','linearinterp'); xwin=ft(h.cfg.study.lat_sim);     % windowing freq waves
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_pli_final = xwin; % x12;
        ft=fit(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi_plv',h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_pli','linearinterp'); xwin=ft(h.cfg.study.lat_sim);     % windowing freq waves
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_pli_final = xwin; % x12;
    else
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.x_roi_plv = [0 0 0 0 0 0];
        %         h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_plv = [0 0 0 0 0 0];
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_plv_final  = zeros(size(h.cfg.study.lat_sim));
        %         h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_pli = [0 0 0 0 0 0];
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.sig_pli_final  = zeros(size(h.cfg.study.lat_sim));
        
        %         h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_plv = [0 0 0 0 0 0];
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_plv_final  = zeros(size(h.cfg.study.lat_sim));
        %         h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_pli = [0 0 0 0 0 0];
        h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_pli_final  = zeros(size(h.cfg.study.lat_sim));
    end
    %     figure(v); clf; plot(h.tfr_ROI(v).h(h.current_tfr_idx).UserData.roi.prepost_plv_final);
end

%% plotting Power
function plot_source_tfr()
global h

for v=1:3
    tfr_data=zeros(length(h.freq_vals),length(h.cfg.study.lat_sim));
    tfr_data1=zeros(length(h.freq_vals),length(h.cfg.study.lat_sim),length(h.tfr_ROI(v).h),2);
    
    for r=1:length(h.tfr_ROI(v).h)
        if isvalid(h.tfr_ROI(v).h(r))
            for f=1:length(h.tfr_ROI(v).h(r).UserData.freq_idx)
                fidx(f)=find(h.tfr_data.s1.YData==h.tfr_ROI(v).h(r).UserData.freq_idx(f));
                tfr_data1(fidx(f),:,r,1)=h.tfr_ROI(v).h(r).UserData.roi.sig_win_final;     % sig_win
                tfr_data1(fidx(f),:,r,2)=h.tfr_ROI(v).h(r).UserData.roi.prepost_win_final;     % prepost_win
            end
            
        end
    end
    
    % summating sig + prepost windows to get envelope
    xwin=squeeze(sum(sum(tfr_data1,3),4))';
    xwin=(xwin-xwin(1,:));
    %     rwin=bsxfun(@rdivide,xwin,xwin(1,:));
    %     tfr_data=(20*log10(rwin))'; tfr_data(isnan(tfr_data))=0;
    tfr_data=100*(xwin)'; tfr_data(isnan(tfr_data))=0;
    
    
    if v==1
        %         h.tfr_data.s1 = surf(h.ax_power(1),h.cfg.study.lat_sim',h.freq_vals',tfr_data); view(h.ax_power(1),0,90); shading(h.ax_power(1),'interp'); caxis(h.ax_power(1),[-100 100])
        h.tfr_data.s1.XData = h.cfg.study.lat_sim;
        h.tfr_data.s1.YData = h.freq_vals;
        h.tfr_data.s1.CData = tfr_data;  h.tfr_data.s1.ZData=zeros(size(tfr_data));
    elseif v==2
        h.tfr_data.s2.XData = h.cfg.study.lat_sim;
        h.tfr_data.s2.YData = h.freq_vals;
        h.tfr_data.s2.CData = tfr_data;  h.tfr_data.s2.ZData=zeros(size(tfr_data));
    elseif v==3
        h.tfr_data.s3.XData = h.cfg.study.lat_sim;
        h.tfr_data.s3.YData = h.freq_vals;
        h.tfr_data.s3.CData = tfr_data;  h.tfr_data.s3.ZData=zeros(size(tfr_data));
    end
    h.ax_power(v).CLim = (str2num(h.edit_plot_caxis.String));
end
plot_plv_tfr();
function plot_source_waves()
global h
tfr_idx=h.current_tfr_idx;
for v=1:3
    
    axes(h.ax_sig_waves(v)); cla; hold on;
    ft=fit(h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi',h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_amp','linearinterp'); xwin=ft(h.cfg.study.lat_sim);     % windowing freq waves
    % ft=fit(h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi',h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_amp','gauss1'); xwin=ft(h.cfg.study.lat_sim);     % windowing freq waves
    % xwin=xwin/max(abs(xwin));
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_win_final=xwin;
    h.tfr_ROI(v).h(tfr_idx).UserData.sig_waves=plot(h.cfg.study.lat_sim,bsxfun(@times,h.tfr_ROI(v).h(tfr_idx).UserData.waves',xwin),'color',h.src_clr2(v,:),'linewidth',1);
    
    h.tfr_ROI(v).h(tfr_idx).UserData.sig_waves_evk=plot(h.cfg.study.lat_sim,bsxfun(@times,nanmean(h.tfr_ROI(v).h(tfr_idx).UserData.waves,1)',xwin*(str2num(h.edit_sig_evoked_perc(v).String)/100) ),'color',h.src_clr(v,:),'linewidth',2);
    
    
    
    %      h.tfr_ROI(v).h(tfr_idx).UserData.sig_waves_evk=plot(h.cfg.study.lat_sim,nanmean(bsxfun(@times,h.tfr_ROI(v).h(tfr_idx).UserData.waves',xwin),2),'color',h.evk_clr(v,:),'linewidth',2);
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_win=plot(h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi,h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_amp,'color',h.src_clr2(v,:)*0.2,'linewidth',2,'MarkerFaceColor',h.src_clr2(v,:),'MarkerEdgeColor','none','marker','o','MarkerSize',5,'hittest','on','buttondownfcn',@clickmarker); axis([h.cfg.study.plot_time_int -1 1]);
    %     h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_win.UserData=src.UserData;
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_win.UserData.source_num = v; h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_win.UserData.tfr_idx = tfr_idx;
    
    if v==1; legend([h.tfr_ROI(1).h(tfr_idx).UserData.sig_waves_evk,h.tfr_ROI(v).h(tfr_idx).UserData.sig_waves(1),h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_win],{'evoked' 'freqs' 'window'},'Location','northwest'); end
    
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_win.UserData.type='signal'; % used to designate which plots when passed to clickmarker
    h.ax_sig_waves(v).YLim = (str2num(h.edit_plot_caxis.String))/100;
    
    axes(h.ax_prepost_waves(v)); cla; hold on;
    ft=fit(h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi',h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_amp','linearinterp'); xwin=ft(h.cfg.study.lat_sim);     % windowing freq waves
    % ft=fit(h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi',h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_amp','gauss1'); xwin=ft(h.cfg.study.lat_sim);     % windowing freq waves
    % xwin=1-xwin;
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_win_final=xwin;
    h.tfr_ROI(v).h(tfr_idx).UserData.prepost_waves=plot(h.cfg.study.lat_sim,bsxfun(@times,h.tfr_ROI(v).h(tfr_idx).UserData.waves',xwin),'color',h.src_clr2(v,:),'linewidth',1);
    h.tfr_ROI(v).h(tfr_idx).UserData.prepost_waves_evk=plot(h.cfg.study.lat_sim,bsxfun(@times,nanmean(h.tfr_ROI(v).h(tfr_idx).UserData.waves,1)',xwin*(str2num(h.edit_prepost_evoked_perc(v).String)/100) ),'color',h.src_clr(v,:),'linewidth',2);
    
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_win=plot(h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi,h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_amp,'color',h.src_clr(v,:)*0.2,'linewidth',2,'MarkerFaceColor','none','MarkerEdgeColor',h.src_clr2(v,:),'marker','o','MarkerSize',5,'hittest','on','buttondownfcn',@clickmarker); axis([h.cfg.study.plot_time_int -1 1]);
    
    %     h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_win.UserData=src.UserData;
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_win.UserData.source_num = v; h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_win.UserData.tfr_idx = tfr_idx;
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_win.UserData.type='prepost';
    h.ax_prepost_waves(v).YLim = (str2num(h.edit_plot_caxis.String))/100;
    
end
h.ax_sig_waves(1).YLabel.String = 'Power (%/100)';
h.ax_prepost_waves(1).YLabel.String = 'Power (%/100)';
h.ax_prepost_waves(1).XLabel.String = 'Time (sec)';

%% plotting PLV
function plot_plv_tfr()
global h
% h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('\nPlotting PLV/PLI TFR ...\n'); drawnow;
for v=1:3
    vx = h.plv_contrasts(v,:);
    tfr_data=zeros(length(h.freq_vals),length(h.cfg.study.lat_sim));
    tfr_data1=zeros(length(h.freq_vals),length(h.cfg.study.lat_sim),length(h.tfr_ROI(v).h),2);
    update_plv_pli();
    
    for r=1:length(h.tfr_ROI(v).h)
        if isvalid(h.tfr_ROI(v).h(r))
            
            for f=1:length(h.tfr_ROI(v).h(r).UserData.freq_idx)
                
                fidx(f)=find(h.tfr_data.s1.YData==h.tfr_ROI(v).h(r).UserData.freq_idx(f));
                if sum(sum(abs([h.tfr_ROI(vx(1)).h(r).UserData.roi.sig_win_final h.tfr_ROI(vx(2)).h(r).UserData.roi.sig_win_final]))>0)==2 || ...
                        sum(sum(abs([h.tfr_ROI(vx(1)).h(r).UserData.roi.prepost_win_final h.tfr_ROI(vx(2)).h(r).UserData.roi.prepost_win_final]))>0)==2
                    tfr_data1(fidx(f),:,r,1)=h.tfr_ROI(v).h(r).UserData.roi.sig_plv_final;     % sig_win
                    tfr_data1(fidx(f),:,r,2)=h.tfr_ROI(v).h(r).UserData.roi.prepost_plv_final;     % prepost_win
                else
                    tfr_data1(fidx(f),:,r,1)=zeros(size(h.tfr_ROI(v).h(r).UserData.roi.sig_plv_final));     % sig_win
                    tfr_data1(fidx(f),:,r,2)=zeros(size(h.tfr_ROI(v).h(r).UserData.roi.prepost_plv_final));     % prepost_win
                end
                
            end
            
        end
    end
    
    % summating sig + prepost windows to get envelope
    xwin=squeeze(sum(sum(tfr_data1,3),4))';
    xwin=(xwin-xwin(1,:));
    %     rwin=bsxfun(@rdivide,xwin,xwin(1,:));
    %     tfr_data=(20*log10(rwin))'; tfr_data(isnan(tfr_data))=0;
    tfr_data=(xwin)'; tfr_data(isnan(tfr_data))=0;
    
    
    if v==1
        %         h.plv_data.s1 = surf(h.ax_power(1),h.cfg.study.lat_sim',h.freq_vals',tfr_data); view(h.ax_power(1),0,90); shading(h.ax_power(1),'interp'); caxis(h.ax_power(1),[-100 100])
        h.plv_data.s1.XData = h.cfg.study.lat_sim;
        h.plv_data.s1.YData = h.freq_vals;
        h.plv_data.s1.CData = tfr_data;  h.plv_data.s1.ZData=zeros(size(tfr_data));
    elseif v==2
        h.plv_data.s2.XData = h.cfg.study.lat_sim;
        h.plv_data.s2.YData = h.freq_vals;
        h.plv_data.s2.CData = tfr_data;  h.plv_data.s2.ZData=zeros(size(tfr_data));
    elseif v==3
        h.plv_data.s3.XData = h.cfg.study.lat_sim;
        h.plv_data.s3.YData = h.freq_vals;
        h.plv_data.s3.CData = tfr_data;  h.plv_data.s3.ZData=zeros(size(tfr_data));
    end
    
end
plot_plv_waves(); plot_pli_tfr();
% h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
function plot_plv_waves()
global h
tfr_idx=h.current_tfr_idx;
for v=1:3
    vx = h.plv_contrasts(v,:);
    if sum(sum(abs([h.tfr_ROI(vx(1)).h(tfr_idx).UserData.roi.sig_win_final h.tfr_ROI(vx(2)).h(tfr_idx).UserData.roi.sig_win_final]))>0)==2 || ...
            sum(sum(abs([h.tfr_ROI(vx(1)).h(tfr_idx).UserData.roi.prepost_win_final h.tfr_ROI(vx(2)).h(tfr_idx).UserData.roi.prepost_win_final]))>0)==2
        sig_data = h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_plv_final;
        prepost_data = h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_plv_final;
    else
        sig_data = zeros(size(h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_plv_final));
        prepost_data = zeros(size(h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_plv_final));
    end
    
    axes(h.ax_sig_plv(v)); cla; hold on;
    plot(h.cfg.study.lat_sim, sig_data,'color',h.plv_clr2(v,:)*0,'linewidth',2);
    axis([h.cfg.study.plot_time_int 0 1]);
    
    axes(h.ax_prepost_plv(v)); cla; hold on;
    plot(h.cfg.study.lat_sim, prepost_data,'color',h.plv_clr2(v,:)*0,'linewidth',2);
    axis([h.cfg.study.plot_time_int 0 1]);
    
    
end
h.ax_sig_plv(1).YLabel.String = 'PLV';
h.ax_prepost_plv(1).YLabel.String = 'PLV';
h.ax_prepost_plv(1).XLabel.String = 'Time (sec)';


%% plotting PLI
function plot_pli_tfr()
global h

for v=1:3
    vx = h.plv_contrasts(v,:);
    tfr_data=zeros(length(h.freq_vals),length(h.cfg.study.lat_sim));
    tfr_data1=zeros(length(h.freq_vals),length(h.cfg.study.lat_sim),length(h.tfr_ROI(v).h),2);
    %     update_plv_pli();
    
    for r=1:length(h.tfr_ROI(v).h)
        if isvalid(h.tfr_ROI(v).h(r))
            
            
            for f=1:length(h.tfr_ROI(v).h(r).UserData.freq_idx)
                fidx(f)=find(h.tfr_data.s1.YData==h.tfr_ROI(v).h(r).UserData.freq_idx(f));
                if sum(sum(abs([h.tfr_ROI(vx(1)).h(r).UserData.roi.sig_win_final h.tfr_ROI(vx(2)).h(r).UserData.roi.sig_win_final]))>0)==2 || ...
                        sum(sum(abs([h.tfr_ROI(vx(1)).h(r).UserData.roi.prepost_win_final h.tfr_ROI(vx(2)).h(r).UserData.roi.prepost_win_final]))>0)==2
                    tfr_data1(fidx(f),:,r,1)=h.tfr_ROI(v).h(r).UserData.roi.sig_pli_final;     % sig_win
                    tfr_data1(fidx(f),:,r,2)=h.tfr_ROI(v).h(r).UserData.roi.prepost_pli_final;     % prepost_win
                else
                    tfr_data1(fidx(f),:,r,1)=zeros(size(h.tfr_ROI(v).h(r).UserData.roi.sig_pli_final));     % sig_win
                    tfr_data1(fidx(f),:,r,2)=zeros(size(h.tfr_ROI(v).h(r).UserData.roi.prepost_pli_final));     % prepost_win
                end
                
            end
            
        end
    end
    
    % summating sig + prepost windows to get envelope
    xwin=squeeze(sum(sum(tfr_data1,3),4))';
    xwin=(xwin-xwin(1,:));
    %     rwin=bsxfun(@rdivide,xwin,xwin(1,:));
    %     tfr_data=(20*log10(rwin))'; tfr_data(isnan(tfr_data))=0;
    tfr_data=(xwin)'; tfr_data(isnan(tfr_data))=0;
    
    
    if v==1
        %         h.pli_data.s1 = surf(h.ax_power(1),h.cfg.study.lat_sim',h.freq_vals',tfr_data); view(h.ax_power(1),0,90); shading(h.ax_power(1),'interp'); caxis(h.ax_power(1),[-100 100])
        h.pli_data.s1.XData = h.cfg.study.lat_sim;
        h.pli_data.s1.YData = h.freq_vals;
        h.pli_data.s1.CData = tfr_data;  h.pli_data.s1.ZData=zeros(size(tfr_data));
    elseif v==2
        h.pli_data.s2.XData = h.cfg.study.lat_sim;
        h.pli_data.s2.YData = h.freq_vals;
        h.pli_data.s2.CData = tfr_data;  h.pli_data.s2.ZData=zeros(size(tfr_data));
    elseif v==3
        h.pli_data.s3.XData = h.cfg.study.lat_sim;
        h.pli_data.s3.YData = h.freq_vals;
        h.pli_data.s3.CData = tfr_data;  h.pli_data.s3.ZData=zeros(size(tfr_data));
    end
    
end
plot_pli_waves();
function plot_pli_waves()
global h
tfr_idx=h.current_tfr_idx;
for v=1:3
    vx = h.plv_contrasts(v,:);
    if sum(sum(abs([h.tfr_ROI(vx(1)).h(tfr_idx).UserData.roi.sig_win_final h.tfr_ROI(vx(2)).h(tfr_idx).UserData.roi.sig_win_final]))>0)==2 || ...
            sum(sum(abs([h.tfr_ROI(vx(1)).h(tfr_idx).UserData.roi.prepost_win_final h.tfr_ROI(vx(2)).h(tfr_idx).UserData.roi.prepost_win_final]))>0)==2
        sig_data = h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_pli_final;
        prepost_data = h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_pli_final;
    else
        sig_data = zeros(size(h.tfr_ROI(v).h(tfr_idx).UserData.roi.sig_pli_final));
        prepost_data = zeros(size(h.tfr_ROI(v).h(tfr_idx).UserData.roi.prepost_pli_final));
    end
    
    axes(h.ax_sig_pli(v)); cla; hold on;
    plot(h.cfg.study.lat_sim,  sig_data,'color',h.plv_clr2(v,:)*0,'linewidth',2);
    axis([h.cfg.study.plot_time_int -1 1]);
    
    axes(h.ax_prepost_pli(v)); cla; hold on;
    plot(h.cfg.study.lat_sim, prepost_data,'color',h.plv_clr2(v,:)*0,'linewidth',2);
    axis([h.cfg.study.plot_time_int -1 1]);
    
end
h.ax_sig_pli(1).YLabel.String = 'dPLI';
h.ax_prepost_pli(1).YLabel.String = 'dPLI';
h.ax_prepost_pli(1).XLabel.String = 'Time (sec)';



%% set TFR ROIs
function set_tfr_roi_cfg(src,~)
global h
h.current_tfr_roi_idx = h.menu_tfr_roi_idx.Value;
src_str = src.String;
try
    if ~isempty(src.Tag)
        sig_type =  src.Tag(1:5);
        box_type =  src.Tag(7:end);
        switch sig_type
            case 'Sig 1'; h.current_source_num = 1; h.menu_source_idx.Value = 1;     update_source_edit_boxes;  src.String = src_str;
            case 'Sig 2'; h.current_source_num = 2; h.menu_source_idx.Value = 2;     update_source_edit_boxes;  src.String = src_str;
            case 'Sig 3'; h.current_source_num = 3; h.menu_source_idx.Value = 3;     update_source_edit_boxes;  src.String = src_str;
            case 'Pre 1'; h.current_source_num = 1; h.menu_source_idx.Value = 1;     update_source_edit_boxes;  src.String = src_str;
            case 'Pre 2'; h.current_source_num = 2; h.menu_source_idx.Value = 2;     update_source_edit_boxes;  src.String = src_str;
            case 'Pre 3'; h.current_source_num = 3; h.menu_source_idx.Value = 3;     update_source_edit_boxes;  src.String = src_str;
        end
        
        %     switch box_type % checking type of edit box being called
        %         case 'Latency Range'
        %         case 'Frequency Range'
        %         case 'Power'
        %             update_source_edit_boxes;
        %         case 'PLV'
        %             update_source_edit_boxes;
        %         case 'PLI'
        %             update_source_edit_boxes;
        %     end
    end
catch
    fprintf('Error in "set_tfr_roi_cfg"\n');
end

if ~isempty(h.menu_tfr_roi_idx.String{1})   % TFR ROIs exist
    v=h.current_source_num;
    tfr_idx = h.menu_tfr_roi_idx.Value;
    
    % making sure windowing is within data limits
    %     v= h.current_source_num; tfr_idx = h.current_tfr_idx;
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi([2 5]) = str2num(h.edit_tfr_roi_xlim.String);
    
    if diff(h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi(1:2))==0 || h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi(2)<h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi(1)
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi(2)=h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi(1)+1/h.cfg.study.srate;
    end
    if diff(h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi(5:6))==0 || h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi(5)>h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi(6)
        h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi(5)=h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi(6)-1/h.cfg.study.srate;
    end
    %    sprintf('%.3f %.3f',h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi([2 5]))
    h.edit_tfr_roi_xlim.String = sprintf('%.3f %.3f',h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi([2 5]));
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.xpos(1:2) = h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi([2 5]);
    
    % setting new ROI limits
    xlim=str2num(h.edit_tfr_roi_xlim.String);
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.xpos = [xlim(1) xlim(2) xlim(2) xlim(1) xlim(1)];
    h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi([2 5])= xlim;
    rise_lat = str2num(h.edit_tfr_roi_risetime.String);
    xrise = [xlim(1)+rise_lat xlim(2)-rise_lat];
    if xrise(1)>=xrise(2); rise_lat=(range(h.tfr_ROI(v).h(tfr_idx).UserData.roi.x_roi([2 5]))/2)-(1/h.cfg.study.srate);
        xrise = [xlim(1)+rise_lat xlim(2)-rise_lat]; h.edit_tfr_roi_risetime.String = sprintf('%.3f',rise_lat); end
    h.tfr_ROI(h.current_source_num).h(tfr_idx).UserData.roi.x_roi([3 4])= xrise;
    
    for vx = 1:3
        ylim=str2num(h.edit_tfr_roi_ylim.String);
        if ylim(1)<=0; ylim(1)=1; elseif ylim(2)<=0 ylim(2)=1; end; if ylim(2)<ylim(1); ylim(2)=ylim(1); elseif ylim(1)>ylim(2); ylim(1)=ylim(2); end
        h.edit_tfr_roi_ylim.String = sprintf('%.f %.f',ylim);
        h.tfr_ROI(vx).h(tfr_idx).UserData.roi.ypos = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)];
        freqs=round(ylim(1)):round(ylim(2));
        h.tfr_ROI(vx).h(tfr_idx).UserData.freq_idx=[]; h.tfr_ROI(vx).h(tfr_idx).UserData.waves=[]; % clearing
        for f=1:length(freqs)
            h.tfr_ROI(vx).h(tfr_idx).UserData.waves(f,:)=sin(2*pi*h.cfg.study.lat_sim*freqs(f)+(rand(1)*2*pi));
            h.tfr_ROI(vx).h(tfr_idx).UserData.freq_idx(f)=freqs(f);
        end
        h.tfr_ROI(vx).h(tfr_idx).UserData.deleted=0;
        h.tfr_ROI(vx).h(tfr_idx).XData = h.tfr_ROI(vx).h(tfr_idx).UserData.roi.xpos;
        h.tfr_ROI(vx).h(tfr_idx).YData = h.tfr_ROI(vx).h(tfr_idx).UserData.roi.ypos;
        
    end
    
    %% setting waveform amplitudes
    check_edit_box_limits(src);
    
    h.tfr_ROI(h.current_source_num).h(h.current_tfr_roi_idx).UserData.roi.sig_amp(3:4) = [str2num(h.edit_sig_power_perc(h.current_source_num).String)/100 str2num(h.edit_sig_power_perc(h.current_source_num).String)/100];
    h.tfr_ROI(h.current_source_num).h(h.current_tfr_roi_idx).UserData.roi.prepost_amp([1 2 5 6]) = repmat(str2num(h.edit_prepost_power_perc(h.current_source_num).String)/100,1,4);
    
    %     update_plv_pli();   % updating plv and pli parameters
    
    plot_source_waves(); plot_source_tfr();
end
update_tfr_roi_cfg([],[]);
function set_current_ROI(src,~,curr_tab)
global h
current_source_num = h.current_source_num;
if isprop(src,'Value')
    h.current_tfr_idx = src.Value;
elseif isfield(src,'UserData')
    h.current_tfr_idx = src.UserData.tfr_idx;
else
    h.current_tfr_idx = h.menu_tfr_roi_idx.Value;
end

% making all ROI menu have same tfr_idx value
h.menu_tfr_roi_idx.Value=h.current_tfr_idx; h.menu_tfr_roi_idx_PLV.Value=h.current_tfr_idx; h.menu_tfr_roi_idx_PLI.Value=h.current_tfr_idx;


for v=1:3   % looping through sources to update edit box info and plots
    h.current_source_num = v;
    
    if (curr_tab == 0)
        h.menu_source_idx.Value = h.current_source_num;
    elseif(curr_tab == 1)
        h.menu_source_idx_PLV.Value = h.current_source_num;
    elseif(curr_tab == 2)
        h.menu_source_idx_PLI.Value = h.current_source_num;
    end
    update_source_edit_boxes;
    
    %     set_current_source([],[],curr_tab)
end
h.current_source_num = current_source_num;  % resetting back to original
h.current_tfr_idx = h.menu_tfr_roi_idx.Value;
if (curr_tab == 0)
    h.menu_source_idx.Value = h.current_source_num;
elseif(curr_tab == 1)
    h.menu_source_idx_PLV.Value = h.current_source_num;
elseif(curr_tab == 2)
    h.menu_source_idx_PLI.Value = h.current_source_num;
end
set_current_source([],[],curr_tab);
function set_current_source(~,~,curr_tab)
global h

if (curr_tab == 0)
    h.current_source_num = h.menu_source_idx.Value;
    h.current_tfr_idx = h.menu_tfr_roi_idx.Value;
elseif(curr_tab == 1)
    h.current_source_num = h.menu_source_idx_PLV.Value;
    h.current_tfr_idx = h.menu_tfr_roi_idx_PLV.Value;
elseif(curr_tab == 2)
    h.current_source_num = h.menu_source_idx_PLI.Value;
    h.current_tfr_idx = h.menu_tfr_roi_idx_PLI.Value;
end

update_tfr_roi_cfg();
if ~isempty(h.menu_tfr_roi_idx.String{1})
    
    h.menu_tfr_roi_idx.Value            = h.current_tfr_idx;
    h.menu_tfr_roi_idx_PLV.Value        = h.current_tfr_idx;
    h.menu_tfr_roi_idx_PLI.Value        = h.current_tfr_idx;
    h.menu_source_idx.Value             =  h.current_source_num;
    h.menu_source_idx_PLV.Value         =  h.current_source_num;
    h.menu_source_idx_PLI.Value         =  h.current_source_num;
    plot_source_waves();  plot_source_tfr()
end
function check_edit_box_limits(src)
global h

switch src.Tag
    case 'Sig 1'; h.current_source_num = 1; h.menu_source_idx.Value = 1;
    case 'Sig 2'; h.current_source_num = 2; h.menu_source_idx.Value = 2;
    case 'Sig 3'; h.current_source_num = 3; h.menu_source_idx.Value = 3;
    case 'Pre 1'; h.current_source_num = 1; h.menu_source_idx.Value = 1;
    case 'Pre 2'; h.current_source_num = 2; h.menu_source_idx.Value = 2;
    case 'Pre 3'; h.current_source_num = 3; h.menu_source_idx.Value = 3;
end

%% Power
% Sig
if str2num(h.edit_sig_power_perc(h.current_source_num).String)<0; h.edit_sig_power_perc(h.current_source_num).String='0';
elseif str2num(h.edit_sig_power_perc(h.current_source_num).String)>100; h.edit_sig_power_perc(h.current_source_num).String='100';
end
% Prepost
if str2num(h.edit_prepost_power_perc(h.current_source_num).String)<0; h.edit_prepost_power_perc(h.current_source_num).String='0';
elseif str2num(h.edit_prepost_power_perc(h.current_source_num).String)>100; h.edit_prepost_power_perc(h.current_source_num).String='100';
end
%% Evoked
% Sig
if str2num(h.edit_sig_evoked_perc(h.current_source_num).String)<0; h.edit_sig_evoked_perc(h.current_source_num).String='0';
elseif str2num(h.edit_sig_evoked_perc(h.current_source_num).String)>100; h.edit_sig_evoked_perc(h.current_source_num).String='100';
end
% Prepost
if str2num(h.edit_prepost_evoked_perc(h.current_source_num).String)<0; h.edit_prepost_evoked_perc(h.current_source_num).String='0';
elseif str2num(h.edit_prepost_evoked_perc(h.current_source_num).String)>100; h.edit_prepost_evoked_perc(h.current_source_num).String='100';
end
%% PLV
% Sig
if str2num(h.edit_sig_phase_locking(h.current_source_num).String)<0; h.edit_sig_phase_locking(h.current_source_num).String='0';
elseif str2num(h.edit_sig_phase_locking(h.current_source_num).String)>1; h.edit_sig_phase_locking(h.current_source_num).String='1';
end
% Prepost
if str2num(h.edit_prepost_phase_locking(h.current_source_num).String)<0; h.edit_prepost_phase_locking(h.current_source_num).String='0';
elseif str2num(h.edit_prepost_phase_locking(h.current_source_num).String)>1; h.edit_prepost_phase_locking(h.current_source_num).String='1';
end
%% PLI
% Sig
% if str2num(h.edit_sig_phase_lag(h.current_source_num).String)<0; h.edit_sig_phase_lag(h.current_source_num).String='0';
if abs(str2num(h.edit_sig_phase_lag(h.current_source_num).String))>1; h.edit_sig_phase_lag(h.current_source_num).String='1';
end
% Prepost
% if str2num(h.edit_prepost_phase_lag(h.current_source_num).String)<0; h.edit_prepost_phase_lag(h.current_source_num).String='0';
if abs(str2num(h.edit_prepost_phase_lag(h.current_source_num).String))>1; h.edit_prepost_phase_lag(h.current_source_num).String='0';
end

% updated TFR data with changes made in edit boxes
h.tfr_ROI(h.current_source_num).h(h.current_tfr_idx).UserData.roi.sig_amp([3 4]) =  round(str2num(h.edit_sig_power_perc(h.current_source_num).String))/100;
h.tfr_ROI(h.current_source_num).h(h.current_tfr_idx).UserData.roi.prepost_amp([1 2 5 6]) =  round(str2num(h.edit_prepost_power_perc(h.current_source_num).String))/100;

h.tfr_ROI(h.current_source_num).h(h.current_tfr_idx).UserData.roi.sig_evk_perc =  str2num(h.edit_sig_evoked_perc(h.current_source_num).String);
h.tfr_ROI(h.current_source_num).h(h.current_tfr_idx).UserData.roi.prepost_evk_perc =  str2num(h.edit_prepost_evoked_perc(h.current_source_num).String);

h.tfr_ROI(h.current_source_num).h(h.current_tfr_idx).UserData.roi.sig_plv([3 4]) =  round(str2num(h.edit_sig_phase_locking(h.current_source_num).String));
h.tfr_ROI(h.current_source_num).h(h.current_tfr_idx).UserData.roi.prepost_plv([1 2 5 6]) =  round(str2num(h.edit_prepost_phase_locking(h.current_source_num).String));

h.tfr_ROI(h.current_source_num).h(h.current_tfr_idx).UserData.roi.sig_pli([3 4]) =  round(str2num(h.edit_sig_phase_lag(h.current_source_num).String));
h.tfr_ROI(h.current_source_num).h(h.current_tfr_idx).UserData.roi.prepost_pli([1 2 5 6]) =  round(str2num(h.edit_prepost_phase_lag(h.current_source_num).String));


%% Dragging Time Windowing Shape
function clickmarker(src,~)
global h
h.current_source_num=src.UserData.source_num;
h.current_tfr_roi=h.tfr_ROI(h.current_source_num).h; %h.current_source_num
set(ancestor(src,'figure'),'windowbuttonmotionfcn',{@dragmarker,src})
set(ancestor(src,'figure'),'windowbuttonupfcn',@stopdragging)
function dragmarker(~,~,src)
global h
% h.current_source_num=src.UserData(1).source_num;

%get current axes and coords
h1=gca;
% h1=h.ax_sig_waves(h.current_source_num);
coords=get(h1,'currentpoint');

%get all x and y data
x=h1.Children.XData;
y=h1.Children.YData;
%check which data point has the smallest distance to the dragged point
x_diff=abs(x-coords(1,1,1));
y_diff=abs(y-coords(1,2,1));
[value index]=min(x_diff+y_diff);
%create new x and y data and exchange coords for the dragged point
x_new=x; y_new=y;

if coords(1,2,1)<0; coords(1,2,1)=0;
elseif coords(1,2,1)>1; coords(1,2,1)=1;
end
min_diff=1e-9;
% update tfr_roi data
switch src.UserData.type
    case 'signal'
        if index==3
            xdiff=(coords(1,1,1)-x_new(2));
            if xdiff<=0
                xdiff=min_diff;
                x_new(3)=x_new(2)+min_diff;
            else
                x_new(3)=coords(1,1,1);
            end
            y_new([3 4])=coords(1,2,1);
            x_new(4)=x_new(5)-xdiff;
        elseif  index==4
            xdiff=(x_new(5)-coords(1,1,1));
            if xdiff<=0
                xdiff=min_diff;
                x_new(4)=x_new(5)-min_diff;
            else
                x_new(4)=coords(1,1,1);
            end
            y_new([3 4])=coords(1,2,1);
            x_new(3)=x_new(2)+xdiff;
        else
            %             % shifting pointd 1 & 2 y-coords together
            %             y_new([1 2 5 6])=coords(1,2,1);
        end
        %update plot
        set(src,'xdata',x_new,'ydata',y_new);
        
        h.tfr_ROI(h.current_source_num).h(src.UserData.tfr_idx).UserData.roi.sig_amp=y_new;
        h.tfr_ROI(h.current_source_num).h(src.UserData.tfr_idx).UserData.roi.x_roi=x_new;
        h.pt_txt.String = sprintf('Signal Rise Lat = %.3f     Amp = %.3f',h.tfr_ROI(h.current_source_num).h(src.UserData.tfr_idx).UserData.roi.x_roi(3),h.tfr_ROI(h.current_source_num).h(src.UserData.tfr_idx).UserData.roi.sig_amp(3));
        
    case 'prepost'
        if index==3
            xdiff=(coords(1,1,1)-x_new(2));
            if xdiff<=0
                xdiff=min_diff;
                x_new(3)=x_new(2)+min_diff;
            else
                x_new(3)=coords(1,1,1);
            end
            y_new([3 4])=0;
            x_new(4)=x_new(5)-xdiff;
        elseif  index==4
            xdiff=(x_new(5)-coords(1,1,1));
            if xdiff<=0
                xdiff=min_diff;
                x_new(4)=x_new(5)-min_diff;
            else
                x_new(4)=coords(1,1,1);
            end
            y_new([3 4])=0;
            x_new(3)=x_new(2)+xdiff;
        else
            % shifting pointd 1 & 2 y-coords together
            y_new([1 2 5 6])=coords(1,2,1);
        end
        %update plot
        set(src,'xdata',x_new,'ydata',y_new);
        
        h.tfr_ROI(h.current_source_num).h(src.UserData.tfr_idx).UserData.roi.prepost_amp=y_new;
        h.tfr_ROI(h.current_source_num).h(src.UserData.tfr_idx).UserData.roi.x_roi=x_new;
        
        h.pt_txt.String = sprintf('PrePost Fall Lat = %.3f     Amp = %.3f',h.tfr_ROI(h.current_source_num).h(src.UserData.tfr_idx).UserData.roi.x_roi(3),h.tfr_ROI(h.current_source_num).h(src.UserData.tfr_idx).UserData.roi.prepost_amp(2));
        
end
h.tfr_ROI(h.current_source_num).h(src.UserData.tfr_idx).UserData.roi.sig_win.XData=x_new;
h.tfr_ROI(h.current_source_num).h(src.UserData.tfr_idx).UserData.roi.prepost_win.XData=x_new;
function stopdragging(fig,~)
plot_source_waves(); plot_source_tfr();
update_tfr_roi_cfg
set(fig,'windowbuttonmotionfcn','')
set(fig,'windowbuttonupfcn','')


%% %%%%%%% "PAC" Tab Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_PAC(varargin)
global h



if h.load_study_flag == 0
    h.cfg.source.phase_amp_contrasts                        = [1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; % These are fixed --> sig_contrasts = cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
    h.cfg.source.sig_phase_amp_freq_idx                     = zeros(size(h.cfg.source.phase_amp_contrasts,1), size(h.cfg.source.sig_freqs,2)); %[0 0; 0 0; 0 0 ; 0 0 ; 0 0 ; 0 0 ]; %(sig_contrasts x Nfreqs);  % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
    h.cfg.source.prepost_phase_amp_freq_idx                 = zeros(size(h.cfg.source.phase_amp_contrasts,1), size(h.cfg.source.sig_freqs,2));   %(PLV_contrasts x Nfreqs);     % cross-frequency contrasts among 3signals % sig_contrasts=[1 2; 1 3; 2 1; 2 3; 3 1; 3 2]; Note: 1st index is the signal amplitude being modulated by the 2nd signal
    h.cfg.source.sig_phase_amp_depth_perc                   = zeros(size(h.cfg.source.phase_amp_contrasts,1), size(h.cfg.source.sig_freqs,2)); % amplitude-modulation depth as a percentage of signal's amplitude (sig_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
    h.cfg.source.prepost_phase_amp_depth_perc               = zeros(size(h.cfg.source.phase_amp_contrasts,1), size(h.cfg.source.sig_freqs,2)); % depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx)
    h.cfg.source.sig_phase_amp_depth_perc_range             = zeros(size(h.cfg.source.phase_amp_contrasts,1), size(h.cfg.source.sig_freqs,2)); % +/- range of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: sig_phase_amp_depth_perc +/- sig_phase_amp_depth_perc_range must be within [0 100]
    h.cfg.source.prepost_phase_amp_depth_perc_range         = zeros(size(h.cfg.source.phase_amp_contrasts,1), size(h.cfg.source.sig_freqs,2)); % +/- range devitaion of depth percentage of prepost's amplitude (prepost_amp_perc) modulated at phase of sig_freq(sig_phase_amp_freq_idx). NOTE: prepost_phase_amp_depth_perc +/- prepost_phase_amp_depth_perc_range must be within [0 100]
    
    %% Updating Carrier & Modulator Menu items based on frequencies selected for TFR ROI
    if sum(sum(sum(h.cfg.source.sig_freqs)))>0
        for a=1:6
            if size(h.cfg.source.sig_freqs,2)>1
                sig_freqs = num2str (squeeze(h.cfg.source.sig_freqs(1,:,:)));
            else
                sig_freqs = num2str (squeeze(h.cfg.source.sig_freqs(1,:,:))');
            end
            
            h.menu_PAC_sig_carrier_source_freq(a).String = sig_freqs;
            h.menu_PAC_sig_modulator_source_freq(a).String = sig_freqs;
            h.menu_PAC_prepost_carrier_source_freq(a).String = sig_freqs;
            h.menu_PAC_prepost_modulator_source_freq(a).String = sig_freqs;
        end
    else
        fprintf('No ROI frequencies selected\n');
    end
    
    %% updating PAC parameters
    for a = 1:6
        % Signal freq
        h.cfg.source.sig_phase_amp_freq_idx(a,h.menu_PAC_sig_carrier_source_freq(a).Value) = h.menu_PAC_sig_modulator_source_freq(a).Value; % a = PAC contrast, h.menu_PAC_sig_carrier_source_freq(a).Value=index to set the number for the modulator index that is set by the h.menu_PAC_sig_modulator_source_freq(a).Value
        % Signal Modulator Depth
        h.cfg.source.sig_phase_amp_depth_perc(a,h.menu_PAC_sig_carrier_source_freq(a).Value) = str2num(h.edit_PAC_modulator_sig_depth(a).String{:}); % index of carrier frequency relative to selected TFR ROI
        % Signal Modulator Depth Range
        h.cfg.source.sig_phase_amp_depth_perc_range(a,h.menu_PAC_sig_carrier_source_freq(a).Value) = str2num(h.edit_PAC_modulator_sig_depth_range(a).String{:}); % index of carrier frequency relative to selected TFR ROI
        
        % PrePost freq
        h.cfg.source.prepost_phase_amp_freq_idx(a,h.menu_PAC_prepost_carrier_source_freq(a).Value) = h.menu_PAC_prepost_modulator_source_freq(a).Value; % index of carrier frequency relative to selected TFR ROI
        % PrePost Modulator Depth
        h.cfg.source.prepost_phase_amp_depth_perc(a,h.menu_PAC_prepost_carrier_source_freq(a).Value) = str2num(h.edit_PAC_modulator_prepost_depth(a).String{:}); % index of carrier frequency relative to selected TFR ROI
        % PrePost Modulator Depth Range
        h.cfg.source.prepost_phase_amp_depth_perc_range(a,h.menu_PAC_prepost_carrier_source_freq(a).Value) = str2num(h.edit_PAC_modulator_prepost_depth_range(a).String{:}); % index of carrier frequency relative to selected TFR ROI
    end
    
    sig_min_depth = find( (h.cfg.source.sig_phase_amp_depth_perc - h.cfg.source.sig_phase_amp_depth_perc_range)<0, 1);
    sig_max_depth = find( (h.cfg.source.sig_phase_amp_depth_perc + h.cfg.source.sig_phase_amp_depth_perc_range)>100, 1);
    prepost_min_depth = find( (h.cfg.source.prepost_phase_amp_depth_perc - h.cfg.source.prepost_phase_amp_depth_perc_range)<0, 1);
    prepost_max_depth = find( (h.cfg.source.prepost_phase_amp_depth_perc + h.cfg.source.prepost_phase_amp_depth_perc_range)>100, 1);
    
    if ~isempty(sig_min_depth); w = warndlg(sprintf('\nSignal: Modulator Depth - Depth Range is < 0%%\n\nPlease make sure values fall within 0 - 100%%\n')); w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left'; end
    if ~isempty(sig_max_depth); w = warndlg(sprintf('\nSignal: Modulator Depth + Depth Range is > 100%%\n\nPlease make sure values fall within 0 - 100%%\n')); w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left'; end
    if ~isempty(prepost_min_depth); w = warndlg(sprintf('\nPrePost: Modulator Depth - Depth Range is < 0%%\n\nPlease make sure values fall within 0 - 100%%\n')); w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left'; end
    if ~isempty(prepost_max_depth); w = warndlg(sprintf('\nPrePost: Modulator Depth + Depth Range is > 100%%\n\nPlease make sure values fall within 0 - 100%%\n')); w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left'; end
    
elseif h.load_study_flag == 1   % loading study
    %% Updating Carrier & Modulator Menu items based on frequencies selected for TFR ROI
    if sum(sum(sum(h.cfg.source.sig_freqs)))>0
        for a=1:6
            if size(h.cfg.source.sig_freqs,2)>1
                sig_freqs = num2str (squeeze(h.cfg.source.sig_freqs(1,:,:)));
            else
                sig_freqs = num2str (squeeze(h.cfg.source.sig_freqs(1,:,:))');
            end
            
            h.menu_PAC_sig_carrier_source_freq(a).String = sig_freqs;
            h.menu_PAC_sig_modulator_source_freq(a).String = sig_freqs;
            h.menu_PAC_prepost_carrier_source_freq(a).String = sig_freqs;
            h.menu_PAC_prepost_modulator_source_freq(a).String = sig_freqs;
        end
    else
        fprintf('No ROI frequencies selected\n');
    end
    for a = 1:6
        % Signal freq
        xidx = find(h.cfg.source.sig_phase_amp_freq_idx(a,:)~=0);
        if ~isempty(xidx); h.menu_PAC_sig_carrier_source_freq(a).Value = xidx; % a = PAC contrast, h.menu_PAC_sig_carrier_source_freq(a).Value=index to set the number for the modulator index that is set by the h.menu_PAC_sig_modulator_source_freq(a).Value
        else; h.menu_PAC_sig_carrier_source_freq(a).Value=1; h.menu_PAC_sig_carrier_source_freq(a).String={'none'};
        end
        
        if ~isempty(h.cfg.source.sig_phase_amp_freq_idx(a,h.menu_PAC_sig_carrier_source_freq(a).Value)) && h.cfg.source.sig_phase_amp_freq_idx(a,h.menu_PAC_sig_carrier_source_freq(a).Value)~=0
            h.menu_PAC_sig_modulator_source_freq(a).Value = h.cfg.source.sig_phase_amp_freq_idx(a,h.menu_PAC_sig_carrier_source_freq(a).Value); % a = PAC contrast, h.menu_PAC_sig_carrier_source_freq(a).Value=index to set the number for the modulator index that is set by the h.menu_PAC_sig_modulator_source_freq(a).Value
        else; h.menu_PAC_sig_modulator_source_freq(a).Value = 1; h.menu_PAC_sig_modulator_source_freq(a).String = {'none'};
        end
        % Signal Modulator Depth
        h.edit_PAC_modulator_sig_depth(a).String{:} = num2str(h.cfg.source.sig_phase_amp_depth_perc(a,h.menu_PAC_sig_carrier_source_freq(a).Value)); % index of carrier frequency relative to selected TFR ROI
        % Signal Modulator Depth Range
        h.edit_PAC_modulator_sig_depth_range(a).String{:} = num2str(h.cfg.source.sig_phase_amp_depth_perc_range(a,h.menu_PAC_sig_carrier_source_freq(a).Value)); % index of carrier frequency relative to selected TFR ROI
        
        % PrePost freq
        xidx = find(h.cfg.source.prepost_phase_amp_freq_idx(a,:)~=0);
        if ~isempty(xidx); h.menu_PAC_prepost_carrier_source_freq(a).Value   = xidx; % index of carrier frequency relative to selected TFR ROI
        else h.menu_PAC_prepost_carrier_source_freq(a).Value=1; h.menu_PAC_prepost_carrier_source_freq(a).String = {'none'};
        end
        if ~isempty(h.cfg.source.prepost_phase_amp_freq_idx(a,h.menu_PAC_prepost_carrier_source_freq(a).Value)) && h.cfg.source.prepost_phase_amp_freq_idx(a,h.menu_PAC_prepost_carrier_source_freq(a).Value)~=0
            h.menu_PAC_prepost_modulator_source_freq(a).Value = h.cfg.source.prepost_phase_amp_freq_idx(a,h.menu_PAC_prepost_carrier_source_freq(a).Value); % index of carrier frequency relative to selected TFR ROI
        else; h.menu_PAC_prepost_modulator_source_freq(a).Value = 1; h.menu_PAC_prepost_modulator_source_freq(a).String = {'none'};
        end
        % PrePost Modulator Depth
        h.edit_PAC_modulator_prepost_depth(a).String{:} = num2str(h.cfg.source.prepost_phase_amp_depth_perc(a,h.menu_PAC_prepost_carrier_source_freq(a).Value)); % index of carrier frequency relative to selected TFR ROI
        % PrePost Modulator Depth Range
        h.edit_PAC_modulator_prepost_depth_range(a).String{:} = num2str(h.cfg.source.prepost_phase_amp_depth_perc_range(a,h.menu_PAC_prepost_carrier_source_freq(a).Value)); % index of carrier frequency relative to selected TFR ROI
    end
    
    
    
end
function btn_update_PAC_waves(varargin)
global h

update_source_cfg();
try
    % % Plot PAC waves
    for a = 1:size(h.PAC_source_contrasts,1)
        %% Signal waves
        % Carrier wave
        fc_idx = find(h.cfg.source.sig_phase_amp_freq_idx(a,:)>0);  % carrier freq index
        fm_idx = h.cfg.source.sig_phase_amp_freq_idx(a,fc_idx);     % modulator freq index
        sig_win = h.tfr_ROI(h.PAC_source_contrasts(a,1)).h(fc_idx).UserData.roi.sig_win_final;
        sig_waves = h.tfr_ROI(h.PAC_source_contrasts(a,1)).h(fc_idx).UserData.waves; % carrier frequency
        fc_waves = sig_waves.*sig_win';
        % Modulator wave
        sig_win = h.tfr_ROI(h.PAC_source_contrasts(a,2)).h(fm_idx).UserData.roi.sig_win_final;
        sig_waves = h.tfr_ROI(h.PAC_source_contrasts(a,2)).h(fm_idx).UserData.waves; % carrier frequency
        fm_waves = sig_waves.*sig_win';
        sig_win_ap = ( fm_waves ./max(max( abs(fm_waves ) ))); % renormalizing back to 0 to 1 for applying phase_amp windowing
        sig_win_ap = (sig_win_ap+1)/2; % reset range 0 to 1 so that troughs =0 not -1
        depth_range = h.cfg.source.sig_phase_amp_depth_perc_range(a,fc_idx);
        rn_range=randi([-depth_range depth_range],1)/100; % randomizing window amplitudes across trials
        rn_perc=(h.cfg.source.sig_phase_amp_depth_perc(a,fc_idx)/100)+rn_range;
        rn_diff=1-rn_perc;
        sig_win_ap = sig_win_ap.*rn_perc;
        sig_win_ap=sig_win_ap+rn_diff;
        % modulating carrier by modulator
        if sum(sum(sum(fc_waves)))==0 % carrier waves = 0
            pac_waves = zeros(size(fc_waves));
            fm_waves = zeros(size(fm_waves));
        elseif sum(sum(isnan(sig_win_ap)))>0 && sum(sum(fc_waves))~=0 % nans in windowing data
            pac_waves = fc_waves;
        elseif sum(sum(size(fc_waves)-size(sig_win_ap)))==0 % fc_waves and fm_waves have different sizes likely in frequency
            pac_waves = fc_waves.*sig_win_ap;
        elseif sum(sum(size(fc_waves)-size(sig_win_ap)))~=0 % fc_waves and fm_waves have different sizes likely in frequency
            w = warndlg(sprintf('\nCarrier and Modulater waves must have same number \nof frequencies (e.g., [10 10])'));
            w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
            pac_waves = zeros(size(fc_waves));
            fm_waves = zeros(size(fm_waves));
        else
            pac_waves = zeros(size(fc_waves));
            fm_waves = zeros(size(fm_waves));
        end  % in case that users zero the activity for this frequency at oneo of the sources
        
        % plotting Signal PAC waves
        h.axes_PAC_waves_sig(a).NextPlot = 'replace';
        plot(h.axes_PAC_waves_sig(a), h.cfg.study.lat_sim, pac_waves,'Color',h.src_clr(h.PAC_source_contrasts(a,1),:) );
        h.axes_PAC_waves_sig(a).NextPlot = 'add';
        plot(h.axes_PAC_waves_sig(a), h.cfg.study.lat_sim, fm_waves,'Color',h.src_clr(h.PAC_source_contrasts(a,2),:) , 'linewidth',1);
        plot(h.axes_PAC_waves_sig(a), h.cfg.study.lat_sim, envelope(pac_waves),'Color',h.src_clr(h.PAC_source_contrasts(a,1),:) , 'linewidth',2);
        title(h.axes_PAC_waves_sig(a),sprintf('Signal: Source %.f modulated by Source %.f ', h.PAC_source_contrasts(a,:)));
        h.axes_PAC_waves_sig(a).XLim = str2num(h.edit_plot_time_int.String); h.axes_PAC_waves_sig(a).YLim=[-1.1 1.1];
        h.axes_PAC_waves_sig(a).YLim = str2num(h.edit_plot_caxis.String)/100;
        %% PrePost waves
        % Carrier wave
        fc_idx = find(h.cfg.source.prepost_phase_amp_freq_idx(a,:)>0);  % carrier freq index
        fm_idx = h.cfg.source.prepost_phase_amp_freq_idx(a,fc_idx);     % modulator freq index
        prepost_win = h.tfr_ROI(h.PAC_source_contrasts(a,1)).h(fc_idx).UserData.roi.prepost_win_final;
        prepost_waves = h.tfr_ROI(h.PAC_source_contrasts(a,1)).h(fc_idx).UserData.waves; % carrier frequency
        fc_waves = prepost_waves.*prepost_win';
        % Modulator wave
        prepost_win = h.tfr_ROI(h.PAC_source_contrasts(a,2)).h(fm_idx).UserData.roi.prepost_win_final;
        prepost_waves = h.tfr_ROI(h.PAC_source_contrasts(a,2)).h(fm_idx).UserData.waves; % carrier frequency
        fm_waves = prepost_waves.*prepost_win';
        prepost_win_ap = ( fm_waves ./max(max( abs(fm_waves ) ))); % renormalizing back to 0 to 1 for applying phase_amp windowing
        prepost_win_ap = (prepost_win_ap+1)/2; % reset range 0 to 1 so that troughs =0 not -1
        depth_range = h.cfg.source.prepost_phase_amp_depth_perc_range(a,fc_idx);
        rn_range=randi([-depth_range depth_range],1)/100; % randomizing window amplitudes across trials
        rn_perc=(h.cfg.source.prepost_phase_amp_depth_perc(a,fc_idx)/100)+rn_range;
        rn_diff=1-rn_perc;
        prepost_win_ap = prepost_win_ap.*rn_perc;
        prepost_win_ap=prepost_win_ap+rn_diff;
        % modulating carrier by modulator
        % modulating carrier by modulator
        if sum(fc_waves)==0 % carrier waves = 0
            pac_waves = zeros(size(fc_waves));
            fm_waves = zeros(size(fm_waves));
        elseif sum(sum(isnan(prepost_win_ap)))>0 && sum(sum(fc_waves))~=0 % nans in windowing data
            pac_waves = fc_waves;
        else
            pac_waves = fc_waves.*prepost_win_ap;
        end  % in case that users zero the activity for this frequency at oneo of the sources
        
        % plotting PrePost PAC waves
        h.axes_PAC_waves_prepost(a).NextPlot = 'replace';
        plot(h.axes_PAC_waves_prepost(a), h.cfg.study.lat_sim, pac_waves,'Color',h.src_clr(h.PAC_source_contrasts(a,1),:) );
        h.axes_PAC_waves_prepost(a).NextPlot = 'add';
        plot(h.axes_PAC_waves_prepost(a), h.cfg.study.lat_sim, fm_waves,'Color',h.src_clr(h.PAC_source_contrasts(a,2),:) , 'linewidth',1);
        plot(h.axes_PAC_waves_prepost(a), h.cfg.study.lat_sim, envelope(pac_waves),'Color',h.src_clr(h.PAC_source_contrasts(a,1),:) , 'linewidth',2);
        title(h.axes_PAC_waves_prepost(a),sprintf('Prepost: Source %.f modulated by Source %.f ', h.PAC_source_contrasts(a,:)));
        h.axes_PAC_waves_prepost(a).XLim = str2num(h.edit_plot_time_int.String); h.axes_PAC_waves_prepost(a).YLim=[-1.1 1.1];
        h.axes_PAC_waves_prepost(a).YLim = str2num(h.edit_plot_caxis.String)/100;
        
    end
    h.axes_PAC_waves_sig(1).YLabel.String = 'Power (%/100)';
    h.axes_PAC_waves_sig(3).YLabel.String = 'Power (%/100)';
    h.axes_PAC_waves_sig(5).YLabel.String = 'Power (%/100)';
catch me
    fprintf('Error in "btn_update_PAC_waves"\n%s\n',me.message);
end


%% %%%%%%% "Simulate M/EEG" Panel Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% Simulate M/EEG Panel
%% Anatomy
function load_default_anatomy(varargin)
global h
% h.anatomy = struct('mri','','sens_eeg','','sens_meg','','mesh_volumes',[],'headmodel_eeg_vol','','headmodel_meg_vol','','headmodel_eeg_cortex','','headmodel_meg_cortex','','leadfield_eeg_vol','','leadfield_eeg_cortex','','leadfield_meg_vol','','leadfield_meg_cortex','','sens','','headmodel','','leadfield','','mesh_cortex','');
if h.start_flag == 1 % starting up SimMEEG -loads default anatomy
    % fpath = 'C:\BRANELab\matlab_progs\general_progs\EEG_sim\SimSignals_GUI\anatomy';
    %     h.anat_path = 'F:\brainstorm_db\SimMEEG_anat\data\S001_MEG_EEG\@default_study\';
    fpath = h.anat_path;
    h.anat_file = 'ANATOMY_DEFAULT_adult_SimMEEG_v19.mat';
    if exist(fullfile(fpath,h.anat_file),'file')
        h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Loading Default Anatomy\n\n%s\n',h.anat_file); %drawnow;
    else
        w = warndlg(sprintf('\nDefault Anatomy File: %s\nis not within selected Anatomy folder:\n%s\n\n',h.anat_file,h.anat_path));
        w.Position(3)=450; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left';
        [fname,fpath]=uigetfile('.mat','Load SimMEEG structured Anatomy');
        h.anat_file = fname;
        h.anat_path = fpath;
    end
elseif h.start_flag == 0
    [fname,fpath]=uigetfile('.mat','Load SimMEEG structured Anatomy');
    if fname~=0
        h.anat_file = fname; h.anat_path = fpath;
        h.pwd = pwd; cd(fpath);
        h.waitfor_panel.Visible='on';  h.waitfor_txt.String = sprintf('Loading Default Anatomy\n\n%s\n',h.anat_file); drawnow;
    else
    end
    
end
% try
fprintf('loading Anatomy File: %s\n',h.anat_file);


try
    h.anatomy=load(fullfile(h.anat_path,h.anat_file),'mri','sens_eeg','sens_meg','sens_meg','mesh_volumes','headmodel_eeg_vol','headmodel_meg_vol',...
        'headmodel_eeg_cortex','headmodel_meg_cortex','leadfield_eeg_vol','leadfield_eeg_cortex','leadfield_meg_vol','leadfield_meg_cortex');
    
    if strcmp(h.anat_file,'ANATOMY_DEFAULT_adult_SimMEEG_v16.mat')  % reordering meshc_volumes because leadfiels in v16 were based on white matter (mesh_volumes(5)) and not pial surface (mesh_volumes(4))
        h.anatomy.mesh_cortex = h.anatomy.mesh_volumes(4);
        h.anatomy.mesh_volumes(4) = h.anatomy.mesh_volumes(5);
    end
    
    for a=1:length(h.anatomy.mesh_volumes); h.anatomy.mesh_volumes(a) = ft_convert_units(h.anatomy.mesh_volumes(a),'mm'); end
    %% convert all units to 'mm'
    if ~isempty(h.anatomy.sens_eeg);  h.anatomy.sens_eeg = ft_convert_units(h.anatomy.sens_eeg,'mm'); end
    if ~isempty(h.anatomy.sens_meg);  h.anatomy.sens_meg = ft_convert_units(h.anatomy.sens_meg,'mm'); end
    if ~isempty(h.anatomy.mri);       h.anatomy.mri = ft_convert_units(h.anatomy.mri,'mm'); end
    
    
    if ~isempty(h.anatomy.headmodel_eeg_cortex);    h.anatomy.headmodel_eeg_cortex = ft_convert_units(h.anatomy.headmodel_eeg_cortex,'mm'); end
    if ~isempty(h.anatomy.headmodel_eeg_vol);       h.anatomy.headmodel_eeg_vol = ft_convert_units(h.anatomy.headmodel_eeg_vol,'mm'); end
    if ~isempty(h.anatomy.headmodel_meg_cortex);    h.anatomy.headmodel_meg_cortex = ft_convert_units(h.anatomy.headmodel_meg_cortex,'mm'); end
    if ~isempty(h.anatomy.headmodel_meg_vol);       h.anatomy.headmodel_meg_vol = ft_convert_units(h.anatomy.headmodel_meg_vol,'mm'); end
    
    
    if ~isempty(h.anatomy.leadfield_eeg_cortex); h.anatomy.leadfield_eeg_cortex = ft_convert_units(h.anatomy.leadfield_eeg_cortex,'mm'); end
    if ~isempty(h.anatomy.leadfield_eeg_vol); h.anatomy.leadfield_eeg_vol = ft_convert_units(h.anatomy.leadfield_eeg_vol,'mm'); end
    if ~isempty(h.anatomy.leadfield_meg_cortex); h.anatomy.leadfield_meg_cortex = ft_convert_units(h.anatomy.leadfield_meg_cortex,'mm'); end
    if ~isempty(h.anatomy.leadfield_meg_vol); h.anatomy.leadfield_meg_vol = ft_convert_units(h.anatomy.leadfield_meg_vol,'mm'); end
    
    %% make sens compatbile with SimMEEG
    if ~isempty(h.anatomy.sens_meg)
        h.anatomy.sens_meg = sm_make_sens_compatible(h.anatomy.sens_meg,'meg');
        cfg.sens_type = 'meg'; cfg.sens_idx = 1:length(h.anatomy.sens_meg.label); h.anatomy.sens_meg = bs_select_meeg_sensors(cfg,h.anatomy.sens_meg);
    end
    if ~isempty(h.anatomy.sens_eeg)
        h.anatomy.sens_eeg = sm_make_sens_compatible(h.anatomy.sens_eeg,'eeg');
        cfg.sens_type = 'eeg'; cfg.sens_idx = 1:length(h.anatomy.sens_eeg.label); h.anatomy.sens_eeg = bs_select_meeg_sensors(cfg,h.anatomy.sens_eeg);
    end
    
    %% add BRANELab leadfield format
    if ~isempty(h.anatomy.leadfield_eeg_cortex)
        x=cell2mat(h.anatomy.leadfield_eeg_cortex.leadfield(h.anatomy.leadfield_eeg_cortex.inside==1));
        h.anatomy.leadfield_eeg_cortex.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
        h.anatomy.leadfield_eeg_cortex.voxel_pos=h.anatomy.leadfield_eeg_cortex.pos(h.anatomy.leadfield_eeg_cortex.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
        h.anatomy.leadfield_eeg_cortex.voxel_res=[]; % resolution not available for cortical surface
    end
    if ~isempty(h.anatomy.leadfield_meg_cortex)
        x=cell2mat(h.anatomy.leadfield_meg_cortex.leadfield(h.anatomy.leadfield_meg_cortex.inside==1));
        h.anatomy.leadfield_meg_cortex.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
        h.anatomy.leadfield_meg_cortex.voxel_pos=h.anatomy.leadfield_meg_cortex.pos(h.anatomy.leadfield_meg_cortex.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
        h.anatomy.leadfield_meg_cortex.voxel_res=[]; % resolution not available for cortical surface
    end
    if ~isempty(h.anatomy.leadfield_eeg_vol)
        x=cell2mat(h.anatomy.leadfield_eeg_vol.leadfield(h.anatomy.leadfield_eeg_vol.inside==1));
        h.anatomy.leadfield_eeg_vol.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
        h.anatomy.leadfield_eeg_vol.voxel_pos=h.anatomy.leadfield_eeg_vol.pos(h.anatomy.leadfield_eeg_vol.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
        vx_res = diff(h.anatomy.leadfield_eeg_vol.voxel_pos(:,1)); vx_res = min(abs(vx_res(vx_res~=0)));
        h.anatomy.leadfield_eeg_vol.voxel_res = vx_res;
    end
    if ~isempty(h.anatomy.leadfield_meg_vol)
        x=cell2mat(h.anatomy.leadfield_meg_vol.leadfield(h.anatomy.leadfield_meg_vol.inside==1));
        h.anatomy.leadfield_meg_vol.H=reshape(x,[size(x,1) 3 size(x,2)/3]);
        h.anatomy.leadfield_meg_vol.voxel_pos=h.anatomy.leadfield_meg_vol.pos(h.anatomy.leadfield_meg_vol.inside,:); % brain voxel locations --> these correspond to the 16008 positions for the leadfield.H
        vx_res = diff(h.anatomy.leadfield_eeg_vol.voxel_pos(:,1)); vx_res = min(abs(vx_res(vx_res~=0)));
        h.anatomy.leadfield_eeg_vol.voxel_res = vx_res;
    end
    
    %% default starting with MEG
    if ~isempty(h.anatomy.sens_meg)
        h.anatomy.sens = h.anatomy.sens_meg; h.menu_sens_type.Value = 1;
        if ~isempty(h.anatomy.headmodel_meg_vol)
            h.anatomy.headmodel = h.anatomy.headmodel_meg_vol;
            h.anatomy.leadfield = h.anatomy.leadfield_meg_vol;
        elseif ~isempty(h.anatomy.headmodel_meg_cortex)
            h.anatomy.headmodel = h.anatomy.headmodel_meg_cortex;
            h.anatomy.leadfield = h.anatomy.leadfield_meg_cortex;
        else
            h.anatomy.headmodel = [];
            h.anatomy.leadfield =[];
        end
    elseif ~isempty(h.anatomy.sens_eeg) && isempty(h.anatomy.sens_meg)
        h.anatomy.sens = h.anatomy.sens_eeg; h.menu_sens_type.Value = 1;
        if ~isempty(h.anatomy.headmodel_eeg_vol)
            h.anatomy.headmodel = h.anatomy.headmodel_eeg_vol;
            h.anatomy.leadfield = h.anatomy.leadfield_eeg_vol;
        elseif ~isempty(h.anatomy.headmodel_eeg_cortex)
            h.anatomy.headmodel = h.anatomy.headmodel_eeg_cortex;
            h.anatomy.leadfield = h.anatomy.leadfield_eeg_cortex;
        else
            h.anatomy.headmodel = [];
            h.anatomy.leadfield =[];
        end
    end
    
    
    if length(h.anatomy.mesh_volumes)>3
        h.anatomy.mesh_cortex = h.anatomy.mesh_volumes(5);
    else
        h.anatomy.mesh_cortex = h.anatomy.mesh_volumes(3);
    end
    update_anatomy_fields;
    
    
catch
    txt = [sprintf('\nPlease create SimMEEG Anatomy Structure, save as .mat file, and then load them into SimMEEG\n\nSee FieldTrip Tutorial: http://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_bem/ \n\n'),...
        sprintf('Make sure the anatomy .mat file has the following variables:\n\n'), ...
        sprintf('          mri\n'), sprintf('          sens_eeg\n'),sprintf('          sens_meg\n'),...
        sprintf('          headmodel_eeg_vol\n'), sprintf('          headmodel_meg_vol\n'), ...
        sprintf('          headmodel_eeg_cortex\n'), sprintf('          headmodel_meg_cortex\n'),...
        sprintf('          mesh_volumes\n'),...
        sprintf('          leadfield_eeg_vol\n'), sprintf('          leadfield_eeg_cortex\n'), sprintf('          leadfield_meg_vol\n'),...
        sprintf('          leadfield_meg_cortex\n')];
    
    %% Clearing Anatomy Strings and Turning off Head Model menu
    h.anatomy_file_txt.String = sprintf('%s',h.anat_file);
    h.mri_txt.String = sprintf('No Anatomy Loaded');
    h.sensors_txt.String = sprintf('Missing Anatomy Fields');
    h.headmodel_txt.String = sprintf('or wrong formats');
    h.leadfield_txt.String = sprintf(' ');
    h.menu_head_model.Enable = 'inactive';
    h.menu_sens_type.Enable = 'inactive';
    h.menu_head_model.Value = 1;
    h.menu_head_model.BackgroundColor = [1 1 1]*.9;  h.menu_head_model.ForegroundColor = [1 1 1]*.4;
    
    w = warndlg(txt);
    w.Position=[h.main_fig.Position(1:2)+100 500 350]; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left';
end
try plot_3D_mri; catch; end
h.listbox_chans.Value=1;

h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
function load_anatomy(varargin)
global h
cd(h.anat_path);
[fname,fpath]=uigetfile('.mat','MRI File');
% [fdir,fname2,fext]=fileparts(fname);
h.anat_file = fname;
h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Loading Anatomy\n\n%s\n',h.anat_file); drawnow;
%% store in SimMEEG naming format
try
    h.pwd = pwd; cd(fpath);
    fprintf('loading Default Anatomy File: %s\n',h.anat_file);
    rmfield(h,'anatomy');
    h.anatomy=load(h.anat_file,'mri','please_read','sens_eeg','sens_meg','headmodel_eeg','headmodel_meg','mesh_volumes',...
        'leadfield_eeg_vol','leadfield_eeg_cortex','leadfield_meg_vol','leadfield_meg_cortex');
    
    % default starting with EEG
    h.anatomy.sens = h.anatomy.sens_eeg;
    h.anatomy.headmodel = h.anatomy.headmodel_eeg;
    h.anatomy.leadfield = h.anatomy.leadfield_eeg_vol;
    h.anatomy.mesh_cortex = h.anatomy.mesh_volumes(5);
    
    
catch
    
    txt = [sprintf('Please create Anatomy Structure in Field Trip, save as mat file, and then load them into SimMEEG\n\nSee FieldTrip Tutorial: http://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_bem/ \n\n'),...
        sprintf('Make sure the anatomy .mat file has the following variables:\n\n'), ...
        sprintf('          mri\n'), sprintf('          sens_eeg\n'),sprintf('          sens_meg\n'),...
        sprintf('          headmodel_eeg_vol\n'), sprintf('          headmodel_meg_vol\n'), ...
        sprintf('          headmodel_eeg_cortex\n'), sprintf('          headmodel_meg_vol\n'),...
        sprintf('          mesh_volumes\n'),...
        sprintf('          leadfield_eeg_vol\n'), sprintf('          leadfield_eeg_cortex\n'), sprintf('          leadfield_meg_vol\n'),...
        sprintf('          leadfield_meg_cortex\n')];
    
    %% Clearing Anatomy Strings and Turning off Head Model menu
    h.mri_txt.String = sprintf('No Anatomy Loaded');
    h.headmodel_txt.String = sprintf(' ');
    h.leadfield_txt.String = sprintf(' ');
    h.sensors_txt.String = sprintf(' ');
    h.menu_head_model.Enable = 'inactive';
    h.menu_sens_type.Enable = 'inactive';
    h.menu_head_model.Value = 1;
    h.menu_head_model.BackgroundColor = [1 1 1]*.9;  h.menu_head_model.ForegroundColor = [1 1 1]*.4;
    
    warndlg(txt);
    w.Position=[h.main_fig.Position(1:2)+100 500 350]; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left';
    
end
update_anatomy_fields;
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
function set_anat_file(varargin)
global h
h.anat_file = h.anatomy_file_txt.String;
function sm_load_real_sensor_file(varargin)
% this program loads in real data from a .mat file with the following variables:
%
%   sensor_file.mat
%       data = [sensor x samples x trials]
%       srate = sample rate
%       lat = latency values of samples
global h

study_trials = str2num(h.edit_num_trials.String);
study_srate = str2num(h.edit_srate.String);

if h.monte_carlo_flag == 1 && h.menu_monte_synthetic_real_data.Value==4  % load Real Sensor Data in Monte-Carlo sims
    [fpath,fname] = fileparts(h.monte_real_dsname);
    h.waitfor_txt.String = sprintf('Loading Real Sensor data from file: %s\n',fullfile(fpath,fname)); drawnow;
else                                                                     % load Real Sensor Data in manual mode
    [fname,fpath]=uigetfile('*.mat','Load Real Sensor Dataset');
    h.waitfor_panel.Visible='on';
    h.waitfor_txt.String = sprintf('Loading Real Sensor data from file: %s\n',fullfile(fpath,fname)); drawnow;
end

x = load(fullfile(fpath,fname),'data','srate','lat');
h.anatomy.sens.bad_sensors = find(isnan(squeeze(x.data(:,1,1)))==1);    % finding bad sensors for data with nans
update_anatomy_fields;
% if ~isempty(bad_chans) % checking for bad channels being set to nans
%     x.data = sm_interpolate_data(x.data,bad_chans);
% end


if x.srate > study_srate  % downsampling the data
    x.data = resample(double(x.data), study_srate , x.srate, 'Dimension',2);
    x.lat = resample(double(x.lat), study_srate , x.srate, 'Dimension',2);
    x.srate = study_srate;
    
    % trying to align timings - if not selecting first samples
    % see if all study.lat_sim are within x.lat values
    lat_sim = h.cfg.study.lat_sim;
    if x.lat(1)<lat_sim(1) && x.lat(end)>lat_sim(end)   % study.lat_sim lies within x.lat so data can be extracted
        xs = find(x.lat<lat_sim(1)); xs1 = xs(end);
        xs = find(x.lat<lat_sim(end)); xs2 = xs(end);
        x.data = x.data(:,xs1:xs2,:);
        x.lat = lat_sim;
    elseif length(x.lat) >= length(lat_sim)
        x.data = x.data(:,1:length(lat_sim),:);
        x.lat = lat_sim;
    end
else
end

%% randomly selecting number of trials
if x.srate==study_srate && size(x.data,2)<=length(h.cfg.study.lat_sim) && size(x.data,3)>=study_trials && size(x.data,1)>=length(h.anatomy.leadfield.label) % >= num of sensors used for leadfields
    
    nt = randperm(size(x.data,3));  % randomly selecting trials
    chan_idx = 1:length(h.anatomy.leadfield.label); % selecting first sensors
    xsamps = 1:length(h.cfg.study.lat_sim);     % selecting first samples
    fprintf('Sensor indices for sensor %s data: ',varargin{3}); fprintf('%.f ',chan_idx); fprintf('\n');
    fprintf('Randomly selecting %.f of %.f trials\n',study_trials,size(x.data,3))
    
    switch varargin{3}
        case 'noise'
            h.sim_data.sens_noise = permute(x.data(chan_idx,xsamps,nt(1:study_trials)),[2 1 3]);
            h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
            h.sim_data.sens_noise_final = h.sim_data.sens_noise_scaled;
            h.sim_data.sens_noise_final_org = h.sim_data.sens_noise_final;
            h.sim_data.sens_noise_type = h.menu_noise_projection.String{h.menu_noise_projection.Value};
        case 'final'
            h.sim_data.sens_final = permute(x.data(chan_idx,xsamps,nt(1:study_trials)),[2 1 3]);
            h.sim_data.sens_sig_data = zeros(size(h.sim_data.sens_final));
            h.sim_data.sens_noise_final = zeros(size(h.sim_data.sens_final));
            
            h.sim_data.sens_final_org = h.sim_data.sens_final;  % setting orignal data to newly load M/EEG dataset
            h.sim_data.sens_noise_final_org = h.sim_data.sens_noise_final;  % removing noise
            h.sim_data.sens_sig_data_org = h.sim_data.sens_sig_data;   % removing sensor data from projected sources
            h.sim_data.sens_noise_type      = 'None - Loaded in Real Full Data';
            h.sim_data.source_waveform_type = 'None - Loaded in Real Full Data';
    end
    
    
else
    s1 = sprintf('Some Data File parameters are NOT compatible with Study Parameters:\n\n                         Data File          Study Parameters\n');
    if x.srate~=study_srate; s2 = sprintf('Sample Rate    %.f       ~=      %s   (No)\n',x.srate,h.edit_srate.String);
    else; s2 = sprintf('Sample Rate     %.f       ==      %s   (Yes)\n',x.srate,h.edit_srate.String); end
    
    if size(x.data,2)~=length(h.cfg.study.lat_sim); s3 = sprintf('Num Samples   %.f       <       %.f   (No)\n',size(x.data,2),length(h.cfg.study.lat_sim));
    else;  s3 = sprintf('Num Samples   %.f       >=      %.f   (Yes)\n',size(x.data,2),length(h.cfg.study.lat_sim));  end
    
    if size(x.data,1)<length(h.anatomy.leadfield.label); s4 = sprintf('Num Sensors   %.f      <       %.f   (No)\n',size(x.data,1),length(h.anatomy.leadfield.label));
    else;  s4 = sprintf('Num Sensors    %.f      >=      %.f   (Yes)\n',size(x.data,1),length(h.anatomy.leadfield.label));   end
    
    if size(x.data,3)<study_trials; s5 = sprintf('Num Trials    %.f       <       %.f   (No)\n',size(x.data,3),study_trials);
    else;  s5 = sprintf('Num Trials        %.f       >=      %.f   (Yes)\n',size(x.data,3),study_trials);   end
    sx = [s1, s2, s3, s4, s5];
    w = warndlg(sprintf('%s\n',sx));
    w.Position(3:4)=[400 200]; htext = findobj(w, 'Type', 'Text'); htext.FontSize = 11; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
end

update_study_cfg();

if h.monte_carlo_flag ~= 1
    h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
end
function sm_load_source_waves(varargin)
% this program loads in source waveform data from a .mat file with the following variables:
%
%   source_wave.mat
%       source_data = [sources x samples x trials]
%       srate = sample rate
%       lat = latency values of samples

global h

study_trials = str2num(h.edit_num_trials.String);
study_srate = str2num(h.edit_srate.String);

if h.monte_carlo_flag == 1 && h.menu_monte_synthetic_real_source.Value==2  % load Real Source Data in Monte-Carlo sims
    [fpath,fname] = fileparts(h.monte_real_dsname);
    h.waitfor_txt.String = sprintf('Loading Source Waveform Data from file: %s\n',fullfile(fpath,fname)); drawnow;
else                                                                     % load Real Source Data in manual mode
    [fname,fpath]=uigetfile('*.mat','Load Source Waveform Data');
    h.waitfor_panel.Visible='on';
    h.waitfor_txt.String = sprintf('Loading Real Sensor data from file: %s\n',fullfile(fpath,fname)); drawnow;
end


x = load(fullfile(fpath,fname),'source_data','srate','lat');
if ~isfield(x,'source_data')
    warndlg(sprintf('Variable "source_data" not found in file: %s',fullfile(fpath,fname)));
else
    
    if x.srate > study_srate  % downsampling the data
        
        x.source_data = resample(double(x.source_data), study_srate , x.srate, 'Dimension',2);
        x.lat = resample(double(x.lat), study_srate , x.srate, 'Dimension',2);
        x.srate = study_srate;
        
        % trying to align timings - if not selecting first samples
        % see if all study.lat_sim are within x.lat values
        lat_sim = h.cfg.study.lat_sim;
        if x.lat(1)<lat_sim(1) && x.lat(end)>lat_sim(end)   % study.lat_sim lies within x.lat so data can be extracted
            xs = find(x.lat<lat_sim(1)); xs1 = xs(end);
            xs = find(x.lat<lat_sim(end)); xs2 = xs(end);
            x.source_data = x.source_data(:,xs1:xs2,:);
            x.lat = lat_sim;
        elseif length(x.lat) >= length(lat_sim)
            x.source_data = x.source_data(:,1:length(lat_sim),:);
            x.lat = lat_sim;
        end
    else
    end
    
    %% randomly selecting number of trials
    if x.srate==study_srate && size(x.source_data,2)<=length(h.cfg.study.lat_sim) && size(x.source_data,3)>=study_trials && size(x.source_data,1)>=size(h.cfg.source.vx_locs,1)
        nt = randperm(size(x.source_data,3));  % randomly selecting trials
        s_idx = 1:size(h.cfg.source.vx_locs,1) ; % selecting first sensors
        xsamps = 1:length(h.cfg.study.lat_sim);     % selecting first samples
        fprintf('Source indices for source waveforms: '); fprintf('%.f ',s_idx); fprintf('\n');
        fprintf('Randomly selecting %.f of %.f trials\n',study_trials,size(x.source_data,3))
        h.sim_data.sig_final = permute(x.source_data(s_idx,xsamps,nt(1:study_trials)),[2 1 3]);
        %     h.sim_data.sig_final = (h.sim_data.sig_final./max(max(max(abs(h.sim_data.sig_final))))); % scaled between -1 to 1
        h.sim_data.sig_wav = zeros(size(h.sim_data.sig_final));
        h.sim_data.sig_win = h.sim_data.sig_wav;
        h.sim_data.prepost_wav = h.sim_data.sig_wav;
        h.sim_data.prepost_win = h.sim_data.sig_wav;
        h.sim_data.source_waveform_type = 'Real Sources';
    else
        s1 = sprintf('Some Data File parameters are NOT compatible with Study Parameters:\n\n                         Data File          Study Parameters\n');
        if x.srate~=study_srate; s2 = sprintf('Sample Rate    %.f       ~=      %s   (No)\n',x.srate,h.edit_srate.String);
        else; s2 = sprintf('Sample Rate     %.f       ==      %s   (Yes)\n',x.srate,h.edit_srate.String); end
        if size(x.source_data,2)~=length(h.cfg.study.lat_sim); s3 = sprintf('Num Samples   %.f       <       %.f   (No)\n',size(x.source_data,2),length(h.cfg.study.lat_sim));
        else;  s3 = sprintf('Num Samples   %.f       >=      %.f   (Yes)\n',size(x.source_data,2),length(h.cfg.study.lat_sim));  end
        if size(x.source_data,1)~=size(h.cfg.source.vx_locs,1); s4 = sprintf('Num Sources   %.f      <      %.f   (No)\n',size(x.source_data,1),size(h.cfg.source.vx_locs,1));
        else;  s4 = sprintf('Num Sensors    %.f      <=      %.f   (Yes)\n',size(x.source_data,1),size(h.cfg.source.vx_locs,1));   end
        if size(x.source_data,3)<study_trials; s5 = sprintf('Num Trials    %.f       <       %.f   (No)\n',size(x.source_data,3),study_trials);
        else;  s5 = sprintf('Num Trials        %.f       >=      %.f   (Yes)\n',size(x.source_data,3),study_trials);   end
        sx = [s1, s2, s3, s4, s5];
        w = warndlg(sprintf('%s\n',sx));
        w.Position(3:4)=[400 200]; htext = findobj(w, 'Type', 'Text'); htext.FontSize = 11; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
    end
    % h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
    update_study_cfg();
end

if h.monte_carlo_flag ~= 1
    h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
end
function sm_fft_randomize_phase(varargin)
% use Alex Moiseev's program that will shuffle the phase using FFT of PCA components then put it back together as [sensors x samples x trials]
global h

if h.monte_carlo_flag == 1
    h.waitfor_txt.String = sprintf('Running FFT randomization of phases for sensor noise data.\n'); drawnow;
else
    h.waitfor_panel.Visible='on';
    h.waitfor_txt.String = sprintf('Running FFT randomization of phases for sensor noise data.\n'); drawnow;
end

xdata = permute(h.sim_data.sens_noise(:,h.anatomy.sens.good_sensors,:),[2 1 3]);    % selecting only good channels
dims = size(xdata);
xdata = reshape(xdata,[dims(1) dims(2)*dims(3)]);

ydata = genSurrogateData(xdata);
h.sim_data.sens_noise(:,h.anatomy.sens.good_sensors,:) = permute(reshape(ydata,dims),[2 1 3]);
h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1

if h.monte_carlo_flag ~= 1
    h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
end
function menu_filt_type_Callback(varargin)
global h

if h.menu_filt_type.Value <=4 && str2num(h.edit_filt_order.String)==0 % FIR filter with no filter order
    h.menu_filt_method.String = {'equiripple'; 'kaiserwin'};
elseif h.menu_filt_type.Value <=4 && str2num(h.edit_filt_order.String)~=0 % FIR filter with no filter order
    h.menu_filt_method.String = {'equiripple'; 'maxflat'; 'window'; 'ls'};
elseif h.menu_filt_type.Value >4 % IIR filter with no filter order
    h.menu_filt_method.String = {'butter';    'cheby1';    'cheby2';    'ellip';};
end
if h.menu_filt_method.Value>length(h.menu_filt_method.String)
    h.menu_filt_method.Value=1;
end
function radio_panel_preprocess_CallBack(varargin)
global h
h.panel_preprocess.Visible = 'off';
if h.radio_panel_preprocess.Value == 1 % Turn on Filtering Panel
    h.panel_preprocess.Visible = 'on';
end
%% %%%% Project Source Data --> Sensors
%% Set Visibility of Brain Noise panel
function menu_noise_projection_CallBack(varargin)
global h
h.cfg.study.noise_projection_type = h.menu_noise_projection.String{h.menu_noise_projection.Value};
h.btn_sim_sens_noise.Visible = 'off';
h.panel_brain_noise.Visible = 'off';
h.btn_get_real_noise_dataset.Visible = 'off';
h.btn_fft_radnomize_phase.Visible = 'off';

if h.menu_noise_projection.Value==1 % Sensor Noise
    h.btn_sim_sens_noise.Visible = 'on';
    h.btn_fft_radnomize_phase.Visible = 'on';
elseif h.menu_noise_projection.Value==2 % Brain Noise
    h.btn_sim_sens_noise.Visible = 'on';
    h.panel_brain_noise.Visible = 'on';
    h.btn_fft_radnomize_phase.Visible = 'on';
elseif h.menu_noise_projection.Value==3 % Synthetic MVAR or GAN (Generative Adversarial Network) - machine learning <-- Not implemented yet
    % add in function in future updates
elseif h.menu_noise_projection.Value==4 % Real Sensor
    h.btn_get_real_noise_dataset.Visible = 'on';
    h.btn_sim_sens_noise.Visible = 'off';
    h.btn_fft_radnomize_phase.Visible = 'on';
end
%% Set Real Noise Directory
function set_real_noise_dir(varargin)
global h
%% Simulate Sensor Noise
function sim_sens_noise(varargin)
global h
if h.monte_carlo_flag == 1
    h.waitfor_txt.String = sprintf('Simulating Sensor Noise.\n\nThis might take some time.\n\nPlease wait ...'); drawnow;
    h.sim_data.sens_noise_type = h.menu_monte_synthetic_real_data.String{h.menu_monte_synthetic_real_data.Value};
else
    h.waitfor_panel.Visible='on';
    h.waitfor_txt.String = sprintf('Simulating Sensor Noise.\n\nThis might take some time.\n\nPlease wait ...'); drawnow;
    h.sim_data.sens_noise_type = h.menu_noise_projection.String{h.menu_noise_projection.Value};
end


if h.menu_noise_projection.Value == 1   % Simulate Sensor noise for all sensors
    %% NOTE: Simulating sensor noise is quite complicated because of spatial and temporal covariance needs to be modeled --> see reading in GUI directory: DeMunckEtAl_StatDip_KronPrud_IEEE_SP_2002_spatiotemporal_noise_covariance_simulation.pdf
    % Currently not implementing covariance structure (see 'underconstruction code' below as a start).
    h.sens_noise_cov_flag=0;
    if h.sens_noise_cov_flag==0 % no shaping by spatial or temporal covariance
        h.sim_data.sens_noise = sim_noise_wav(h.cfg,size(h.anatomy.leadfield.H,1));
        h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
        %         for t=randperm(h.cfg.study.num_trials); figure(3); clf; surf(cov(squeeze(h.sim_data.sens_noise_scaled(:,:,t(1))))); view(0,90); shading interp; axis tight; caxis([-.045 .045]); pause; end
        %% NOT IMPLEMENTED simulating sensor noise with spatiotemporal covariance --> "spatio-temporal noise covariance is a Kronecker product of a spatial and a temporal covariance matrix"
    elseif h.sens_noise_cov_flag==1     % Spatial Covariance shaping only
        
        noiseVec2 = sim_noise_wav(h.cfg,size(h.anatomy.leadfield.H,1));
        
        h.sensor_noise_cov_exp = (1/3.57^1.02);
        
        %         noiseVec2 = noiseVec2/max(max(max(abs(noiseVec2))));
        noiseVec = permute(noiseVec2,[2 1 3]);
        for t=1:size(noiseVec,3)     % randomizing covariance structure across trials
            c_cov = random('exp',h.sensor_noise_cov_exp,size(noiseVec,1),size(noiseVec,1)); % spatial covariances 'exp' (1/3.57^1.02) not accounting for distance between channels so loosely based on Huizenga et al., 2002 IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. 49, NO. 6, pp. 533:539; see Huizenga_2002_noise_cov_randomization_exp_values.pdf
            noiseVec(:,:,t) = c_cov*squeeze(noiseVec(:,:,t));
        end
        h.sim_data.sens_noise = permute(noiseVec,[2 1 3]);
        h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
        %         for t=randperm(h.cfg.study.num_trials); figure(3); clf; surf(cov(squeeze(h.sim_data.sens_noise_scaled(:,:,t(1))))); view(0,90); shading interp; axis tight; caxis([-55 55]); pause; end
        %% NOT IMPLEMENTED simulating sensor noise with spatiotemporal covariance --> "spatio-temporal noise covariance is a Kronecker product of a spatial and a temporal covariance matrix"
    elseif h.sens_noise_cov_flag==2     % Spatial and Temporal Covariance shaping <-- NOT IMPLEMENTED yet
        
    end
    
    
elseif h.menu_noise_projection.Value == 2   % Simulate Brain noise
    if ~isfield(h.cfg.source,'brain_noise_idx') % in case "Random locs" button has yet to be pressed
        btn_rand_brain_noise_locs_Callback([],[]);
    end
    sim_brain_noise_wav;
    %% projecting brain source locations across trials
    h.sim_data.sens_noise=[];
    for t=1:h.cfg.study.num_trials
        h.sim_data.sens_noise(:,:,t) = project_SimSignals(h.sim_data.brain_noise(:,:,t),h.anatomy.leadfield,h.cfg.source.brain_noise_idx(:,t),h.cfg.source.brain_noise_amp(:,t),h.cfg.source.brain_noise_ori(:,:,t));
    end
    h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
    % validating code
    %     figure(1); clf; hold on; for t=1:h.cfg.study.num_trials; plot(h.cfg.study.lat_sim,squeeze(h.sim_data.sens_noise_scaled(:,:,t)),'k'); end
    %     figure(2); clf; calc_fft(squeeze(h.sim_data.sens_noise(:,48,:)),256,1,'k');
    %         for t=randperm(h.cfg.study.num_trials); figure(3); clf; surf(cov(squeeze(h.sim_data.sens_noise_scaled(:,:,t(1))))); view(0,90); shading interp; axis tight; caxis([-.025 .025]); pause; end
    %     scatter3(h.axes_anatomy,h.anatomy.leadfield.voxel_pos(h.cfg.source.brain_noise_idx,1),h.anatomy.leadfield.voxel_pos(h.cfg.source.brain_noise_idx,2),h.anatomy.leadfield.voxel_pos(h.cfg.source.brain_noise_idx,3),'k.');
elseif h.menu_noise_projection.Value == 3       % Generative Adversarial Networks (GAN) generated noise <-- currently must have srate=256 and duration -1 to 1 sec with chans = 1:66
    implemented_flag = 0;    % Not implmented yet
    if implemented_flag ==1
        rest_dir = 'C:\BRANELab\matlab_progs\general_progs\EEG_sim\SimSignals_GUI\GANsimEEG\Trained_GANs';
        sname = fullfile(rest_dir,sprintf('Tujillo_Resting_State_%s_GAN_trained.mat',blk_name{h.menu_resting_state.Value}));
        fprintf('GAN generated noise using:  %s\n',sname);
        load(sname,'dlnetGenerator');
        
        %% Generate EEG sensor noise using pre-trained GAN
        
        executionEnvironment = "auto";
        ZZNew = randn(1,1,1000,h.cfg.study.num_trials,'single');
        dlZZNew = dlarray(ZZNew,'SSCB');
        if (executionEnvironment == "auto" && canUseGPU) || executionEnvironment == "gpu"
            dlZZNew = gpuArray(dlZZNew);
        end
        
        dlXGeneratedNew = predict(dlnetGenerator,dlZZNew);
        y = gather(extractdata(dlXGeneratedNew));
        ydata = permute(squeeze(y(:,:,:,:)),[2 1 3]);
        ydata = ydata(1:size(h.anatomy.sens.label),:,:);
        ydata = bl_reref(ydata,[1:size(ydata,1)],[],1);    % average referencing
        h.sim_data.sens_noise = permute(ydata,[2 1 3]);
        h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
    end
elseif h.menu_noise_projection.Value == 4        % Real Sensor Noise
    
end
if h.monte_carlo_flag ~= 1
    h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
end
fprintf('Finished Simulating "%s" Noise\n',h.menu_noise_projection.String{h.menu_noise_projection.Value});
%% simulate Brain Noise
function sim_brain_noise_wav
global h
num_chans = str2num(h.edit_num_noise_sources.String);
% broadband noise
h.cfg.study.noise_flag = h.edit_noise_flag.Value;    % 1=Broadband, 2=NarrowBand, 3=Notched Band, 4=Pink
h.sim_data.brain_noise = sim_noise_wav(h.cfg,num_chans);
h.sim_data.brain_noise = h.sim_data.brain_noise/max(max(max(abs(h.sim_data.brain_noise))));
% figure(1); clf; plot(h.cfg.study.lat_sim,squeeze(h.sim_data.brain_noise(:,:,1))); figure(2); clf; calc_fft(squeeze(h.sim_data.brain_noise(:,:,1)),h.cfg.study.srate,1,'k');
%% Randomly set Brain Noise Locs
function btn_rand_brain_noise_locs_Callback(varargin)
global h
% Fully randomize brain noise locations, orientations, and amplitudes for each trial
num_noise = str2num(h.edit_num_noise_sources.String);
brn_vx = setdiff(1:size(h.anatomy.leadfield.voxel_pos,1),h.cfg.source.vx_idx); % removing true source voxels from being included in possible brain noise voxels
% random noise amps
noise_amp = str2num(h.edit_brain_noise_amp.String);
h.cfg.source.brain_noise_idx = zeros(num_noise,h.cfg.study.num_trials);
h.cfg.source.brain_noise_amp = h.cfg.source.brain_noise_idx;
h.cfg.source.brain_noise_ori=ones(num_noise,3,h.cfg.study.num_trials);

rn = randperm(length(brn_vx)); % same locs across trials
%     az = deg2rad(normrnd(0,180,num_noise,1));
%     el = deg2rad(normrnd(0,180,num_noise,1));
for t=1:h.cfg.study.num_trials % new set of brain noise sources each trial
    % random locs across trials
    %      rn = randperm(length(brn_vx));
    h.cfg.source.brain_noise_idx(:,t) = brn_vx(rn(1:num_noise)); % brain noise voxel locations
    h.cfg.source.brain_noise_amp(:,t) = abs(normrnd(noise_amp(1),noise_amp(2),num_noise,1));  % randomly sampling ampltidues for noise sources from a mean = 20 nA and stdev = 10 <-- This si completely arbitirary right now.
    
    % random noise orientations
    if h.menu_ori_normal.Value==1  % Volume
        az = deg2rad(normrnd(0,180,num_noise,1));
        el = deg2rad(normrnd(0,180,num_noise,1));
        [x,y,z]=sph2cart(az,el,ones(size(az)));
        h.cfg.source.brain_noise_ori(:,:,t) = [x,y,z];
    elseif h.menu_ori_normal.Value==2  % Cortically Constrained
        h.cfg.source.brain_noise_ori(:,:,t) = h.anatomy.leadfield_eeg_cortex.ori(h.cfg.source.brain_noise_idx(:,t),:);
    end
    
    
    
    
    
end
%% Load Brain Noise Locs from preset file
function btn_load_brain_noise_locs_Callback(varargin)
global h
msgbox('Not Implemented Yet');
%% Use preset Default mode network locations - randomized slightly across trials
function btn_load_default_mode_noise_locs_Callback(varargin)
global h
msgbox('Not Implemented Yet');
%% plot noise locations on Anatomy
function btn_plot_noise_locs(varargin)
global h
% plot_3D_mri('off'); % reploting anatomy + sources without scalp
scatter3(h.axes_anatomy,h.anatomy.leadfield.voxel_pos(h.cfg.source.brain_noise_idx,1),h.anatomy.leadfield.voxel_pos(h.cfg.source.brain_noise_idx,2),h.anatomy.leadfield.voxel_pos(h.cfg.source.brain_noise_idx,3),'k.');
function menu_head_model_CallBack(varargin)
global h

h.listbox_chans.Value = 1; % resetting when swithcing between MEG and EEG sensors beccause of different # sensors

if h.menu_ori_normal.Value == 2 && h.menu_head_model.Value == 1
    h.menu_ori_normal.Value = 1; % Volume Headmodel can only have random orientations
    for v=1:3; h.edit_source_ori(v).ForegroundColor = h.src_clr(v,:); h.edit_source_ori(v).Enable = 'on'; end
end

menu_sens_type_CallBack();

if h.menu_head_model.Value == 1     % Volume
    if h.menu_sens_type.Value == 1 && ~isempty(h.anatomy.leadfield_meg_vol) && ~isempty(h.anatomy.headmodel_meg_vol)  % MEG
        h.anatomy.leadfield = h.anatomy.leadfield_meg_vol;
        h.anatomy.headmodel = h.anatomy.headmodel_meg_vol;
    elseif h.menu_sens_type.Value == 2 && ~isempty(h.anatomy.leadfield_eeg_vol) && ~isempty(h.anatomy.headmodel_eeg_vol) % EEG
        h.anatomy.leadfield = h.anatomy.leadfield_eeg_vol;
        h.anatomy.headmodel = h.anatomy.headmodel_eeg_vol;
    elseif h.menu_sens_type.Value == 3  % M/EEG
        h.anatomy.leadfield = [];  % will need to change this in the future for adding M/EEG
        h.anatomy.headmodel = [];
    end
    %     h.menu_ori_normal.Enable = 'inactive';  h.menu_ori_normal.ForegroundColor=[1 1 1]*0.5;   h.menu_ori_normal.Value = 1;    % Random orientations locked
elseif h.menu_head_model.Value == 2 && ...
        ( ~isempty(h.anatomy.leadfield_meg_cortex) && ~isempty(h.anatomy.headmodel_meg_cortex) || ... % Cortical Surface
        ~isempty(h.anatomy.leadfield_eeg_cortex) && ~isempty(h.anatomy.headmodel_eeg_cortex) )
    if h.menu_sens_type.Value == 1  % MEG
        h.anatomy.leadfield = h.anatomy.leadfield_meg_cortex;
        h.anatomy.headmodel = h.anatomy.headmodel_meg_cortex;
    elseif h.menu_sens_type.Value == 2  % EEG
        h.anatomy.leadfield = h.anatomy.leadfield_eeg_cortex;
        h.anatomy.headmodel = h.anatomy.headmodel_eeg_cortex;
        
    elseif h.menu_sens_type.Value == 3  % M/EEG
        h.anatomy.leadfield = [];  % will need to change this in the future for adding M/EEG
        h.anatomy.headmodel = [];
    end
    if ~isfield(h.anatomy.leadfield,'ori')
        h.anatomy.leadfield.ori = normals(h.anatomy.leadfield.voxel_pos, h.anatomy.mesh_volumes(4).tri);
    end
    %     h.menu_ori_normal.Enable = 'inactive';  h.menu_ori_normal.ForegroundColor=[1 1 1]*.5; h.menu_ori_normal.Value = 2;    % Random orientations locked
end
menu_sens_type_CallBack();
try h.leadfield_txt.String = sprintf('%.f sensors x  %.f orientations x %.f dipoles',size(h.anatomy.leadfield.H)); catch; h.leadfield_txt.String = sprintf('No Lead Fields'); end


if  ~isempty(h.anatomy.leadfield)
    h.listbox_chans.String = h.anatomy.leadfield.label;
else
    w = warndlg(sprintf('\nPlease Create Lead Fields\n'));
    w.Position(3)=350; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left';
    if h.menu_sens_type.Value==1    % MEG
        h.listbox_chans.String = h.anatomy.sens_meg.label; h.listbox_chans.Value=1;
    elseif h.menu_sens_type.Value==2    % EEG
        h.listbox_chans.String = h.anatomy.sens_eeg.label; h.listbox_chans.Value=1;
    end
end

edit_source_CallBack;
update_anatomy_fields;
function menu_sens_type_CallBack(varargin)
global h
if h.menu_sens_type.Value == 1
    if ~isempty(h.anatomy.leadfield)    % MEG
        h.sensors_txt.String = sprintf('%.f sensors (%s) with %.f used for leadfield',size(h.anatomy.sens_meg.chanpos,1),h.anatomy.sens_meg.type,size(h.anatomy.leadfield.H,1));
        h.edit_leadfield_gain.String = '1e6';
        % find sensor labels of Leadfields that were used from sens_eeg
        cfg.sens_idx=[];
        for v=1:length(h.anatomy.sens_meg.label)
            if sum(strcmp(h.anatomy.sens_meg.label(v),h.anatomy.leadfield.label))==1
                cfg.sens_idx = [cfg.sens_idx v];
            end
        end
    else
        cfg.sens_idx = find ( startsWith(h.anatomy.sens_meg.label,'M')==1) ;
        h.sensors_txt.String = sprintf('%.f sensors (%s) -- No Lead Field',size(h.anatomy.sens_meg.chanpos,1),h.anatomy.sens_meg.type);
    end
    h.edit_yscale.String = '-500 500'; h.edit_yscale_txt.String = 'Y Scale (fT)';
    cfg.sens_type = 'meg'; h.anatomy.sens = bs_select_meeg_sensors(cfg,h.anatomy.sens_meg);    % selecting only channels used in Leadfields for forward projecting data
    h.listbox_chans.String = h.anatomy.sens.label;
    h.listbox_chans.ForegroundColor = h.sens_clr;
elseif h.menu_sens_type.Value == 2     % EEG
    if ~isempty(h.anatomy.leadfield)   % MEG
        h.sensors_txt.String = sprintf('%.f sensors (%s) with %.f used for leadfield',size(h.anatomy.sens_eeg.chanpos,1),h.anatomy.sens_eeg.type,size(h.anatomy.leadfield.H,1));
        h.edit_yscale.String = '-100 100'; h.edit_yscale_txt.String = ['Y Scale (' char(181) 'V):'];
        
        % find sensor labels of Leadfields that were used from sens_eeg
        cfg.sens_idx=[];
        for v=1:length(h.anatomy.sens_eeg.label)
            if sum(strcmp(h.anatomy.sens_eeg.label(v),h.anatomy.leadfield.label))==1
                cfg.sens_idx = [cfg.sens_idx v];
            end
        end
    else
        cfg.sens_idx = find( strcmpi(h.anatomy.sens_eeg.chantype,'eeg')==1 );
        h.sensors_txt.String = sprintf('%.f sensors (%s) -- No Lead Field',size(h.anatomy.sens_eeg.chanpos,1),h.anatomy.sens_eeg.type);
    end
    
    h.edit_leadfield_gain.String = '1e-2';
    cfg.sens_type = 'eeg'; h.anatomy.sens = bs_select_meeg_sensors(cfg,h.anatomy.sens_eeg);    % selecting only channels used in Leadfields for forward projecting data
    h.listbox_chans.String = h.anatomy.sens.label;
    h.listbox_chans.ForegroundColor = h.chan_clr;
elseif h.menu_sens_type.Value == 3      % MEEG  = combined MEG and EEG
    h.anatomy.sens = [];
end
update_scale_texts();
function menu_ori_normal_CallBack(varargin)
global h
if h.menu_ori_normal.Value == 1 && (h.menu_head_model.Value == 1 || h.menu_head_model.Value == 2 ) % Random
    for v=1:3; h.edit_source_ori(v).ForegroundColor = h.src_clr(v,:); h.edit_source_ori(v).Enable = 'on'; end
elseif h.menu_ori_normal.Value == 2 && h.menu_head_model.Value == 2  % Cortically Constrained
    if ~isfield(h.anatomy.leadfield,'ori')
        h.anatomy.leadfield.ori = normals(h.anatomy.leadfield.voxel_pos, h.anatomy.mesh_volumes(4).tri);
    end
    for v=1:3; h.edit_source_ori(v).ForegroundColor = [1 1 1]*.5; h.edit_source_ori(v).Enable = 'inactive'; end
elseif h.menu_ori_normal.Value == 2 && h.menu_head_model.Value == 1
    h.menu_ori_normal.Value = 1; % Volume Headmodel can only have random orientations
    for v=1:3; h.edit_source_ori(v).ForegroundColor = h.src_clr(v,:); h.edit_source_ori(v).Enable = 'on'; end
end
update_source_data;
function edit_source_CallBack(varargin)
update_source_data([],[]); plot_3D_mri('on');
function plot_3D_mri(varargin)
global h
try
    src.Tag = ''; update_cfg(src,[]); % updating cfg based on targets set on screen tab's for ERP, PLV, PLI
catch me
end
update_source_cfg;
% h.axes_anatomy.clo;
% h.axes_anatomy.reset;
h.axes_anatomy.NextPlot = 'add'; axes(h.axes_anatomy); axis tight; 
if ~isempty(h.anatomy.mesh_volumes)
    opt.axes_h=h.axes_anatomy;
    opt.source_clr=h.src_clr;
    opt.vol_nums=1:length(h.anatomy.mesh_volumes);
    
    if length(h.anatomy.mesh_volumes)>=4 && h.menu_head_model.Value==2
        vol.FaceAlpha=h.slider_transparency_brain.Value;   vol.FaceColor=h.brain_clr; vol.EdgeAlpha=0; vol.EdgeColor='none'; opt.vol_nums=1;
        vol.tri = h.anatomy.mesh_volumes(4).tri; vol.pos = h.anatomy.mesh_volumes(4).pos;   % cortical surface
        try
            [h.brain_plot_patch,h.source_locs_patch] = bl_plot_source_locs(vol,h.cfg.source.vx_locs,opt);
        catch
            [h.brain_plot_patch,h.source_locs_patch] = bl_plot_source_locs(vol,h.cfg.study.source_locs_mm,opt);
        end
    else
        vol.FaceAlpha=h.slider_transparency_brain.Value;   vol.FaceColor=[1 1 1]*.6; vol.EdgeAlpha=0; vol.EdgeColor='none';
        vol.FaceColor=h.brain_clr; % [.6 .8 1];
        vol.tri =h.anatomy.mesh_volumes(3).tri; vol.pos = h.anatomy.mesh_volumes(3).pos;
        opt.vol_nums=1;
        [h.brain_plot_patch,h.source_locs_patch] = bl_plot_source_locs(vol,h.cfg.source.vx_locs,opt);
    end
    
    % plot scalp
    if strcmp(varargin,'off')
    else
        vol.FaceAlpha=h.slider_transparency_scalp.Value; vol.FaceColor=h.scalp_clr; vol.tri = h.anatomy.mesh_volumes(1).tri; vol.pos = h.anatomy.mesh_volumes(1).pos;
        opt.vol_nums=1;  h.scalp_plot_patch = bl_plot_mesh(vol,opt);
        vol = h.anatomy.mesh_volumes(1);
    end
    
    if h.menu_sens_type.Value == 1      %MEG
        ft_plot_sens(h.anatomy.sens_meg,'label','on','coilshape','circle','coilsize',12,'fontsize',8,'fontcolor',h.sens_clr,'facecolor',h.sens_clr,'facealpha',0,'edgecolor',h.sens_clr,'chantype', h.anatomy.sens_meg.chantype{1});
    elseif h.menu_sens_type.Value == 2      %EEG
        ft_plot_sens(h.anatomy.sens_eeg,'label','on','elecshape','circle','elecsize',8,'fontsize',8,'fontcolor',h.chan_clr,'facecolor',h.chan_clr,'facealpha',0,'edgecolor',h.chan_clr,'chantype', h.anatomy.sens_eeg.chantype{1});
    elseif h.menu_sens_type.Value == 3      %MEEG
        ft_plot_sens(h.anatomy.sens_meg,'label','on','coilshape','circle','coilsize',12,'fontsize',8,'fontcolor',h.sens_clr,'facecolor',h.sens_clr,'facealpha',0,'edgecolor',h.sens_clr,'chantype', h.anatomy.sens_meg.chantype{1});
        ft_plot_sens(h.anatomy.sens_eeg,'label','on','elecshape','circle','elecsize',8,'fontsize',8,'fontcolor',h.chan_clr,'facecolor',h.chan_clr,'facealpha',0,'edgecolor',h.chan_clr,'chantype', h.anatomy.sens_eeg.chantype{1});
    end
    
    
    %     legend({'Source 1' 'Source 2' 'Source 3'},'Location','northwest'); axis tight;
    % plotting orientations
    try
        for v=1:3
            amp_gain = (max(max(abs([h.cfg.source.sig_amp_perc(v,:),h.cfg.source.prepost_amp_perc(v,:)])))/100)*str2num(h.edit_source_amp(v).String)/2;
            plot3(h.axes_anatomy, [h.cfg.source.vx_locs(v,1) h.cfg.source.vx_locs(v,1)+(h.cfg.source.vx_ori(v,1)*amp_gain) ],...
                [h.cfg.source.vx_locs(v,2) h.cfg.source.vx_locs(v,2)+(h.cfg.source.vx_ori(v,2)*amp_gain) ],...
                [h.cfg.source.vx_locs(v,3) h.cfg.source.vx_locs(v,3)+(h.cfg.source.vx_ori(v,3)*amp_gain) ],...
                'Color',h.src_clr(v,:),'linewidth',2)
        end
    catch
        if h.start_flag == 0
            fprintf('ERROR! Source data have not been simulated yet\n');
        end
    end
    
else
    msgbox('Please Load in HeadModel Volume');
end
get_patch_axes_anatomy(); % get objects for Topo, brain, and scalp for slider transparency
chan_text_OnOff_Callback;
toggle_light_OnOff_axes_anatomy;
function chan_text_OnOff_Callback(varargin)
global h
try
    if h.toggle_chan_text_OnOff.Value == 1   % labels On
        h.toggle_chan_text_OnOff.ForegroundColor = [0 .6 0];    h.toggle_chan_text_OnOff.String = 'Labels On';
        set(findobj(h.axes_anatomy.Children,'Type','Text'),'visible','on')
    elseif h.toggle_chan_text_OnOff.Value == 0   % labels Off
        h.toggle_chan_text_OnOff.ForegroundColor = [1 0 0];    h.toggle_chan_text_OnOff.String = 'Labels Off';
        set(findobj(h.axes_anatomy.Children,'Type','Text'),'visible','off')
    end
    
    if h.toggle_sens_OnOff.Value == 1   % labels On
        h.toggle_sens_OnOff.ForegroundColor = [0 .6 0];    h.toggle_sens_OnOff.String = 'Sensors On';
        h.sens_plot_patch.Visible = 'on';
    elseif h.toggle_sens_OnOff.Value == 0   % labels Off
        h.toggle_sens_OnOff.ForegroundColor = [1 0 0];    h.toggle_sens_OnOff.String = 'Sensors Off';
        h.sens_plot_patch.Visible = 'off';
    end
catch
end
function update_source_data(varargin)
global h
% try

% update source info based on Anatomy panel values
if ~isempty(h.anatomy.leadfield)
    for v=1:3
        h.cfg.source.vx_locs(v,:)                               = str2num(h.edit_source_locs(v).String);  % source locations X, Y, Z (mm)
        h.cfg.source.vx_idx(v)                                  = find_nearest_voxel(h.cfg.source.vx_locs(v,:),h.anatomy.leadfield.voxel_pos);    % source's voxel index from leadfield positions
        h.cfg.source.vx_locs(v,:)                               = h.anatomy.leadfield.voxel_pos(h.cfg.source.vx_idx(v),:);
        
        if h.menu_ori_normal.Value == 1     % Random
            az_el = deg2rad(str2num(h.edit_source_ori(v).String)); [x,y,z] = sph2cart(az_el(1),az_el(2),1);
            h.cfg.source.vx_ori(v,:)                                = [x y z];  % source orientations (X, Y, Z)
        elseif h.menu_ori_normal.Value == 2 && h.menu_head_model.Value == 2  % Cortically Constrained
            %             h.cfg.source.vx_ori(v,:) = h.anatomy.leadfield_eeg_cortex.ori(h.cfg.source.vx_idx(v),:);
            h.cfg.source.vx_ori(v,:) = h.anatomy.leadfield.ori(h.cfg.source.vx_idx(v),:);
            
            [az,el] = cart2sph(h.cfg.source.vx_ori(v,1),h.cfg.source.vx_ori(v,2),h.cfg.source.vx_ori(v,3));
            h.edit_source_ori(v).String = sprintf('%.f %.f',rad2deg(az),rad2deg(el));
        end
        h.cfg.source.vx_amp(v)= str2num(h.edit_source_amp(v).String);   % nAmps
    end
    h.cfg.study.source_locs_mm = h.anatomy.leadfield.voxel_pos(h.cfg.source.vx_idx,:);
    h.cfg.study.source_locs = ft_warp_apply(inv(h.anatomy.mri.transform),h.cfg.study.source_locs_mm);
    
    % [v_idx]=find_nearest_voxel(h.cfg.source.vx_locs,h.anatomy.leadfield.voxel_pos);
    for v=1:3; h.edit_source_locs(v).String = sprintf('%.f %.f %.f',h.anatomy.leadfield.voxel_pos(h.cfg.source.vx_idx(v),:)); end
else
    warning('No Lead fields found');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% Monte Carlo Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_monte_source_modeling_CallBack(varargin)
global h
if h.radio_monte_source_modeling.Value==1
    h.panel_monte_inverse_model.Visible = 'on';
elseif h.radio_monte_source_modeling.Value==0
    h.panel_monte_inverse_model.Visible = 'off';
end
function update_monte_carlo(varargin)
global h

if h.menu_monte_synthetic_real_data.Value >= 2  % Real Sensor Data
    if exist(h.real_datadir_txt.String,'dir')
        xdir = dir(fullfile(h.real_datadir_txt.String,'*.mat'));
        h.monte_real_sensor_files = {xdir.name};
    else
        h.monte_real_sensor_files = [];
    end
    h.edit_monte_num_sims.String = length(h.monte_real_sensor_files);
end

%% initializing variables
h.cfg.study.study_name = h.edit_study_name.String;

% Checking that sensors exist
if strcmp(h.menu_sens_type.Enable,'inactive')   % means that one of the sensor types does not exist
    if any(h.listbox_monte_sens_type.Value ~= h.menu_sens_type.Value)
        h.listbox_monte_sens_type.Value = h.menu_sens_type.Value;
        w = warndlg(sprintf('\nOnly one sensor type exists.\n\nSetting Monte-Carlo forward model\nSensors to: %s\n', h.menu_sens_type.String{h.menu_sens_type.Value}),'Only one Sensor Type!');
        w.Position(3)=250; htext = findobj(w, 'Type', 'Text'); htext.FontSize = 11; htext.HorizontalAlignment = 'left'; % setting fontsize to being readable
    end
else
end
%% Ranges
try
    h.monte_params.num_sims = h.edit_monte_num_sims.String;
    h.monte_params.source_amp_range=[];
    h.monte_params.source_loc_range_X=[];  h.monte_params.source_loc_range_Y=[];  h.monte_params.source_loc_range_Z=[];
    h.monte_params.source_ori_range_Az=[]; h.monte_params.source_ori_range_El=[];
    for v=1:3
        % Amp
        h.monte_params.source_amp_range(v,:) = str2num(h.edit_monte_source_amp_range(v).String);
        % Loc
        h.monte_params.source_loc_range_X(v,:) = str2num(h.edit_monte_source_loc_range_X(v).String);
        h.monte_params.source_loc_range_Y(v,:) = str2num(h.edit_monte_source_loc_range_Y(v).String);
        h.monte_params.source_loc_range_Z(v,:) = str2num(h.edit_monte_source_loc_range_Z(v).String);
        % Ori
        h.monte_params.source_ori_range_Az(v,:) = str2num(h.edit_monte_source_ori_range_Az(v).String);
        h.monte_params.source_ori_range_El(v,:) = str2num(h.edit_monte_source_ori_range_El(v).String);
        
    end
catch
    msgbox(sprintf('Make sure that Ranges for all three sources have the same number of steps\nwithin each category:\n          Amplitude\n          Location X\n          Location Y\n          Location Z\n          Orientation Az\n          Orientation El'))
end
% SNR
%% Std Devs
% Amp
h.monte_params.source_amp_StdDev = str2num(h.edit_monte_source_amp_StdDev.String);
% Locs
h.monte_params.source_loc_StdDev_X = str2num(h.edit_monte_source_loc_StdDev_X.String);
h.monte_params.source_loc_StdDev_Y = str2num(h.edit_monte_source_loc_StdDev_Y.String);
h.monte_params.source_loc_StdDev_Z = str2num(h.edit_monte_source_loc_StdDev_Z.String);
% Ori
h.monte_params.source_ori_StdDev_Az = str2num(h.edit_monte_source_ori_StdDev_Az.String);
h.monte_params.source_ori_StdDev_El = str2num(h.edit_monte_source_ori_StdDev_El.String);
% SNR
h.monte_params.SNR_StdDev = str2num(h.edit_monte_SNR_StdDev.String);
h.monte_params.plv_StdDev = str2num(h.edit_monte_plv_StdDev.String);
%% other params
h.monte_params.SNR_range = str2num(h.edit_monte_SNR_range.String); h.monte_params.SNR_StdDev = str2num(h.edit_monte_SNR_StdDev.String);
h.monte_params.plv_range = str2num(h.edit_monte_plv_range.String); h.monte_params.plv_StdDev = str2num(h.edit_monte_plv_StdDev.String);

h.monte_params.inv_soln = h.listbox_monte_inv_soln.Value;
h.monte_params.FC_analysis = h.listbox_monte_inv_analyses.Value;
h.monte_params.save_FT_data = h.radio_monte_save_FT_data.Value;
h.monte_params.sens_type =h.listbox_monte_sens_type.Value;
function menu_monte_synthetic_real_data_Callback(varargin)
global h
h.btn_set_real_datadir.Visible = 'off';
h.real_datadir_txt.Visible = 'off';

if h.menu_monte_synthetic_real_data.Value == 4  % Real Sensor Noise
    h.btn_set_real_datadir.Visible = 'on';
    h.real_datadir_txt.Visible = 'on';
end
if h.menu_monte_synthetic_real_source.Value == 2  % Real Source Data
    h.btn_set_real_datadir.Visible = 'on';
    h.real_datadir_txt.Visible = 'on';
end
function btn_set_real_datadir(varargin)
global h

h.real_datadir = uigetdir(h.real_datadir,'Set Real Data Directory');
h.real_datadir_txt.String = h.real_datadir;
function run_monte_carlo(varargin)
global h

update_monte_carlo; % just in case user changed something before clicking on run
% h.monte_real_file_idx = 0; % index of files within the "Set Real Dir" for running "Real Sensor" Monte-Carlo simulations
if exist(h.real_datadir_txt.String,'dir')
    xdir = dir(fullfile(h.real_datadir_txt.String,'*.mat'));
    h.monte_real_sensor_files = {xdir.name};
else
    h.monte_real_sensor_files = [];
end

h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Monte Carlo Simulation\n\n     Running');
h.monte_carlo_flag = 1;    % turning on flag for when Monte Carlo Sims are running
src.Tag = '';
update_cfg(src,[]); % updating cfg based on targets set on screen tab's for ERP, PLV, PLI

%% get original Source Amplitudes, Locs, and Ori
clear init_source_amps init_source_locs init_source_ori
for v = 1:3
    init_source_amps(v) = str2num(h.edit_source_amp(v).String);
    init_source_locs(v,:) = str2num(h.edit_source_locs(v).String);
    init_source_ori(v,:) = str2num(h.edit_source_ori(v).String);
end
snr_org = h.edit_sens_SNR.String;
%% checking loc_range sizes are equal across X, Y, Z
dimsX = size( h.monte_params.source_loc_range_X );
dimsY = size( h.monte_params.source_loc_range_Y );
dimsZ = size( h.monte_params.source_loc_range_Z );
run_flag=0;
if isequal(dimsX,dimsY,dimsZ)
    run_flag = 1;
else
    w=msgbox(sprintf('Make sure that all the Locations Change Ranges\nhave the same dimensions'));
    w.Position(3)=300; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left';
end
%% checking ori_range sizes are equal across X, Y, Z
dims1 = size( h.monte_params.source_ori_range_Az );
dims2 = size( h.monte_params.source_ori_range_El);
run_flag=0;
if isequal(dims1,dims2)
    run_flag = 1;
else
    w=msgbox(sprintf('Make sure that all the Orientation Ranges\nhave the same dimensions'));
    w.Position(3)=300; htext = findobj(w, 'Type', 'Text'); htext.FontSize = h.font_size_warndlg; htext.HorizontalAlignment = 'left';
end
%% saving anatomy so it isn't stored redundantly
xsname = sprintf('%s_Anatomy',h.cfg.study.study_name);
sname = fullfile(h.data_dir,[xsname '.mat']);
if exist(sname,'file')
    fprintf('Anatomy File Already Exists: %s\n',sname);
else
    fprintf('Saving Anatomy File: %s\n',sname);
    anatomy = h.anatomy;
    save(sname,'-struct','anatomy');
end

%% Run Monte Loop
hw = waitbar(0,'Running Monte-Carlo Simulations');
% h.monte_params=[];
h.monte_params.sim_run_num = str2num(h.edit_monte_run_num.String)-1; % iterating through monte carlo simulation

if h.menu_monte_synthetic_real_data.Value == 1  % 'Synthetic Data' - simulate data from user-defined singals+prepost
    num_runs = str2num(h.edit_monte_num_sims.String);
elseif h.menu_monte_synthetic_real_data.Value >= 2  % 'Real Sensors + Sources' - load real sensors data and randomize phase + load in source waveform data
    num_runs = length(h.monte_real_sensor_files); h.edit_monte_num_sims.String = num_runs;
end

num_trials = str2num(h.edit_monte_num_trials.String);
num_reruns = str2num(h.edit_monte_num_reruns.String);
nn=0; % number of simulations completed as iterating through all possibile parameters

for nr = 1:num_reruns
%     h.monte_params.sim_run_num = h.monte_params.sim_run_num+1;
    for r = 1:num_runs     % Looping through number of runs
    h.monte_params.sim_run_num =  h.monte_params.sim_run_num + 1;

        for nt = 1:length(num_trials) % numbr of trials
            h.edit_num_trials.String = num2str(num_trials(nt));
            
            for p = 1:length(h.monte_params.plv_range)  % this level needs to run "Sim Source Data" btn and simulate new source waveforms
                fprintf('Sim Run #%.f: \n',h.monte_params.sim_run_num);
                %     fprintf('                    Source 1          Source 2           Source 3\n');
                %             fprintf('Amplitudes(nA):      %.1f              %.1f               %.1f               \n',source_amps);
                
                h.monte_params.cfg = h.cfg; % this will be used for looping through the source parameters on the Monte Carlo tab
                
                %% Update PLV/PLI
                plv_change = h.monte_params.plv_range(p) + (randn*h.monte_params.plv_StdDev);   % making it consistent change for all PLV/PLI so that find_plv function doesn't fail.
                % PLV
                h.monte_params.cfg.source.sig_PLV_targets(h.monte_params.cfg.source.sig_PLV_targets~=0) = h.monte_params.cfg.source.sig_PLV_targets(h.monte_params.cfg.source.sig_PLV_targets~=0) + plv_change;
                h.monte_params.cfg.source.prepost_PLV_targets(h.monte_params.cfg.source.prepost_PLV_targets~=0) = h.monte_params.cfg.source.prepost_PLV_targets(h.monte_params.cfg.source.prepost_PLV_targets~=0) + plv_change;
                % PLI
                h.monte_params.cfg.source.sig_PLI_targets(h.monte_params.cfg.source.sig_PLI_targets~=0) = h.monte_params.cfg.source.sig_PLI_targets(h.monte_params.cfg.source.sig_PLI_targets~=0) + plv_change;
                h.monte_params.cfg.source.prepost_PLI_targets(h.monte_params.cfg.source.prepost_PLI_targets~=0) = h.monte_params.cfg.source.prepost_PLI_targets(h.monte_params.cfg.source.prepost_PLI_targets~=0) + plv_change;
                
                % making sure that PLV >0 and PLV<1
                h.monte_params.cfg.source.sig_PLV_targets(h.monte_params.cfg.source.sig_PLV_targets < 0) = 0;    h.monte_params.cfg.source.sig_PLV_targets(h.monte_params.cfg.source.sig_PLV_targets > 1) = 1;
                h.monte_params.cfg.source.prepost_PLV_targets(h.monte_params.cfg.source.prepost_PLV_targets < 0) = 0;    h.monte_params.cfg.source.prepost_PLV_targets(h.monte_params.cfg.source.prepost_PLV_targets > 1) = 1;
                % making sure that PLI >-1 and PLI<1
                h.monte_params.cfg.source.sig_PLV_targets(h.monte_params.cfg.source.sig_PLV_targets < -1) = -1;    h.monte_params.cfg.source.sig_PLV_targets(h.monte_params.cfg.source.sig_PLV_targets > 1) = 1;
                h.monte_params.cfg.source.prepost_PLV_targets(h.monte_params.cfg.source.prepost_PLV_targets < -1) = -1;    h.monte_params.cfg.source.prepost_PLV_targets(h.monte_params.cfg.source.prepost_PLV_targets > 1) = 1;
                
                
                %% Run a new Sim --> Simulate New Source Waveforms
                %% removing original sens data
                if isfield(h.sim_data,'sens_final_org')
                    h.sim_data = rmfield(h.sim_data,'sens_final_org');
                    h.sim_data = rmfield(h.sim_data,'sens_noise_final_org');
                    h.sim_data = rmfield(h.sim_data,'sens_sig_data_org');
                    h.sim_data = rmfield(h.sim_data,'sig_final_org');
                end
                
                %% Source Waveforms: Simulate "Synthetic" or Load "Real" - same source waveforms for all following changes in parameters such as location, orientation, sensor type, amplitudes, SNRs
                if h.menu_monte_synthetic_real_source.Value == 1  % 'Synthetic Source Waveforms' - simulate source waveforms from user-defined singals+prepost
                    [h.sim_data.sig_final,h.sim_data.sig_wav,h.sim_data.prepost_wav,h.sim_data.noise_wav,h.sim_data.cfg,h.sim_data.prepost_win,h.sim_data.sig_win] = SimSignals(h.monte_params.cfg);
                elseif h.menu_monte_synthetic_real_source.Value == 2  % 'Real Sensors + Sources' - load real sensors data and randomize phase + load in source waveform data
                    %             h.monte_real_file_idx =  h.monte_real_file_idx + 1; % index of files within the "Set Real Dir" for running "Real Sensor" Monte-Carlo simulations
                    if exist(h.real_datadir_txt.String,'dir')
                        h.monte_real_dsname = fullfile(h.real_datadir_txt.String,h.monte_real_sensor_files{r});
                        if exist(h.monte_real_dsname,'file')
                            sm_load_source_waves;          % load source data
                        end
                    end
                end
                
                %% %%%%%% LOOPS for Sensor Type, Amplitude, Locations, Orientations, and SNR
                %% LOCATION Change from initial source locations based on random norm distribution for parameters set in Monte Carlo Tab by "Location Range (min:step:max)" and "Std Dev"
                for lc = 1:size(h.monte_params.source_loc_range_X,2) % looping through "Location Range (min:step:max)"
                    % X Location
                    locX_mu    = h.monte_params.source_loc_range_X(:,lc);   % "Location Range (min:step:max)"
                    sigma = h.monte_params.source_loc_StdDev_X;       % "Std Dev"
                    change_source_locs_X = round( (locX_mu + (sigma * randn(3,1)))*10)/10 ;    % rounded to 1st decimal place
                    % Y Location
                    locY_mu    = h.monte_params.source_loc_range_Y(:,lc);   % "Location Range (min:step:max)"
                    sigma = h.monte_params.source_loc_StdDev_Y;       % "Std Dev"
                    change_source_locs_Y = round( (locY_mu + (sigma * randn(3,1)))*10)/10 ;    % rounded to 1st decimal place
                    % Z Location
                    locZ_mu    = h.monte_params.source_loc_range_Z(:,lc);   % "Location Range (min:step:max)"
                    sigma = h.monte_params.source_loc_StdDev_Z;       % "Std Dev"
                    change_source_locs_Z = round( (locZ_mu + (sigma * randn(3,1)))*10)/10 ;    % rounded to 1st decimal place
                    source_locs = init_source_locs + [change_source_locs_X, change_source_locs_Y, change_source_locs_Z];
                    
                    % setting locs before simulating M/EEG
                    for v=1:3; h.edit_source_locs(v).String = num2str(source_locs(v,:)); end
                    
                    %% ORIENTATION Change from initial source locations based on random norm distribution for parameters set in Monte Carlo Tab by "Orientation Range (min:step:max)" and "Std Dev"
                    for m = 1:size(h.monte_params.source_ori_range_Az,2) % looping through "Location Range (min:step:max)"
                        
                        % Orientation Az
                        oriAz_mu    = h.monte_params.source_ori_range_Az(:,m);   % "Orientation Range (min:step:max)"
                        sigma = h.monte_params.source_ori_StdDev_Az;       % "Std Dev"
                        change_source_ori_Az = round( (oriAz_mu + (sigma * randn(3,1)))*10)/10 ;    % rounded to 1st decimal place
                        % Orientation El
                        oriEl_mu    = h.monte_params.source_ori_range_El(:,m);   % "Orientation (min:step:max)"
                        sigma = h.monte_params.source_ori_StdDev_El;       % "Std Dev"
                        change_source_ori_El = round( (oriEl_mu + (sigma * randn(3,1)))*10)/10 ;    % rounded to 1st decimal place
                        
                        source_ori = init_source_ori + [change_source_ori_Az, change_source_ori_El];
                        
                        % setting orientations before simulating M/EEG
                        for v=1:3; h.edit_source_ori(v).String = num2str(source_ori(v,:)); end
                        
                        %% AMPLITUDE Values based on random norm distribution for parameters set in Monte Carlo Tab by "Amplitude Range (min:step:max)" and "Std Dev"
                        for a=1:size(h.monte_params.source_amp_range,2) % looping through "Amplitude Range (min:step:max)"
                            
                            amp_mu    = h.monte_params.source_amp_range(:,a);   % "Amplitude Range (min:step:max)"
                            sigma = h.monte_params.source_amp_StdDev;       % "Std Dev"
                            source_amps = round( (amp_mu + (sigma * randn(3,1)))*10)/10 ;    % rounded to 1st decimal place
                            % setting amplitudes before simulating M/EEG
                            for v=1:3; h.edit_source_amp(v).String = num2str(source_amps(v)); end
                            
                            
                            %% SNR based on random norm distribution for parameters set in Monte Carlo Tab by "Location Range (min:step:max)" and "Std Dev"
                            for s = 1:length(h.monte_params.SNR_range) % looping through "SNR Range (min:step:max)"
                                %% Orientation Az
                                snr_mu    = h.monte_params.SNR_range(s);   % "SNR Range (min:step:max)"
                                sigma = h.monte_params.SNR_StdDev;       % "Std Dev"
                                snr_val = snr_mu + (sigma * randn(1,1));
                                
                                h.edit_sens_SNR.String = sprintf('%.3f',snr_val);
                                
                                h.cfg.study.SNR = str2num(h.edit_sens_SNR.String);
                                h.sim_data.cfg.study.SNR = str2num(h.edit_sens_SNR.String);
                                
                                update_source_cfg;  % updating source configuration information before simulating M/EEG
                                update_source_edit_boxes;
                                %                         plot_3D_mri;
                                %                         drawnow
                                
                                %% SENSOR Type
                                for st = 1:length(h.monte_params.sens_type) % SENSOR Type
                                    %% SENSOR Type
                                    h.menu_sens_type.Value = h.monte_params.sens_type(st);
                                    
                                    edit_source_CallBack([],[]);
                                    
                                    menu_head_model_CallBack();
                                    
                                    %% Simulate Sensor Noise: Simulate "Synthetic" of Load "Real" Sensor Noise
                                    if h.menu_monte_synthetic_real_data.Value < 3  % 'Synthetic Data' - simulate data from user-defined singals+prepost
                                        %% Simulating New Noise based on
                                        sim_sens_noise;
                                    elseif h.menu_monte_synthetic_real_data.Value == 4  % 'Real Sensors + Sources' - load real sensors data and randomize phase + load in source waveform data
                                        %             h.monte_real_file_idx =  h.monte_real_file_idx + 1; % index of files within the "Set Real Dir" for running "Real Sensor" Monte-Carlo simulations
                                        if exist(h.real_datadir_txt.String,'dir')
                                            h.monte_real_dsname = fullfile(h.real_datadir_txt.String,h.monte_real_sensor_files{r});
                                            if exist(h.monte_real_dsname,'file')
                                                sm_load_real_sensor_file([],[],'noise');   % load real sensor data to h.sim_data.sens_noise as generated sensor noise
                                                sm_fft_randomize_phase;     % randomize phase of sensor data
                                            end
                                            
                                        end
                                    elseif h.menu_monte_synthetic_real_data.Value == 5  % 'Real Full Data' - load real sensors data and analyze - no simulation conducted
                                        sm_load_real_sensor_file([],[],'final');   % load real sensor data to h.sim_data.sens_final as generated sensor data for analyses
                                        if isfield(h,'filt_design'); if ~isempty(h.filt_design);sm_filter_data_Callback; end; end    % filtering data
                                    end
                                    
                                    %% simulating M/EEG data
                                    sim_meeg;    % combines the sensor data and the forward-projected source waveforms for respective SNR values
                                    
                                    %% Filter Data
                                    if h.radio_monte_filter_flag.Value == 1
                                        varin.String ='Filter Data';
                                        if isfield(h,'filt_design'); if ~isempty(h.filt_design); sm_filter_data_Callback(varin); end; end    % filtering data
                                    end
                                    
                                    %% running batch source modeling
                                    if ~isempty(h.listbox_monte_inv_soln.Value)
                                        %% clearing all inverse solutions to start fresh
                                        try
                                            h.listbox_inv_solns.Value = 1:length(h.listbox_inv_solns.String);
                                            delete_inv_solns;
                                        catch
                                        end
                                    end
                                    if h.radio_monte_source_modeling.Value==1
                                        for inv = h.listbox_monte_inv_soln.Value
                                            h.current_inv_soln=inv;
                                            h.menu_inv_soln.Value = inv;
                                            menu_inv_soln_txt_Callback;
                                            
                                            if  strcmp(h.menu_inv_soln.String{h.menu_inv_soln.Value},'dics') || ...
                                                    strcmp(h.menu_inv_soln.String{h.menu_inv_soln.Value},'pcc')
                                                sm_ft_freqanalysis; % perform FT timefreq for dics and pcc inverse modeling
                                                h.menu_inv_datatype.Value = 2;  % doing time-freq modeling
                                            end
                                            
                                            run_source_modeling;
                                            h.slider_3D_image_thresh.Value = 0;
                                            h.inv_soln(h.current_inv_soln).soln.plot_thresh =0;
                                            h.slider_3D_image_thresh.Min = 0;
                                            bs_plot_inv_soln;
                                            if ~isempty(h.current_3D_peak_idx)
                                                h.inv_soln(h.current_inv_soln).soln.peak_idx_found = h.current_3D_peak_idx(1:3); % nearest 3 found peaks relative to true source peaks
                                                h.inv_soln(h.current_inv_soln).soln.true_source_idx = h.sim_data.cfg.source.vx_idx; % nearest 3 found peaks relative to true source peaks
                                                
                                                h.inv_soln(h.current_inv_soln).soln.error_locs = h.current_mse_locs;
                                                h.inv_soln(h.current_inv_soln).soln.error_ori = h.current_mse_ori;
                                                h.inv_soln(h.current_inv_soln).soln.error_waves = h.current_mse_evk_waves;
                                            else
                                                h.inv_soln(h.current_inv_soln).soln.peak_idx_found = []; % nearest 3 found peaks relative to true source peaks
                                                h.inv_soln(h.current_inv_soln).soln.true_source_idx = []; % nearest 3 found peaks relative to true source peaks
                                                
                                                h.inv_soln(h.current_inv_soln).soln.error_locs = [];
                                                h.inv_soln(h.current_inv_soln).soln.error_ori = [];
                                                h.inv_soln(h.current_inv_soln).soln.error_waves = [];
                                            end
                                            if h.monte_params.FC_analysis == 2 && ~isempty(h.current_3D_peak_idx) % Perform functional connectivity analyses
                                                sm_downsample_leadfield;
                                                sm_menu_inv_analyses_CallBack;
                                                sm_create_plv_contrasts;
                                                sm_calc_PLV_PLI;
                                                if isfield(h.cfg.source,'TFR_results')
                                                    if isempty(h.cfg.source.TFR_results)
                                                        sm_calc_true_PLV_PLI;
                                                    end
                                                else
                                                    sm_calc_true_PLV_PLI;
                                                end
                                                
                                                %                                         bs_calc_FC;
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.TFR_freqs = h.current_TFR_freqs;
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.coi_wt2 = h.current_coi_wt2;
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.PLV_freqs = h.current_PLV_freqs;
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.pli_lat = h.current_pli_lat;
                                                %
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.avg_wt = h.current_avg_wt;
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.avg_wt_evk = h.current_avg_wt_evk;
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.avg_wt_ind = h.current_avg_wt_ind;
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.plv_data = h.current_plv_data;
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.pli_data = h.current_plv_data;
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.plv_based = h.current_plv_based;
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.pli_based = h.current_pli_based;
                                                %                                         h.inv_soln(h.current_inv_soln).FC_results.dpli_based = h.current_dpli_based;
                                                
                                            elseif h.monte_params.FC_analysis == 2 && isempty(h.current_3D_peak_idx)
                                                h.inv_soln(h.current_inv_soln).TFR_results = [];
                                            end
                                        end
                                    else
                                    end
                                    
                                    if h.menu_monte_synthetic_real_data.Value == 1  % 'Synthetic Data' - simulate data from user-defined singals+prepost
                                        xsname = sprintf('%s_%s_Amps_%.f_%.f_%.f_LocX_%.f_%.f_%.f_LocY_%.f_%.f_%.f_LocZ_%.f_%.f_%.f_Az_%.f_%.f_%.f_El_%.f_%.f_%.f_SNR_%.3f_PLV_%.3f_Trials_%s_Run_%.f.mat',...
                                            h.cfg.study.study_name,h.menu_sens_type.String{h.menu_sens_type.Value},amp_mu,locX_mu,locY_mu,locZ_mu,oriAz_mu,oriEl_mu,snr_mu,h.monte_params.plv_range(p),h.edit_num_trials.String,h.monte_params.sim_run_num);
                                        sname = fullfile(h.data_dir,xsname);
%                                         while exist(sname,'file') % do not overwrite existing file but add run number
%                                             fprintf('File Exist!  %s\n',sname);
%                                             h.monte_params.sim_run_num=h.monte_params.sim_run_num+1;
%                                             xsname = sprintf('%s_%s_Amps_%.f_%.f_%.f_LocX_%.f_%.f_%.f_LocY_%.f_%.f_%.f_LocZ_%.f_%.f_%.f_Az_%.f_%.f_%.f_El_%.f_%.f_%.f_SNR_%.3f_PLV_%.3f_Trials_%s_Run_%.f.mat',...
%                                                 h.cfg.study.study_name,h.menu_sens_type.String{h.menu_sens_type.Value},amp_mu,locX_mu,locY_mu,locZ_mu,oriAz_mu,oriEl_mu,snr_mu,h.monte_params.plv_range(p),h.edit_num_trials.String,h.monte_params.sim_run_num);
%                                             sname = fullfile(h.data_dir,xsname);
%                                             fprintf('New File Name = %s\n',sname);
%                                         end
                                    elseif h.menu_monte_synthetic_real_data.Value >= 2  % 'Real Sensors or Sources' - load real sensors data and randomize phase + load in source waveform data
                                        [~,real_fname,~]= fileparts(h.monte_real_sensor_files{r});
                                        xsname = sprintf('%s_%s_Amps_%.f_%.f_%.f_LocX_%.f_%.f_%.f_LocY_%.f_%.f_%.f_LocZ_%.f_%.f_%.f_Az_%.f_%.f_%.f_El_%.f_%.f_%.f_SNR_%.3f_PLV_%.3f_Trials_%s_Run_%.f.mat',...
                                            real_fname,h.menu_sens_type.String{h.menu_sens_type.Value},amp_mu,locX_mu,locY_mu,locZ_mu,oriAz_mu,oriEl_mu,snr_mu,h.monte_params.plv_range(p),h.edit_num_trials.String,h.monte_params.sim_run_num);
                                        sname = fullfile(h.data_dir,xsname);
%                                         while exist(sname,'file') % do not overwrite existing file but add run number
%                                             fprintf('File Exist!  %s\n',sname);
%                                             h.monte_params.sim_run_num = h.monte_params.sim_run_num+1;
%                                             xsname = sprintf('%s_%s_Amps_%.f_%.f_%.f_LocX_%.f_%.f_%.f_LocY_%.f_%.f_%.f_LocZ_%.f_%.f_%.f_Az_%.f_%.f_%.f_El_%.f_%.f_%.f_SNR_%.3f_PLV_%.3f_Trials_%s_Run_%.f.mat',...
%                                                 real_fname,h.menu_sens_type.String{h.menu_sens_type.Value},amp_mu,locX_mu,locY_mu,locZ_mu,oriAz_mu,oriEl_mu,snr_mu,h.monte_params.plv_range(p),h.edit_num_trials.String,h.monte_params.sim_run_num);
%                                             sname = fullfile(h.data_dir,xsname);
%                                             fprintf('New File Name = %s\n',sname);
%                                         end
                                    end
                                    
                                    fprintf('Saving:    %s\n',sname);
                                    
                                    cfg = h.monte_params.cfg;
                                    monte_params = h.monte_params;
                                    study_name = h.cfg.study.study_name;
                                    sim_data = h.sim_data;
                                    
                                    % converting to single-point precision for storage
                                    sim_data = double2single(sim_data);
                                    
                                    
                                    % Save Field Trip format along with the BS data format
                                    if h.radio_monte_save_FT_data.Value==1
                                        ft_data = convert_bs2ft_data(h.sim_data.sens_final,h.anatomy,h.cfg);
                                        ft_data = h.ft_data;
                                        save(sname,'cfg','sim_data','ft_data','monte_params','study_name');
                                    else
                                        save(sname,'cfg','sim_data','monte_params','study_name');
                                    end
                                    
                                    if ~isempty(h.listbox_monte_inv_soln.Value)
                                        inv_soln = h.inv_soln;
%                                         inv_soln = double2single(inv_soln);
                                        
                                        save(sname,'inv_soln','-append');
                                        %% clearing all inverse solutions to start fresh for next run
                                        h.listbox_inv_solns.Value = 1:length(h.listbox_inv_solns.String);
                                        try delete_inv_solns; catch; end
                                    end
                                    fprintf('Finished Saving\n',sname);
                                    % clearing true source TFR data
                                    if isfield(h.cfg.source,'TFR_results')
                                        h.cfg.source = rmfield(h.cfg.source,'TFR_results');
                                    end
                                    
                                    nn = nn+1; % number of simulations completed
                                    tot_sims = str2num(h.edit_monte_num_sims.String)*length(h.monte_params.plv_range)*size(h.monte_params.source_loc_range_X,2)*size(h.monte_params.source_ori_range_Az,2)*size(h.monte_params.source_amp_range,2)*length(h.monte_params.SNR_range)*length(h.monte_params.sens_type);
                                    wh_prog = nn/tot_sims;
                                    hw = waitbar(wh_prog,hw,'Running Monte-Carlo Simulations');
                                    
                                end
                            end
                        end
                        
                    end
                end
            end
            
            
            
        end     % num_trials
    end
end

close(hw);

%% resetting original Amps, Locs, and Oris
for v = 1:3
    h.edit_source_amp(v).String = num2str(init_source_amps(v));
    h.edit_source_locs(v).String = num2str(init_source_locs(v,:));
    h.edit_source_ori(v).String  = num2str(init_source_ori(v,:));
end
h.edit_sens_SNR.String = snr_org;

update_source_cfg;  % updating source configuration information before simulating M/EEG
update_source_edit_boxes;
plot_3D_mri;
% drawnow

h.monte_carlo_flag = 0;    % turing off flag for when Monte Carlo Sims are running
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');


%% %%%%%%% Misc Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hide_waitfor_panel(varargin)
global h
h.waitfor_panel.Visible = 'off';


%% %%% Projecting Sim and Noise data to Sensors
function sim_meeg(varargin)
% project simulated data
global h

if h.monte_carlo_flag == 1
    h.waitfor_txt.String = sprintf('Simulating M/EEG data. Please wait ...\n'); drawnow;
else
    h.waitfor_panel.Visible='on';
    h.waitfor_txt.String = sprintf('Simulating M/EEG data. Please wait ...\n'); drawnow;
end

%% removing original sens data
if isfield(h.sim_data,'sens_final_org')
    h.sim_data = rmfield(h.sim_data,'sens_final_org');
    h.sim_data = rmfield(h.sim_data,'sens_noise_final_org');
    h.sim_data = rmfield(h.sim_data,'sens_sig_data_org');
end

%% checking if source data exist
h.sim_data.cfg = h.cfg;
if ~isfield(h.sim_data,'sig_final') % source data have not yet been simulated
    src.Tag=''; hobj=''; update_cfg(src,hobj); run_sim(src,hobj);    % simulating source data
end

%% projecting source sig_final to sensors
leadfield = h.anatomy.leadfield;
% leadfield.H = leadfield.H(h.anatomy.sens.good_sensors,:,:);
% h.sim_data.sens_sig_data(:,h.anatomy.sens.good_sensors,:) = project_SimSignals(h.sim_data.sig_final,leadfield,h.cfg.source.vx_idx,h.cfg.source.vx_amp,h.cfg.source.vx_ori);
h.sim_data.sens_sig_data = project_SimSignals(h.sim_data.sig_final,leadfield,h.cfg.source.vx_idx,h.cfg.source.vx_amp,h.cfg.source.vx_ori);

combine_sens_data;

if h.monte_carlo_flag ~= 1
    h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
end
function plot_sens_data(varargin)
global h

if h.menu_sens_type.Value == 1  %MEG
    wave_clr = h.sens_clr;
elseif h.menu_sens_type.Value == 2  %EEG
    wave_clr = h.chan_clr;
elseif h.menu_sens_type.Value == 2  %MEEG
    wave_clr = h.meeg_clr;
end

try
    
    h.edit_yscale_txt.String = ['Y Scale (' char(181) 'V):'];
    
    ss = find(h.cfg.study.lat_sim<=h.cfg.study.base_int(1)); bs1=ss(end);
    ss = find(h.cfg.study.lat_sim<=h.cfg.study.base_int(2)); bs2=ss(end); h.cfg.study.base_samps = bs1:bs2;
    v = h.listbox_chans.Value;
    
    if length(v)>1 %  plotting buterfly plot
        %% plot noise
        data = squeeze(h.sim_data.sens_noise_final(:,v,:)); data = bsxfun(@minus,data,nanmean(data(h.cfg.study.base_samps,:,:)));
        h.axes_sens_noise.NextPlot='replace';
        h.plot_noise2 = plot(h.axes_sens_noise,h.cfg.study.lat_sim,squeeze(nanmean(data,3)),'color',wave_clr,'linewidth',1);
        h.axes_sens_noise.YLim = str2num(h.edit_yscale.String); h.axes_sens_noise.XLim = str2num(h.edit_plot_time_int.String);
        h.axes_sens_noise.Title.String = sprintf('Noise'); h.axes_sens_noise.Title.Color = wave_clr*1;
        
        % plot signal
        data = squeeze(h.sim_data.sens_sig_data(:,v,:)); data = bsxfun(@minus,data,nanmean(data(h.cfg.study.base_samps,:,:)));
        h.axes_sens_signal.NextPlot='replace'; h.plot_signal2 = plot(h.axes_sens_signal,h.cfg.study.lat_sim,squeeze(nanmean(data,3)),'color',wave_clr,'linewidth',1);
        h.axes_sens_signal.YLim = str2num(h.edit_yscale.String); h.axes_sens_signal.XLim = str2num(h.edit_plot_time_int.String);
        h.axes_sens_signal.Title.String = sprintf('Signal');  h.axes_sens_signal.Title.Color = wave_clr*1;
        
        % plot final
        data = squeeze(h.sim_data.sens_final(:,v,:)); data = bsxfun(@minus,data,nanmean(data(h.cfg.study.base_samps,:,:)));
        h.axes_sens_final.NextPlot='replace'; h.plot_final2 = plot(h.axes_sens_final,h.cfg.study.lat_sim,squeeze(nanmean(data,3)),'color',wave_clr,'linewidth',1);
        h.axes_sens_final.YLim = str2num(h.edit_yscale.String); h.axes_sens_final.XLim = str2num(h.edit_plot_time_int.String);
        h.axes_sens_final.Title.String = sprintf('Noise+Signal'); h.axes_sens_final.Title.Color = wave_clr*1;
        
        h.axes_sens_final.XLabel.String = 'Time (sec)'; h.axes_sens_final.YLabel.String = 'Amp (microV)';
        
        if ~isempty(h.topo_lat_samp); plot_topo_line; end
        
    elseif length(v)==1 %  plotting single channel
        
        %% plot noise
        data = squeeze(h.sim_data.sens_noise_final(:,v,:)); data = bsxfun(@minus,data,nanmean(data(h.cfg.study.base_samps,:,:)));
        h.axes_sens_noise.NextPlot='replace'; h.plot_noise1 = plot(h.axes_sens_noise,h.cfg.study.lat_sim,data,'color',h.trial_wave_clr);
        h.axes_sens_noise.NextPlot='add'; h.plot_noise2 = plot(h.axes_sens_noise,h.cfg.study.lat_sim,squeeze(nanmean(data,2)),'color',wave_clr,'linewidth',2);
        h.axes_sens_noise.YLim = str2num(h.edit_yscale.String); h.axes_sens_noise.XLim = str2num(h.edit_plot_time_int.String);
        h.axes_sens_noise.Title.String = sprintf('Noise at %s', h.listbox_chans.String{h.listbox_chans.Value}); h.axes_sens_noise.Title.Color = wave_clr*1;
        
        % plot signal
        data = squeeze(h.sim_data.sens_sig_data(:,v,:)); data = bsxfun(@minus,data,nanmean(data(h.cfg.study.base_samps,:,:)));
        h.axes_sens_signal.NextPlot='replace'; h.plot_signal1 = plot(h.axes_sens_signal,h.cfg.study.lat_sim,data,'color',h.trial_wave_clr);
        h.axes_sens_signal.NextPlot='add'; h.plot_signal2 = plot(h.axes_sens_signal,h.cfg.study.lat_sim,squeeze(nanmean(data,2)),'color',wave_clr,'linewidth',2);
        h.axes_sens_signal.YLim = str2num(h.edit_yscale.String); h.axes_sens_signal.XLim = str2num(h.edit_plot_time_int.String);
        h.axes_sens_signal.Title.String = sprintf('Signal at %s', h.listbox_chans.String{h.listbox_chans.Value});  h.axes_sens_signal.Title.Color = wave_clr*1;
        
        % plot final
        data = squeeze(h.sim_data.sens_final(:,v,:)); data = bsxfun(@minus,data,nanmean(data(h.cfg.study.base_samps,:,:)));
        h.axes_sens_final.NextPlot='replace'; h.plot_final1 = plot(h.axes_sens_final,h.cfg.study.lat_sim,data,'color',h.trial_wave_clr);
        h.axes_sens_final.NextPlot='add'; h.plot_final2 = plot(h.axes_sens_final,h.cfg.study.lat_sim,squeeze(nanmean(data,2)),'color',wave_clr,'linewidth',2);
        h.axes_sens_final.YLim = str2num(h.edit_yscale.String); h.axes_sens_final.XLim = str2num(h.edit_plot_time_int.String);
        h.axes_sens_final.Title.String = sprintf('Noise+Signal at %s', h.listbox_chans.String{h.listbox_chans.Value}); h.axes_sens_final.Title.Color = wave_clr*1;
        
        h.axes_sens_final.XLabel.String = 'Time (sec)'; h.axes_sens_final.YLabel.String = 'Amp (microV)';
        
        if ~isempty(h.topo_lat_samp); plot_topo_line; end
        
    end
catch me
end
function plot_source_data(varargin)
global h
try
    h.edit_yscale_txt.String = ['Y Scale (nA):'];
    
    ss = find(h.cfg.study.lat_sim<=h.cfg.study.base_int(1)); bs1=ss(end);
    ss = find(h.cfg.study.lat_sim<=h.cfg.study.base_int(2)); bs2=ss(end); h.cfg.study.base_samps = bs1:bs2;
    v = h.listbox_sources.Value;
    % plot prepost
    data = squeeze(sum(h.sim_data.prepost_wav(:,v,:,:),4)); data = bsxfun(@minus,data,nanmean(data(h.cfg.study.base_samps,:,:)))*str2num(h.edit_source_amp(v).String);
    h.axes_sens_noise.NextPlot='replace'; h.plot_noise1 = plot(h.axes_sens_noise,h.cfg.study.lat_sim,data,'color',h.trial_wave_clr);
    h.axes_sens_noise.NextPlot='add'; h.plot_noise2 = plot(h.axes_sens_noise,h.cfg.study.lat_sim,squeeze(nanmean(data,2)),'color',h.src_clr(v,:),'linewidth',2);
    h.axes_sens_noise.YLim = str2num(h.edit_yscale.String); h.axes_sens_noise.XLim = str2num(h.edit_plot_time_int.String);
    h.axes_sens_noise.Title.String = sprintf('Prepost at Simulated Source %.f', v); h.axes_sens_noise.Title.Color = h.src_clr(v,:);
    
    
    % plot signal
    data = squeeze(sum(h.sim_data.sig_wav(:,v,:,:),4)); data = bsxfun(@minus,data,nanmean(data(h.cfg.study.base_samps,:,:)))*str2num(h.edit_source_amp(v).String);
    h.axes_sens_signal.NextPlot='replace'; h.plot_signal1 = plot(h.axes_sens_signal,h.cfg.study.lat_sim,data,'color',h.trial_wave_clr);
    h.axes_sens_signal.NextPlot='add'; h.plot_signal2 = plot(h.axes_sens_signal,h.cfg.study.lat_sim,squeeze(nanmean(data,2)),'color',h.src_clr(v,:),'linewidth',2);
    h.axes_sens_signal.YLim = str2num(h.edit_yscale.String); h.axes_sens_signal.XLim = str2num(h.edit_plot_time_int.String);
    h.axes_sens_signal.Title.String = sprintf('Signal at Simulated Source %.f', v); h.axes_sens_signal.Title.Color = h.src_clr(v,:);
    
    % plot final
    data = squeeze(h.sim_data.sig_final(:,v,:)); data = bsxfun(@minus,data,nanmean(data(h.cfg.study.base_samps,:,:)))*str2num(h.edit_source_amp(v).String);
    h.axes_sens_final.NextPlot='replace'; h.plot_final1 = plot(h.axes_sens_final,h.cfg.study.lat_sim,data,'color',h.trial_wave_clr);
    h.axes_sens_final.NextPlot='add'; h.plot_final2 = plot(h.axes_sens_final,h.cfg.study.lat_sim,squeeze(nanmean(data,2)),'color',h.src_clr(v,:),'linewidth',2);
    h.axes_sens_final.YLim = str2num(h.edit_yscale.String); h.axes_sens_final.XLim = str2num(h.edit_plot_time_int.String);
    h.axes_sens_final.Title.String = sprintf('PrePost+Signal at Simulated Source %.f', v); h.axes_sens_final.Title.Color = h.src_clr(v,:);
    
    h.axes_sens_final.XLabel.String = 'Time (sec)'; h.axes_sens_final.YLabel.String = 'Amp (microV)';
    if ~isempty(h.topo_lat_samp); plot_topo_line; end
catch me
    
end
function update_y_scale(varargin)
global h

% updating sensor plots
y_scale = str2num(h.edit_yscale.String);
h.axes_sens_noise.YLim = y_scale;
h.axes_sens_signal.YLim = y_scale;
h.axes_sens_final.YLim = y_scale;
update_scale_texts()
set_topo_caxis;
function update_scale_texts(varargin)
global h
% updating topoplot scale
y_scale = str2num(h.edit_yscale.String);
if h.menu_sens_type.Value==1 % MEG
    h.slider_topo_scale_text_max.String = [num2str(y_scale(2)) ' ' 'fT'];
    h.axes_sens_final.YLabel.String = 'Amp (fT)';
    h.edit_yscale_txt.String = ['Y Scale (fT)'];
elseif h.menu_sens_type.Value==2 % EEG
    h.slider_topo_scale_text_max.String = [num2str(y_scale(2)) ' ' char(181) 'V'];
    h.axes_sens_final.YLabel.String = ['Amp (' char(181) 'V)'];
    h.edit_yscale_txt.String = ['Y Scale (' char(181) 'V)'];
end
function save_ft_data(varargin)       % save
global h

ft_data = convert_bs2ft_data(h.sim_data.sens_final,h.anatomy,h.cfg);
h.ft_data = ft_data;
anatomy = h.anatomy;
cfg = h.cfg;
[fname,fpath]=uiputfile(sprintf('%s_FieldTrip_data_%s.mat',h.edit_study_name.String,h.menu_sens_type.String{h.menu_sens_type.Value}) ,'Save BRANE Lab Format');
h.waitfor_panel.Visible='on'; h.waitfor_txt.String =sprintf('Saving FieldTrip format in %s \n Please wait...',fname); drawnow;

h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Saving %s\n Please wait...',fname);
save(fullfile(fpath,fname),'ft_data','anatomy','cfg');
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
function save_bs_data(varargin)
global h

[fname,fpath]=uiputfile('bs_data.mat','Save BRANE Lab Format');
h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Saving %s\n Please wait...',fname);
sim_data = h.sim_data;
anatomy = h.anatomy;
cfg = h.cfg;
save(fullfile(fpath,fname),'sim_data','anatomy','cfg');
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
function save_bl_data(varargin)     % Save BRANE Lab format
global h
msgbox('Under Construction!');
return
%%
[fname,fpath]=uiputfile('bl_data.mat','Save BRANE Lab Format');
h.waitfor_panel.Visible='on'; h.waitfor_txt.String = sprintf('Saving %s\n Please wait...',fname);
% creating BRANELab format
bl_data.hdr.Fs = h.cfg.study.srate;                                                 %	sampling frequency
bl_data.hdr.nChans = size(h.anatomy.sens_eeg.chanpos,1);                                %	number of channels
bl_data.hdr.nSamples = h.cfg.study.num_samps;                                       %	number of samples per trial
bl_data.hdr.nSamplesPre = abs(round(h.cfg.study.lat_sim(1)*h.cfg.study.srate));     %   number of pre-trigger samples in each trial
bl_data.hdr.nTrials = h.cfg.study.num_trials;                                       %   number of trials
bl_data.hdr.label = h.anatomy.sens_eeg.label;                                           %   cell-array with labels of each channel
bl_data.hdr.orig =  'Simulated M/EEG from SimMEEG.m';                                %   detailled EDF header information

%% Need to make compatible with BRANELab_v4.m but currerntly just making it useful for the SimMEEG format
% bl_data.lat = h.cfg.study.lat_sim;
% bl_data.data = read_biosemi_bdf(bl_data.dsname, bl_data.hdr, 1, bl_data.hdr.nSamples, 1:bl_data.hdr.nChans);
% bl_data.data_raw = bl_data.data;
% bl_data.chan_names=bl_data.hdr.label;
% bl_data.srate=bl_data.hdr.Fs;
% bl_data.num_chans=bl_data.hdr.nChans;
% bl_data.num_samples=bl_data.hdr.nSamples;
% bl_data.num_trials=bl_data.hdr.nTrials;
% bl_data.good_bad_chans=ones(1,bl_data.num_chans);
% bl_data.ref_chans=zeros(1,bl_data.num_chans); %bl_data.ref_chans(1)=1; % re-referencing to channel #1 to start
% bl_data.exclude_chans=zeros(1,bl_data.num_chans);
% bl_data.good_bad_chans(strcmpi('Status',bl_data.chan_names))=0; % excluding 'STATUS' channels
% bl_data.artefact_data=zeros(size(bl_data.data));


save(fullfile(fpath,fname));  % saves all figure data to be opened as a new study
h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
h.main_fig.Name = fname;


%% Plot Topo
function plot_topo_data(varargin)
global h

% getting latency for topo
y_scale = str2num(h.edit_yscale.String);
axes(h.axes_sens_noise); [x,y,btn] = ginput(1);
xs = find(h.cfg.study.lat_sim <= x); xs=xs(end)+2;

if h.menu_sens_type.Value==1  % MEG
    %     pos = h.anatomy.sens_meg.coilpos(:,:);
    %     pos = h.anatomy.sens.coilpos(:,:);
    pos = h.anatomy.sens.chanpos(:,:);
    val = nanmean(h.sim_data.sens_final(xs,:,:),3);
elseif h.menu_sens_type.Value==2  % EEG
    %     pos = h.anatomy.sens_eeg.elecpos(:,:);
    pos = h.anatomy.sens.elecpos(:,:);
    % expanding positions slightly so scalp doesn't affect transparency
    pos = pos*1.08;
    val = nanmean(h.sim_data.sens_final(xs,:,:),3);
end


h.topo_lat_samp = xs;
if isfield(h,'topo_lat_line')
    if isvalid(h.topo_lat_line(1))
        h.topo_lat_line(1).XData = [h.cfg.study.lat_sim(xs) h.cfg.study.lat_sim(xs)];
        h.topo_lat_line(2).XData = [h.cfg.study.lat_sim(xs) h.cfg.study.lat_sim(xs)];
        h.topo_lat_line(3).XData = [h.cfg.study.lat_sim(xs) h.cfg.study.lat_sim(xs)];
        
        h.topo_lat_line(1).YData = y_scale;
        h.topo_lat_line(2).YData = y_scale;
        h.topo_lat_line(3).YData = y_scale;
        
        
        h.topo_lat_line_txt(3).Position(1:2) = [h.cfg.study.lat_sim(xs+2) y_scale(1)*0.85]; h.topo_lat_line_txt(3).String = sprintf('%.3f',h.cfg.study.lat_sim(xs));
    else
        plot_topo_line;
    end
else
    plot_topo_line;
end

plot_3D_mri('on'); % reploting anatomy + sources without scalp 'off' or with scalp 'on'
axes(h.axes_anatomy);
h.axes_anatomy.SortMethod='depth';
ft_plot_topo3d(pos, val,'facealpha',.65);
get_patch_axes_anatomy(); % get objects for Topo, brain, and scalp for slider transparency

set_topo_caxis;

if isfield(h,'colorbar_anatomy')
    if isvalid(h.colorbar_anatomy)
        delete(h.colorbar_anatomy); h.colorbar_anatomy = colorbar(h.axes_anatomy,'Location','westoutside','Position',[h.axes_anatomy.Position(1)-.01 .075 .03 .2]);
    else
        h.colorbar_anatomy = colorbar(h.axes_anatomy,'Location','westoutside','Position',[h.axes_anatomy.Position(1)-.01 .075 .03 .2]);
    end
else
    h.colorbar_anatomy = colorbar(h.axes_anatomy,'Location','westoutside','Position',[h.axes_anatomy.Position(1)-.01 .075 .03 .2]);
end

h.axes_anatomy.SortMethod='depth';
% h.axes_anatomy.SortMethod='childorder';
x=opengl('data');    % cheking which renderer
if x.Software == 1
    xp = findobj(h.axes_anatomy.Children,'type', 'patch');
    x1 = findobj(xp,'FaceAlpha',.65); x1.FaceAlpha = 1;
end

h.axes_anatomy.Colormap = jet(255);
function plot_topo_line
global h;
xs = h.topo_lat_samp;
y_scale = str2num(h.edit_yscale.String);
if isfield(h,'topo_lat_line') && isfield(h,'topo_lat_line_txt')
    if all(isvalid(h.topo_lat_line)) && all(isvalid(h.topo_lat_line_txt))
        for a=1:3
            h.topo_lat_line(a).XData = [h.cfg.study.lat_sim(xs) h.cfg.study.lat_sim(xs)]; h.topo_lat_line(a).XData = y_scale;
        end
        h.topo_lat_line_txt(a).String = sprintf('%.3f',h.cfg.study.lat_sim(xs)); h.topo_lat_line_txt(a).Position(1:2) = [h.cfg.study.lat_sim(xs+2) y_scale(1)*.85];
    else
        delete(h.topo_lat_line); delete(h.topo_lat_line_txt);
        h.axes_sens_noise.NextPlot='add'; h.topo_lat_line(1) = plot(h.axes_sens_noise,[h.cfg.study.lat_sim(xs) h.cfg.study.lat_sim(xs)],y_scale ,'--','color','r');
        h.axes_sens_signal.NextPlot='add'; h.topo_lat_line(2) = plot(h.axes_sens_signal,[h.cfg.study.lat_sim(xs) h.cfg.study.lat_sim(xs)],y_scale ,'--','color','r');
        h.axes_sens_final.NextPlot='add'; h.topo_lat_line(3) = plot(h.axes_sens_final,[h.cfg.study.lat_sim(xs) h.cfg.study.lat_sim(xs)],y_scale ,'--','color','r');
        h.topo_lat_line_txt(3) = text(h.axes_sens_final,h.cfg.study.lat_sim(xs+2),y_scale(1)*.85,sprintf('%.3f',h.cfg.study.lat_sim(xs)),'Color','r');
    end
else
    %     delete(h.topo_lat_line); delete(h.topo_lat_line_txt);
    h.axes_sens_noise.NextPlot='add'; h.topo_lat_line(1) = plot(h.axes_sens_noise,[h.cfg.study.lat_sim(xs) h.cfg.study.lat_sim(xs)],y_scale ,'--','color','r');
    h.axes_sens_signal.NextPlot='add'; h.topo_lat_line(2) = plot(h.axes_sens_signal,[h.cfg.study.lat_sim(xs) h.cfg.study.lat_sim(xs)],y_scale ,'--','color','r');
    h.axes_sens_final.NextPlot='add'; h.topo_lat_line(3) = plot(h.axes_sens_final,[h.cfg.study.lat_sim(xs) h.cfg.study.lat_sim(xs)],y_scale ,'--','color','r');
    h.topo_lat_line_txt(3) = text(h.axes_sens_final,h.cfg.study.lat_sim(xs+2),y_scale(1)*.85,sprintf('%.3f',h.cfg.study.lat_sim(xs)),'Color','r');
end
function set_topo_caxis(varargin)
global h
y_scale = str2num(h.edit_yscale.String);

s_pos = h.slider_topo_scale.Position;
h.slider_topo_scale_text_val.Position(1) = s_pos(1) - h.slider_topo_scale_text_val.Position(3) - .01;
h.slider_topo_scale_text_val.Position(2) = s_pos(2)+(s_pos(4)*h.slider_topo_scale.Value*.7)+(s_pos(4)*.1);

if h.menu_sens_type.Value==1 % MEG
    h.slider_topo_scale_text_val.String = [sprintf('%.f ',h.slider_topo_scale.Value*y_scale(2)) 'fT'];
elseif h.menu_sens_type.Value==2 % EEG
    h.slider_topo_scale_text_val.String = [sprintf('%.f ',h.slider_topo_scale.Value*y_scale(2)) char(181) 'V'];
end


h.topo_scale = y_scale*h.slider_topo_scale.Value;
h.axes_anatomy.CLim = h.topo_scale;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% ANALYSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%% MRI Program for sleecting sources %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% SELECT Source Locations %%%%%%%%%%%%
function select_source_locs(varargin)
global h
% bl_set_source_locs_v2(h.anatomy.mri,h.cfg.study.source_locs,h.anatomy.mesh_cortex);
bl_set_source_locs_v2(h.anatomy.mri,h.cfg.study.source_locs);
function bl_set_source_locs_v2(mri,source_locs)

%%   This program will plot mri and allow for selecting 3 sources to be entered into SimSignals.m software.
%   Source locations can either be set by moving the cursors to a specific location and then clicking on "Set Source 1"
%   or by clicking and dragging any of the 3 sources. Remember to click on "Save to Workspace" once the 3 sources are in
%   the locations you want to simulate.
%
% INPUT
%   mri = Field Trip's structure of loading mri. Note, that the transformation matrix mri.transform should be the same one used for creating the headmodel and leadfields.
%           must contain:
%               mri.dim = dimensions of mri.anatomy
%               mri.hdr = header information read in using Field Trip
%               mri.anatomy = in voxel space coordinates
%               mri.transform = 4x4 tranformation matrix to convert voxel to mm coordinates.
%   source_locs = location in mri.anatomy space coordinates
%   vol = (optional) volume data for source model to plot source_locs_mm in 3D (see bl_plot_source_locs.m for more options)
%
% OUTPUT
%   Click on "Save to Workspace" to load the source locations into matlab's base workspace.
%       source_locs = source locations in voxel coordinates
%       source_locs_mm = source locations in mm after transformation using ft_warp_apply.m (see Field Trip for more information)
%
%% initialize figure
global h2 h
h2=figure(101); h2.Color='w'; h2.Position=[300 100 850 700];
% h2 = h.tab_anatomy; h2.BackgroundColor='w';

if nargin<2 || isempty(source_locs)
    h2.UserData.source_locs=[91 48 78; 35 109 78; 151 109 78]; %ones(3,3); %[];
else
    h2.UserData.source_locs=source_locs;
end
h2.UserData.x_sagittal_slice=h2.UserData.source_locs(1,1); % set starting cursor to first source loc
h2.UserData.y_coronal_slice=h2.UserData.source_locs(1,2);
h2.UserData.z_axial_slice=h2.UserData.source_locs(1,3);

% h2.UserData.fiducial=mri.hdr.fiducial;  % Not implementing plotting of fiducials yet
h2.UserData.transform=mri.transform;
h2.UserData.source_clr=h.src_clr;
h2.UserData.source_size=60;
h2.UserData.anatomy=mri.anatomy;
h2.UserData.dims=size(h2.UserData.anatomy);
h2.UserData.drag_axial_source=0; h2.UserData.drag_sagittal_source=0; h2.UserData.drag_coronal_source=0;

%% Panel 1: MRI information - Upper Left Panel
hp1 = uipanel(h2,'Title','MRI Information','FontSize',12,...
    'BackgroundColor','w','Foregroundcolor','k',...
    'Position',[.01 .51 .47 .47],'Units','normalize');
mri_txt=uicontrol(hp1,'Style','text',...
    'BackgroundColor','w','Foregroundcolor','k',...
    'Position',[.01 .51 .45 .45],'Units','normalize',...
    'String',''); mri_txt.Position=[.01 .81 .9 .15];
try
    h2.UserData.mri_str=sprintf('MRI File: %s\n',mri.hdr.mrifile);
catch
    h2.UserData.mri_str='';
end
mri_txt.String=h2.UserData.mri_str; mri_txt.HorizontalAlignment='left';
%% text output
loc_txt=uicontrol(hp1,'Style','text',...
    'BackgroundColor','w','Foregroundcolor','k',...
    'Position',[.01 .61 .55 .55],'Units','normalize',...
    'String','');loc_txt.Position=[.01 .31 .9 .55];
xyz=[h2.UserData.x_sagittal_slice,h2.UserData.y_coronal_slice,h2.UserData.z_axial_slice];
xyz_w=ft_warp_apply(h2.UserData.transform,[1 1 1]); h2.UserData.current_loc=xyz_w;
h2.UserData.source_locs_mm=ft_warp_apply(h2.UserData.transform,h2.UserData.source_locs); h2.UserData.current_loc=xyz_w;
h2.UserData.loc_str=sprintf('pixel =   %.f   %.f   %.f\n mm =   %.f   %.f   %.f\n\nSource 1 =    %.f   %.f   %.f\nSource 2 =    %.f   %.f   %.f\nSource 3 =    %.f   %.f   %.f\n\nSource 1 =    %.f   %.f   %.f\nSource 2 =    %.f   %.f   %.f\nSource 3 =    %.f   %.f   %.f',...
    xyz,xyz_w,h2.UserData.source_locs',h2.UserData.source_locs_mm');
loc_txt.String=h2.UserData.loc_str; loc_txt.HorizontalAlignment='left';
%% btn for source 1 selection locations
btn1 = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(1,:),'ForegroundColor',[1 1 1]*1,'Style','pushbutton','String','Set Source 1','Position',[.15 .05 .15 .25],'units','normalize','Callback',@select_location); btn1.Position=[.02 .23 .41 .1];
btn2 = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(2,:),'ForegroundColor',[1 1 1]*1,'Style','pushbutton','String','Set Source 2','Position',[.45 .05 .15 .25],'units','normalize','Callback',@select_location); btn2.Position=[.02 .12 .41 .1];
btn3 = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(3,:),'ForegroundColor',[1 1 1]*1,'Style','pushbutton','String','Set Source 3','Position',[.75 .05 .15 .25],'units','normalize','Callback',@select_location); btn3.Position=[.02 .01 .41 .1];
%% btn to go to source 1 selection locations
btn1_loc = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(1,:),'ForegroundColor',[1 1 1]*1,'Style','pushbutton','String','Go to Source 1','Position',[.15 .05 .15 .25],'units','normalize','Callback',@goto_location); btn1_loc.Position=[.52 .23 .41 .1];
btn2_loc = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(2,:),'ForegroundColor',[1 1 1]*1,'Style','pushbutton','String','Go to Source 2','Position',[.45 .05 .15 .25],'units','normalize','Callback',@goto_location); btn2_loc.Position=[.52 .12 .41 .1];
btn3_loc = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(3,:),'ForegroundColor',[1 1 1]*1,'Style','pushbutton','String','Go to Source 3','Position',[.75 .05 .15 .25],'units','normalize','Callback',@goto_location); btn3_loc.Position=[.52 .01 .41 .1];
%% btn save source locs to workspace
btn_save = uicontrol(hp1,'BackgroundColor',[1 .6 0],'Style','pushbutton','String','Save Locations','units','normalize','Callback',@save_location); btn_save.Position=[.652 .88 .3 .1];
%% btn default source locs to workspace
btn_defaults = uicontrol(hp1,'BackgroundColor',[1 1 1]*.8,'Style','pushbutton','String','Default Locs','units','normalize','Callback',@default_location); btn_defaults.Position=[.652 .68 .3 .1];
%% btn for plotting source locs in 3D volume
% btn_3D = uicontrol(hp1,'BackgroundColor',[1 .6 .6],'Style','pushbutton','String','Plot 3D','units','normalize','Callback',@plot_source_locs_3D); btn_3D.Position=[.652 .48 .3 .1];
% if nargin>=3
%     h2.UserData.vol=vol;
%     h2.UserData.opt=[];
% elseif nargin<3
%     btn_3D.Enable='inactive'; btn_3D.BackgroundColor=[1 1 1]*.9;
%     btn_3D.ForegroundColor=[1 1 1]*.6;
%     btn_3D.String='No 3D plot ';
% end

h2.UserData.z_axial_slice=round(h2.UserData.z_axial_slice);
h2.UserData.y_coronal_slice=round(h2.UserData.y_coronal_slice);
h2.UserData.x_sagittal_slice=round(h2.UserData.x_sagittal_slice);

%% Axes 4: Axial Slices - Upper Right Panel
% hp2 = uipanel(h,'Title','Axial','FontSize',12,...
%     'BackgroundColor','w','Foregroundcolor','k',...
%     'Position',[.51 .51 .48 .48],'Units','normalize');
h2.UserData.ax(1)=axes('Position',[.51 .51 .45 .45]); box off; %h2.ax(1).Toolbar.Visible = 'off';
imagesc(h2.UserData.ax(1),squeeze(mri.anatomy(:,:,h2.UserData.z_axial_slice))); axis tight; view(-90,90); axis off;
hold on; p1=plot([0 mri.dim(2)],[h2.UserData.x_sagittal_slice h2.UserData.x_sagittal_slice],'r');
p2=plot([h2.UserData.y_coronal_slice h2.UserData.y_coronal_slice],[0 mri.dim(1)],'r');
vidx=h2.UserData.source_locs(:,3)==h2.UserData.z_axial_slice;
for v=fliplr(1:3)
    if vidx(v)==1
        a4(v)=scatter(h2.UserData.source_locs(v,2),h2.UserData.source_locs(v,1),'o','SizeData',h2.UserData.source_size,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    else    % plot off the image
        a4(v)=scatter(0,0,'o','SizeData',1,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    end
end
axis([0 h2.UserData.dims(2) 0 h2.UserData.dims(1)]);
% % set views
% h2.UserData.edit_axial_view_txt = uicontrol(h2,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
%     'Position',[h2.UserData.ax(1).Position(1)+.1 sum(h2.UserData.ax(1).Position([2 4]))+.001 .05 .03],...
%     'FontSize',10,'HorizontalAlignment','left','String','View');
% h2.UserData.edit_axial_view = uicontrol(h2,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
%     'Position',[sum(h2.UserData.edit_axial_view_txt.Position([1 3]))+.01 h2.UserData.edit_axial_view_txt.Position(2) .1 .03],...
%     'FontSize',10,'HorizontalAlignment','center','String','0 90','Callback',@update_slices);


%% Axes 3:  Sagittal Slices - Lower Left Panel
% hp3 = uipanel(h,'Title','Sagittal','FontSize',12,...
%     'BackgroundColor','w','Foregroundcolor','k',...
%     'Position',[.01 .01 .48 .48],'Units','normalize');
h2.UserData.ax(2)=axes('Position',[.03 .03 .45 .45]); box off; %h2.ax(2).Toolbar.Visible = 'off';
imagesc(h2.UserData.ax(2),squeeze(mri.anatomy(h2.UserData.x_sagittal_slice,:,:))); axis tight; view(-90,90); axis off;
hold on; p1=plot([0 mri.dim(3)],[h2.UserData.y_coronal_slice h2.UserData.y_coronal_slice],'r');
p2=plot([h2.UserData.z_axial_slice h2.UserData.z_axial_slice],[0 mri.dim(2)],'r');
vidx=h2.UserData.source_locs(:,1)==h2.UserData.x_sagittal_slice;
for v=fliplr(1:3)
    if vidx(v)==1
        a3(v)=scatter(h2.UserData.source_locs(v,3),h2.UserData.source_locs(v,2),'o','SizeData',h2.UserData.source_size,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    else    % plot off the image
        a3(v)=scatter(0,0,'o','SizeData',1,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    end
end
axis([0 h2.UserData.dims(3) 0 h2.UserData.dims(2)]);
% % set views
% h2.UserData.edit_sagittal_view_txt = uicontrol(h2,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
%     'Position',[h2.UserData.ax(2).Position(1)+.1 sum(h2.UserData.ax(2).Position([2 4]))+.001 .05 .03],...
%     'FontSize',10,'HorizontalAlignment','left','String','View');
% h2.UserData.edit_sagittal_view = uicontrol(h2,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
%     'Position',[sum(h2.UserData.edit_sagittal_view_txt.Position([1 3]))+.01 h2.UserData.edit_sagittal_view_txt.Position(2) .1 .03],...
%     'FontSize',10,'HorizontalAlignment','center','String','0 90','Callback',@update_slices);


%% Axes 2:  Coronal Slices - Lower Right Panel
% hp4 = uipanel(h,'Title','Coronal','FontSize',12,...
%     'BackgroundColor','w','Foregroundcolor','k',...
%     'Position',[.51 .01 .48 .48],'Units','normalize');
h2.UserData.ax(3)=axes('Position',[.51 .03 .45 .45]); box off; %h2.ax(3).Toolbar.Visible = 'off';
imagesc(h2.UserData.ax(3),squeeze(mri.anatomy(:,h2.UserData.y_coronal_slice,:))); axis tight; view(-90,90); axis off;
hold on; p1=plot([0 mri.dim(3)],[h2.UserData.x_sagittal_slice h2.UserData.x_sagittal_slice],'r');
p2=plot([h2.UserData.z_axial_slice h2.UserData.z_axial_slice],[0 mri.dim(1)],'r');
vidx=h2.UserData.source_locs(:,2)==h2.UserData.y_coronal_slice;
for v=fliplr(1:3)
    if vidx(v)==1
        a2(v)=scatter(h2.UserData.source_locs(v,3),h2.UserData.source_locs(v,1),'o','SizeData',h2.UserData.source_size,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    else    % plot off the image
        a2(v)=scatter(0,0,'o','SizeData',1,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    end
end
axis([0 h2.UserData.dims(1) 0 h2.UserData.dims(3)]);
% % set views
% h2.UserData.edit_coronal_view_txt = uicontrol(h2,'Style','text', 'BackgroundColor',h.UserData.bkg_clr,'Foregroundcolor','k','Units','normalize',...
%     'Position',[h2.UserData.ax(3).Position(1)+.1 sum(h2.UserData.ax(3).Position([2 4]))+.001 .05 .03],...
%     'FontSize',10,'HorizontalAlignment','left','String','View');
% h2.UserData.edit_coronal_view = uicontrol(h2,'BackgroundColor',h.UserData.bkg_clr,'ForegroundColor',[1 1 1]*0,'Style','edit','Units','normalize',...
%     'Position',[sum(h2.UserData.edit_coronal_view_txt.Position([1 3]))+.01 h2.UserData.edit_coronal_view_txt.Position(2) .1 .03],...
%     'FontSize',10,'HorizontalAlignment','center','String','0 90','Callback',@update_slices);

%%
colormap(bone);

% h2.Children(4).Children(end).Tag='axial'; h2.Children(4).Children(end).ButtonDownFcn=@startDragFcn_axial;
% h2.Children(3).Children(end).Tag='sagittal'; h2.Children(3).Children(end).ButtonDownFcn=@startDragFcn_sagittal;
% h2.Children(2).Children(end).Tag='coronal'; h2.Children(2).Children(end).ButtonDownFcn=@startDragFcn_coronal;

h2.UserData.ax(1).Children(end).Tag='axial'; h2.UserData.ax(1).Children(end).ButtonDownFcn=@startDragFcn_axial;
h2.UserData.ax(2).Children(end).Tag='sagittal'; h2.UserData.ax(2).Children(end).ButtonDownFcn=@startDragFcn_sagittal;
h2.UserData.ax(3).Children(end).Tag='coronal'; h2.UserData.ax(3).Children(end).ButtonDownFcn=@startDragFcn_coronal;

h2.UserData.ax(1).Children(end).Tag='axial'; h2.UserData.ax(1).Children(end).ButtonDownFcn=@startDragFcn_axial;
h2.UserData.ax(2).Children(end).Tag='sagittal'; h2.UserData.ax(2).Children(end).ButtonDownFcn=@startDragFcn_sagittal;
h2.UserData.ax(3).Children(end).Tag='coronal'; h2.UserData.ax(3).Children(end).ButtonDownFcn=@startDragFcn_coronal;



% for m=2:4; h2.Children(m).Children(3).PickableParts='all'; end

%% making source locs click & drag
for v=1:3
    h2.UserData.ax(1).Children(v).Tag=sprintf('axial Source %.f',v);
    h2.UserData.ax(1).Children(v).ButtonDownFcn=@startDragFcn_source;
    h2.UserData.ax(2).Children(v).Tag=sprintf('sagittal Source %.f',v);
    h2.UserData.ax(2).Children(v).ButtonDownFcn=@startDragFcn_source;
    h2.UserData.ax(3).Children(v).Tag=sprintf('coronal Source %.f',v);
    h2.UserData.ax(3).Children(v).ButtonDownFcn=@startDragFcn_source;
end

%% set up btn Down/Up funcs
h2.WindowButtonUpFcn=@stopDragFcn; % start ginput to get cursor location on axes
h2.WindowButtonDownFcn=@startDragFcn; % start ginput to get cursor location on axes
%% update image data for selected slice
function update_slices(varargin)
global h2
h2.UserData.z_axial_slice=round(h2.UserData.z_axial_slice);
h2.UserData.y_coronal_slice=round(h2.UserData.y_coronal_slice);
h2.UserData.x_sagittal_slice=round(h2.UserData.x_sagittal_slice);

% view(h2.UserData.ax(1),str2num(h2.UserData.edit_axial_view.String))
% view(h2.UserData.ax(2),str2num(h2.UserData.edit_sagittal_view.String))
% view(h2.UserData.ax(3),str2num(h2.UserData.edit_coronal_view.String))

if h2.UserData.x_sagittal_slice>0 && h2.UserData.x_sagittal_slice<=size(h2.UserData.anatomy,1) && ...
        h2.UserData.y_coronal_slice>0 && h2.UserData.y_coronal_slice<=size(h2.UserData.anatomy,2) && ...
        h2.UserData.z_axial_slice>0 && h2.UserData.z_axial_slice<=size(h2.UserData.anatomy,3)
    
    % update loc_text
    xyz=[h2.UserData.x_sagittal_slice,h2.UserData.y_coronal_slice,h2.UserData.z_axial_slice];
    xyz_w=ft_warp_apply(h2.UserData.transform,xyz); h2.UserData.current_loc=xyz_w;
    h2.UserData.source_locs_mm=ft_warp_apply(h2.UserData.transform,h2.UserData.source_locs); h2.UserData.current_loc=xyz_w;
    h2.UserData.loc_str=sprintf('pixel =   %.f   %.f   %.f\n mm =   %.f   %.f   %.f\n\nSource 1 =    %.f   %.f   %.f   pixel\nSource 2 =    %.f   %.f   %.f   pixel\nSource 3 =    %.f   %.f   %.f   pixel\n\nSource 1 =    %.f   %.f   %.f   mm\nSource 2 =    %.f   %.f   %.f   mm\nSource 3 =    %.f   %.f   %.f   mm',...
        xyz,xyz_w,h2.UserData.source_locs',h2.UserData.source_locs_mm');
    h2.Children(1).Children(end-1).String=h2.UserData.loc_str; loc_txt.HorizontalAlignment='left';
    
    % updating images
    h2.UserData.ax(1).Children(end).CData=squeeze(h2.UserData.anatomy(:,:,h2.UserData.z_axial_slice)); % updating axial slice
    h2.UserData.ax(3).Children(end).CData=squeeze(h2.UserData.anatomy(:,h2.UserData.y_coronal_slice,:)); % updating axial slice
    h2.UserData.ax(2).Children(end).CData=squeeze(h2.UserData.anatomy(h2.UserData.x_sagittal_slice,:,:)); % updating axial slice
    
    % updating cursor lines
    % axial slice
    h2.UserData.ax(1).Children(end-2).XData=[0 size(h2.UserData.anatomy,2)];
    h2.UserData.ax(1).Children(end-2).YData=[h2.UserData.x_sagittal_slice h2.UserData.x_sagittal_slice];
    h2.UserData.ax(1).Children(end-1).XData=[h2.UserData.y_coronal_slice h2.UserData.y_coronal_slice];
    h2.UserData.ax(1).Children(end-1).YData=[0 size(h2.UserData.anatomy,1)];
    % sagital slice
    h2.UserData.ax(2).Children(end-2).XData=[h2.UserData.z_axial_slice h2.UserData.z_axial_slice];
    h2.UserData.ax(2).Children(end-2).YData=[0 size(h2.UserData.anatomy,2)];
    h2.UserData.ax(2).Children(end-1).YData=[h2.UserData.y_coronal_slice h2.UserData.y_coronal_slice];
    h2.UserData.ax(2).Children(end-1).XData=[0 size(h2.UserData.anatomy,3)];
    % coronal slice
    h2.UserData.ax(3).Children(end-2).XData=[h2.UserData.z_axial_slice h2.UserData.z_axial_slice];
    h2.UserData.ax(3).Children(end-2).YData=[0 size(h2.UserData.anatomy,1)];
    h2.UserData.ax(3).Children(end-1).YData=[h2.UserData.x_sagittal_slice h2.UserData.x_sagittal_slice];
    h2.UserData.ax(3).Children(end-1).XData=[0 size(h2.UserData.anatomy,3)];
    
    % updating Source locations when slice is selected
    % find all on axial slice
    vidx=h2.UserData.source_locs(:,3)>=h2.UserData.z_axial_slice-2 & h2.UserData.source_locs(:,3)<=h2.UserData.z_axial_slice+2;
    for v=1:3
        if vidx(v)==1
            h2.UserData.ax(1).Children(v).XData=[h2.UserData.source_locs(v,2) h2.UserData.source_locs(v,2)];
            h2.UserData.ax(1).Children(v).YData=[h2.UserData.source_locs(v,1) h2.UserData.source_locs(v,1)];
            h2.UserData.ax(1).Children(v).SizeData=h2.UserData.source_size;
        else
            h2.UserData.ax(1).Children(v).XData=[0 0];
            h2.UserData.ax(1).Children(v).YData=[0 0];
            h2.UserData.ax(1).Children(v).SizeData=1;
        end
    end
    
    % find all on sagittal slice
    vidx=h2.UserData.source_locs(:,1)>=h2.UserData.x_sagittal_slice-2 & h2.UserData.source_locs(:,1)<=h2.UserData.x_sagittal_slice+2;
    for v=1:3
        if vidx(v)==1
            h2.UserData.ax(2).Children(v).XData=[h2.UserData.source_locs(v,3) h2.UserData.source_locs(v,3)];
            h2.UserData.ax(2).Children(v).YData=[h2.UserData.source_locs(v,2) h2.UserData.source_locs(v,2)];
            h2.UserData.ax(2).Children(v).SizeData=h2.UserData.source_size;
        else
            h2.UserData.ax(2).Children(v).XData=[0 0];
            h2.UserData.ax(2).Children(v).YData=[0 0];
            h2.UserData.ax(2).Children(v).SizeData=1;
        end
    end
    % find all on coronal slice
    vidx=h2.UserData.source_locs(:,2)>=h2.UserData.y_coronal_slice-2 & h2.UserData.source_locs(:,2)<=h2.UserData.y_coronal_slice+2;
    for v=1:3
        if vidx(v)==1
            h2.UserData.ax(3).Children(v).XData=[h2.UserData.source_locs(v,3) h2.UserData.source_locs(v,3)];
            h2.UserData.ax(3).Children(v).YData=[h2.UserData.source_locs(v,1) h2.UserData.source_locs(v,1)];
            h2.UserData.ax(3).Children(v).SizeData=h2.UserData.source_size;
        else
            h2.UserData.ax(3).Children(v).XData=[0 0];
            h2.UserData.ax(3).Children(v).YData=[0 0];
            h2.UserData.ax(3).Children(v).SizeData=1;
        end
    end
    
    
    
else
    h2.UserData.x_sagittal_slice=round(size(h2.UserData.anatomy,1)/2);
    h2.UserData.y_coronal_slice=round(size(h2.UserData.anatomy,2)/2);
    h2.UserData.z_axial_slice=round(size(h2.UserData.anatomy,3)/2);
end
function startDragFcn(varargin)
global h2
if strcmpi(varargin{2}.Source.CurrentAxes.Children(end).Tag,'axial')
    xyz=get(gca,'CurrentPoint');
    h2.UserData.x_sagittal_slice=round(xyz(3));
    h2.UserData.y_coronal_slice=round(xyz(1));
    update_slices;
    startDragFcn_axial;
elseif strcmpi(varargin{2}.Source.CurrentAxes.Children(end).Tag,'sagittal')
    xyz=get(gca,'CurrentPoint');
    h2.UserData.z_axial_slice=round(xyz(1));
    h2.UserData.y_coronal_slice=round(xyz(3));
    update_slices;
    startDragFcn_sagittal;
elseif strcmpi(varargin{2}.Source.CurrentAxes.Children(end).Tag,'coronal')
    xyz=get(gca,'CurrentPoint');
    h2.UserData.z_axial_slice=round(xyz(1));
    h2.UserData.x_sagittal_slice=round(xyz(3));
    update_slices;
    startDragFcn_coronal;
end
function stopDragFcn(varargin)
global h2
h2.WindowButtonMotionFcn='';   % stops dragging func
h2.UserData.drag_axial_source=0; % resetting so that sources are not selected to move during dragging
h2.UserData.drag_sagittal_source=0;
h2.UserData.drag_coronal_source=0;
function startDragFcn_source(varargin)
global h2
if strfind(varargin{2}.Source.Tag,'axial Source')==1
    h2.UserData.drag_axial_source=str2num(varargin{2}.Source.Tag(end));
    startDragFcn_axial;
elseif strfind(varargin{2}.Source.Tag,'sagittal Source')==1
    h2.UserData.drag_sagittal_source=str2num(varargin{2}.Source.Tag(end));
    startDragFcn_sagittal;
elseif strfind(varargin{2}.Source.Tag,'coronal Source')==1
    h2.UserData.drag_coronal_source=str2num(varargin{2}.Source.Tag(end));
    startDragFcn_coronal;
end
%% %%%%%%%%  Slice dragging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Axial Slice Dragging
function startDragFcn_axial(varargin)
global h2
h2.WindowButtonMotionFcn=@DragFcn_axial; % start ginput to get cursor location on axes
function DragFcn_axial(varargin)
global h2
xyz=get(gca,'CurrentPoint');
h2.UserData.x_sagittal_slice=round(xyz(3));
h2.UserData.y_coronal_slice=round(xyz(1));
if h2.UserData.drag_axial_source~=0
    v=h2.UserData.drag_axial_source;
    h2.UserData.source_locs(v,1)=round(xyz(3));
    h2.UserData.source_locs(v,2)=round(xyz(1));
end
update_slices;
%% Sagittal Slice Dragging
function startDragFcn_sagittal(varargin)
global h2
h2.WindowButtonMotionFcn=@DragFcn_sagittal; % start ginput to get cursor location on axes
function DragFcn_sagittal(varargin)
global h2
xyz=get(gca,'CurrentPoint');
h2.UserData.z_axial_slice=round(xyz(1));
h2.UserData.y_coronal_slice=round(xyz(3));
if h2.UserData.drag_sagittal_source~=0
    v=h2.UserData.drag_sagittal_source;
    h2.UserData.source_locs(v,3)=round(xyz(1));
    h2.UserData.source_locs(v,2)=round(xyz(3));
end
update_slices;
%% Coronal Slice Dragging
function startDragFcn_coronal(varargin)
global h2
h2.WindowButtonMotionFcn=@DragFcn_coronal; % start ginput to get cursor location on axes
function DragFcn_coronal(varargin)
global h2
xyz=get(gca,'CurrentPoint');
h2.UserData.z_axial_slice=round(xyz(1));
h2.UserData.x_sagittal_slice=round(xyz(3));
if h2.UserData.drag_coronal_source~=0
    v=h2.UserData.drag_coronal_source;
    h2.UserData.source_locs(v,3)=round(xyz(1));
    h2.UserData.source_locs(v,1)=round(xyz(3));
end
update_slices;
%% Select Source location
function select_location(varargin)
global h2
if strcmpi(varargin{2}.Source.String,'Set Source 1')
    h2.UserData.source_locs(1,:)=[h2.UserData.x_sagittal_slice h2.UserData.y_coronal_slice h2.UserData.z_axial_slice];
elseif strcmpi(varargin{2}.Source.String,'Set Source 2')
    h2.UserData.source_locs(2,:)=[h2.UserData.x_sagittal_slice h2.UserData.y_coronal_slice h2.UserData.z_axial_slice];
elseif strcmpi(varargin{2}.Source.String,'Set Source 3')
    h2.UserData.source_locs(3,:)=[h2.UserData.x_sagittal_slice h2.UserData.y_coronal_slice h2.UserData.z_axial_slice];
end
update_slices;
function goto_location(varargin)
global h2
if strcmpi(varargin{2}.Source.String,'Go to Source 1')
    h2.UserData.x_sagittal_slice=h2.UserData.source_locs(1,1);
    h2.UserData.y_coronal_slice=h2.UserData.source_locs(1,2);
    h2.UserData.z_axial_slice=h2.UserData.source_locs(1,3);
elseif strcmpi(varargin{2}.Source.String,'Go to Source 2')
    h2.UserData.x_sagittal_slice=h2.UserData.source_locs(2,1);
    h2.UserData.y_coronal_slice=h2.UserData.source_locs(2,2);
    h2.UserData.z_axial_slice=h2.UserData.source_locs(2,3);
elseif strcmpi(varargin{2}.Source.String,'Go to Source 3')
    h2.UserData.x_sagittal_slice=h2.UserData.source_locs(3,1);
    h2.UserData.y_coronal_slice=h2.UserData.source_locs(3,2);
    h2.UserData.z_axial_slice=h2.UserData.source_locs(3,3);
end
update_slices;
%% save source locs to workspace
function save_location(varargin)
global h2 h
h.cfg.source.vx_idx = find_nearest_voxel(h2.UserData.source_locs_mm,h.anatomy.leadfield.voxel_pos);
h.cfg.study.source_locs_mm = h.anatomy.leadfield.voxel_pos(h.cfg.source.vx_idx,:);
h.cfg.study.source_locs = ft_warp_apply(inv(h.anatomy.mri.transform),h.cfg.study.source_locs_mm);
% h.cfg.study.source_locs=h2.UserData.source_locs; h.cfg.study.source_locs_mm=h2.UserData.source_locs_mm;
h.cfg.source.vx_locs = h.cfg.study.source_locs_mm;
h2.UserData.source_locs = h.cfg.study.source_locs;
h2.UserData.source_locs_mm = h.cfg.study.source_locs_mm;
for v=1:3; h.edit_source_locs(v).String = sprintf('%.1f %.1f %.1f', h.cfg.study.source_locs_mm(v,:)); end
update_slices; update_source_data(); %plot_3D_mri();

% h.tab_power.Visible='on';
%% Plot source_locs_mm in 3D volume "vol"
function plot_source_locs_3D(varargin)
global h2
if ~isfield(h2.UserData.opt,'axes_h')
    figure; set(gcf,'color','w'); cla; ax=gca; axis off;
    h2.UserData.opt.axes_h=ax;
else
    if ~isvalid(h2.UserData.opt.axes_h)
        figure; set(gcf,'color','w'); cla; ax=gca;  axis off;
        h2.UserData.opt.axes_h=ax;
    end
end
h2.UserData.opt.source_clr=h2.UserData.source_clr;
h2.UserData.opt.vol_nums=1:length(h2.UserData.vol);
bl_plot_source_locs(h2.UserData.vol,h2.UserData.source_locs_mm,h2.UserData.opt);
function default_location(varargin)
global h2
h2.UserData.source_locs = [88 57 84; 47 114 84; 146 112 85]; %ones(3,3); %[];
update_slices();


