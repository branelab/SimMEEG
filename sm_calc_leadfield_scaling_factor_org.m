function sm_calc_leadfield_scaling_factor(varargin)
% This function loads in a predetermined wts matrix and source location from BESA to 
%   calculate the leadfield scaling factor of FieldTrip's leadfield calculations or any loaded in Anatomy files with lead fields

global h
h.edit_leadfield_gain.String = '1'; 

if h.menu_sens_type.Value == 1 % MEG
    %% MEG scaling factor
    sens_file = fullfile(h.simmeeg_dir,'G003_SESAEF_tone_on_av_MEG_data.avr');  % sensor data file
    source_file = fullfile(h.simmeeg_dir,'G003_SESAEF_tone_on_av_source_data.avr'); % source data file
    dip_file = fullfile(h.simmeeg_dir,'G003_SESAEF_tone_on_av.bsa');
    sens_type = 'MEG';
    sens_convert = 1e15;    % convert to femtoT
elseif h.menu_sens_type.Value == 2   % EEG
    sens_file = fullfile(h.simmeeg_dir,'GT_1000Hz_av_EEG_data.avr');  % sensor data file
    source_file = fullfile(h.simmeeg_dir,'GT_1000Hz_av_source_data.avr'); % source data file
    dip_file = fullfile(h.simmeeg_dir,'GT_1000Hz_av.bsa');
    sens_type = 'EEG';
    sens_convert = 1e6;    % convert to microV
end

%% get BESA Leadfield wts 
[wts_besa, sens_data, source_data, sens_names] = bl_create_wts_from_BESA_files(sens_file, source_file, sens_type, 0); % source data are in Amps

% [sens_data,sens_names,hdr_line]=read_BESA_files(sens_file); % if MEG then units = femtoT; if EEG then units = microV
% [source_data,source_names,hdr_line]=read_BESA_files(source_file); % units = nano-Amp


%% Load dip locations
[dip_locs,source_names,hdr_line]=read_BESA_files(dip_file);
% location axes might be rotated 90 degrees
loc_order = [2 1 3];


h.leadfield_wts_ratio =[];
for v=1:size(dip_locs,1)
    % dip_ori = [-158.8 84.9; 140.6 -87.6];
    [vidx(v)] = find_nearest_voxel(dip_locs(v,loc_order),h.anatomy.leadfield.voxel_pos);

%     
%     %% need to sum the wts_lf
%     wts_lf2(:,:,v) = h.anatomy.leadfield.H(:,:,vidx(v));
%     wts_lf(:,v) = sum(wts_lf2(:,:,v),2)/3;
%     
%     %% find sensor with maximum wts
%     [~,max_idx] = max(abs(wts_besa(:,v)));
%     h.leadfield_wts_ratio(v) = wts_besa(max_idx,v)./wts_lf(max_idx,v); % divind by 3 because there are 3 summed vectors
end

    source_amp = [1 1]; % units = nAmps
%     data = repmat(source_data',[1 1 2])*1e9; % need trials in data 
    data = repmat(source_data',[1 1 2]); % need trials in data 
    proj_data = project_SimSignals(data,h.anatomy.leadfield,vidx,source_amp,dip_locs(:,4:6));   % source data units must be in Amps to return sensor data in Tesla (MEG) or Volts (EEG)
    [q,v]= max(max(abs(sens_data'))); % sensor with maximum amplitude across time
    [q,s]= max(abs(sens_data(v,:))); % time sample with maximum amplitude across time
    
    %% match sensor names in case of using different channel sizes
    v_idx = [];
    %     for v=1:length(sens_names)
    if ~isempty(find(strcmpi(h.anatomy.leadfield.label,  sens_names(v))==1, 1))
        v_lf = find(strcmpi(h.anatomy.leadfield.label,  sens_names(v))==1);
        v_idx = v;
%         lf_gain = 3* squeeze(proj_data(s,v_lf,1)) ./ (squeeze(sens_data(v_idx,s))') ;
        lf_gain = (squeeze(sens_data(v_idx,s))') ./ squeeze(proj_data(s,v_lf,1)) ;
    else
        lf_gain = 1;
    end
%     end

% h.edit_leadfield_gain.String = num2str(round(max(lf_gain)));
h.edit_leadfield_gain.String = num2str(round(abs(lf_gain)));
    
end

