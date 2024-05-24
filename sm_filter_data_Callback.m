function sm_filter_data_Callback(varargin)
global h

%% Need to check to make sure FieldTrip's external directory with filtfilt messes up BRANE Lab's filter_data_new.m function need to remove this directory when filtering
% rmpath(h.FieldTrip_dir_external);
h.btn_reref_EEG.Value = 0; % turn this off because filtering is done on original sens_final data

% combine_sens_data;
if h.monte_carlo_flag == 1
    h.waitfor_txt.String = sprintf('Filtering Data\n\nThis might take some time\n\nPlease wait ...'); drawnow;
else
    h.waitfor_panel.Visible='on';
    h.waitfor_txt.String = sprintf('Filtering Data\n\nThis might take some time\n\nPlease wait ...'); drawnow;
end

switch varargin{1}.String
    case 'Create Filter'
        hm = msgbox('Creating Filter Design...');
        freqs = str2num(h.edit_filt_band.String);
        if freqs(1)==0 && freqs(2)~=0   % cannot have zero for starting frequency
            freqs(1) = .01; h.edit_filt_band.String = sprintf('%.2f %.1f %.1f %.1f', freqs');
        end
        
        %% set filter order for IIR filters if originally set to "0"
        filt_order = str2num(h.edit_filt_order.String);
        if filt_order<2 && h.menu_filt_type.Value>4         % IIR filters
            filt_order = 2;     % set to lowest order which is not optimal so throw a warning
            h.edit_filt_order.String = filt_order;
            warndlg(sprintf('Automatic calculation for IIR filter order is not implemented.\n   Currently set to lowest filter order = 2\n   Please manually set for different filter order\n'),'IIR filter order');
        elseif filt_order==0 && h.menu_filt_type.Value<=4     % FIR filters
            filt_order = [];
        end
        if any(freqs==0); freqs(freqs==0)=nan; end      % setting nans for freqs for filt_type 'lowpass' and 'highpass'
        
        filt_type = h.menu_filt_type.String{h.menu_filt_type.Value};
        filt_method = h.menu_filt_method.String{h.menu_filt_method.Value};
        filt_dB = str2num(h.edit_filt_StopbandAttenuation.String);
        
        [h.filt_design] = filter_data_new([],[],h.cfg.study.srate,freqs,filt_type,filt_method,1,filt_dB,filt_order,h.radio_plot_filter_design.Value);
        
    case 'Filter Data'
        hm = msgbox('Filtering Data ...');
        if isfield(h,'filt_design')
            if ~isempty(h.filt_design)
                if ~isfield(h.sim_data,'sens_final_org'); h.sim_data.sens_final_org = h.sim_data.sens_final; end
                if ~isfield(h.sim_data,'sens_noise_final_org'); h.sim_data.sens_noise_final_org = h.sim_data.sens_noise_final; end
                if ~isfield(h.sim_data,'sens_sig_data_org'); h.sim_data.sens_sig_data_org = h.sim_data.sens_sig_data; end
                if ~isfield(h.sim_data,'sig_final_org'); h.sim_data.sig_final_org = h.sim_data.sig_final; end
                
                h.sim_data.sens_final_org(isnan(h.sim_data.sens_final_org)) = 0; 
                
                [~,h.sim_data.sens_final(:,h.anatomy.sens.good_sensors,:)] = filter_data_new(h.filt_design,double(h.sim_data.sens_final_org(:,h.anatomy.sens.good_sensors,:)));
                [~,h.sim_data.sens_noise_final(:,h.anatomy.sens.good_sensors,:)] = filter_data_new(h.filt_design,double(h.sim_data.sens_noise_final_org(:,h.anatomy.sens.good_sensors,:)) );
                [~,h.sim_data.sens_sig_data(:,h.anatomy.sens.good_sensors,:)] = filter_data_new(h.filt_design,double(h.sim_data.sens_sig_data_org(:,h.anatomy.sens.good_sensors,:)));
                [~,h.sim_data.sig_final] = filter_data_new(h.filt_design,double(h.sim_data.sig_final_org));
            else
                warndlg('Please Create Filter in Preprocessing Panel','No Filter Created');
            end
        else
            warndlg('Please Create Filter in Preprocessing Panel','No Filter Created');
        end
end

if h.monte_carlo_flag ~= 1
    h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
end

% h.filt_design =[];
% if lpfreq>0 && hpfreq>0     % bandpass
%     h.sim_data.sens_final = filter_data(double(h.sim_data.sens_final),h.sim_data.cfg.study.srate,[hpfreq lpfreq],'butter','bandpass',[]);
% elseif lpfreq<=0 && hpfreq>0     % highpass
%     h.sim_data.sens_final = filter_data(double(h.sim_data.sens_final),h.sim_data.cfg.study.srate,[hpfreq],'butter','high',[]);
% elseif lpfreq>0 && hpfreq<=0     % highpass
%     h.sim_data.sens_final = filter_data(double(h.sim_data.sens_final),h.sim_data.cfg.study.srate,[lpfreq],'butter','low',[]);
% end


close(hm);
