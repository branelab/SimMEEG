function bs_filter_data_Callback(varargin)
global h

combine_sens_data;

lpfreq = str2num(h.edit_filt_lowpass.String);
hpfreq = str2num(h.edit_filt_highpass.String);

TB(1) = hpfreq*.2;  % transition band in 20% of freq
TB(2) = lpfreq*1.1;  % transition band in 20% of freq
freqs = [TB(1) hpfreq lpfreq TB(2)];

filt_type = h.menu_filt_type.String{h.menu_filt_type.Value};

if h.menu_filt_type.Value <=3   % FIR
    filt_method = 'kaiserwin';
elseif h.menu_filt_type.Value >=4   % IIR
    filt_method = 'butter';
end

%% Field Trip's external directory with filtfilt messes up BRANE Lab's filter_data_new.m function need to remove this directory when filtering

hm = msgbox('Filtering Data ...');
%% this filter seems to reduce the channel rank
if sum(h.current_filt_freqs==0)<4 && sum(h.current_filt_freqs-freqs)==0 &&  strcmp(h.current_filt_filt_type,filt_type) % filter options have not changed
    [h.filt_design,h.sim_data.sens_final]=filter_data_new(h.filt_design,double(h.sim_data.sens_final)); % filtering using alread set filter design
else        % update filter design and fitler the data
    fprintf('recalculating filter design ...\n');
    [h.filt_design,h.sim_data.sens_final]=filter_data_new([],double(h.sim_data.sens_final),h.sim_data.cfg.study.srate,freqs,filt_type,filt_method,1,12,[],1);
    h.current_filt_freqs = freqs;
    h.current_filt_filt_type = filt_type;
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
