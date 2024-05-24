function [hwin] = sm_create_amp_window(win_type, lat, srate, sig_time, rise_time, fall_time, pre_perc, sig_perc, post_perc)
% function [hwin] = sm_create_amp_window(win_type, lat, srate, sig_time, rise_time, fall_time, pre_perc, sig_perc, post_perc);
%   Calculates window for ampltidue modulation/shaping
%
% INPUT
%   win_type = type of window: 'Hanning' 'Gaussian' 'Triangular','Blackmann'
%   lat = laetncy values of samples for the window
%   srate = sample rate
%   sigtime = signal times [start end]
%   risetime = rise time
%   falltime = rise time
%   pre_perc = percent amplitude for pre-signal event
%   signal_perc = percent amplitude for signal plateau event
%   pre_perc = percent amplitude for post-signal event
%

plat_time = sig_time(2)-sig_time(1) - rise_time - fall_time;
if plat_time<0   % if rise_time + fall_time exceed sig_time the make rise/fall times half of sig_time.
    rise_time = (sigtime(2)-sigtime(1))/2;
    fall_time = rise_time;
end

rise_samps=round(2*rise_time*srate);
fall_samps=round(2*fall_time*srate);
plat_samps = round(plat_time*srate); % plateau time

switch win_type
    case 'Hanning'
        rise_win = hanning(rise_samps);
        fall_win = hanning(fall_samps);
    case 'Triangular'
        rise_win = triang(rise_samps);
        fall_win = triang(fall_samps);
    case 'Gaussian'
        rise_win = gausswin(rise_samps);
        fall_win = gausswin(fall_samps);
    case 'Blackmann'
        rise_win = blackman(rise_samps);
        fall_win = blackman(fall_samps);
end
%% taking half of rise/fall windows for edges
rise_win2 = rise_win(1:floor(length(rise_win)/2));
fall_win2 = fall_win(end-floor(length(fall_win)/2):end);
plat_win2 = ones(plat_samps,1)*sig_perc/100;

if pre_perc>sig_perc % ERD -- need to flip rise_win2 upside down
    rise_win2 = flipud(rise_win2);
end
if sig_perc<post_perc % ERD -- need to flip rise_win2 upside down
    fall_win2 = flipud(fall_win2);
end

% scaling windows from 0-1
rise_win2=rise_win2-min(rise_win2); rise_win2=rise_win2/max(rise_win2);
rise_win2 = ( (rise_win2*(sig_perc-pre_perc)) + pre_perc )/100;

fall_win2=fall_win2-min(fall_win2); fall_win2=fall_win2/max(fall_win2);
fall_win2 = ( (fall_win2*(sig_perc-post_perc)) + post_perc )/100;

hwin2 = [rise_win2; plat_win2; fall_win2];

ss = find(lat<=sig_time(1)); sx(1)=ss(end);

pre_win = ones(sx(1),1)*pre_perc/100; 
post_win = ones(length(lat)-length(hwin2)-length(pre_win),1)*post_perc/100; 

hwin = [pre_win; hwin2; post_win];


