function sm_fig_fontsize_gain(varargin)
global h

FontSize_gain = str2num(h.edit_fontsize_gain.String); 
htext = findobj(h.main_fig, '-property', 'FontSize'); 
for f = 1:length(htext); htext(f).FontSize = htext(f).FontSize*FontSize_gain; end

