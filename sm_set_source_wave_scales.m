function sm_set_source_wave_scales(varargin)
% Function updates axis scales for source waveform plots
global h
switch varargin{end}
    case 'Xscale'
        h.axes_source_waves.XLim = str2num(h.edit_inv_source_waves_xscale.String); 
    case 'Yscale'
        h.axes_source_waves.YLim = str2num(h.edit_inv_source_waves_yscale.String); 
    case 'Both'
        h.axes_source_waves.XLim = str2num(h.edit_inv_source_waves_xscale.String); 
        h.axes_source_waves.YLim = str2num(h.edit_inv_source_waves_yscale.String); 
end
