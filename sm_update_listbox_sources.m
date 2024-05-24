function sm_update_listbox_sources(varargin)
% update source listbox with all sources created
global h

if max(h.listbox_sources.Value) > length(h.listbox_sources.String); h.listbox_sources.Value = 1; end

pre = '<HTML><FONT color="'; post = '</FONT></HTML>';
h.listbox_sources.UserData = rmfield(h.listbox_sources.UserData,'clr_str'); 

for v=1:length(h.cfg.source.vx_idx)
    
    % changing color of list box items for sources 1, 2, and 3
    Data(v).name = sprintf('Source %.f',v); Data(v).Color = round(h.cfg.source.src_clr(v,:)*255);
    
    clr_str = reshape( dec2hex( Data(v).Color, 2 )',1, 6);
    h.listbox_sources.UserData.clr_str{v} = [pre clr_str '">' Data(v).name post];
end
    h.listbox_sources.String = h.listbox_sources.UserData.clr_str;
    
    