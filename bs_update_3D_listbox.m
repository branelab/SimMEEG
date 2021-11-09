function bs_update_3D_listbox(varargin)
global h

if ~isempty(h.inv_soln)
    h.listbox_inv_solns.String = {h.inv_soln.ListBox_name};
    h.listbox_inv_solns.Value = length(h.inv_soln); 
else
     h.listbox_inv_solns.String = '';
    h.listbox_inv_solns.Value = 1;
end
