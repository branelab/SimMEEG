function listbox_inv_solns_Callback(varargin)
global h

h.current_inv_soln = h.listbox_inv_solns.Value;
if length(h.current_inv_soln)==1
    bs_plot_inv_soln;
else
end

