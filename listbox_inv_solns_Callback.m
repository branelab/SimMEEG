function listbox_inv_solns_Callback(varargin)
global h

h.current_inv_soln = h.listbox_inv_solns.Value;
if length(h.current_inv_soln)==1
    h.find_inv_hits_flag = 0; % setting this so that switch from one InvSoln to the next doesn't always search for hits again. reducing redundancies.
    bs_plot_inv_soln;
else
end

