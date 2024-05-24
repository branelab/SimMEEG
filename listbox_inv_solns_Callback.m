function listbox_inv_solns_Callback(varargin)
global h

if ischar(varargin{end})
    if contains(varargin{end},'dipole'); h = evalin("base",'h'); end % for calling in from base when Brainstorm started SimMEEG
end
h.current_inv_soln = h.listbox_inv_solns.Value;
if isscalar(h.current_inv_soln)
    h.find_inv_hits_flag = 0; % setting this so that switch from one InvSoln to the next doesn't always search for hits again. reducing redundancies.
    bs_plot_inv_soln;
else
end

