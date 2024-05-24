function h = sm_batch_sm_plot_replace_3Dmap(h)


%% Getting selected vx_idx from listbox_peaks_found
for v=1:length(h.listbox_peaks_found.Value)
    vx_str = h.listbox_peaks_found.String{h.listbox_peaks_found.Value(v)};
    st1 = findstr(vx_str,'">');
    st2 = findstr(vx_str,'@');
    if isempty(str2num(vx_str(st1+2:st2-2)))
        vx_idx(v) = nan;
    else
        vx_idx(v) = str2num(vx_str(st1+2:st2-2));
    end
end
vx_idx = vx_idx(~isnan(vx_idx));

%% 

switch  h.menu_invsoln_map_type.String{h.menu_invsoln_map_type.Value}
    case 'image'
    h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).org_img; 
    case 'normalized image to max(abs)'
        h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).org_img./ max(abs(h.inv_soln(h.current_inv_soln).org_img));
    case 'Point Spread Function (PSF)'
        h.inv_soln(h.current_inv_soln).soln.P.img = nansum(h.inv_soln(h.current_inv_soln).Rm(:,vx_idx),2);
    case 'Cross Talk Function (CTF)'
        h.inv_soln(h.current_inv_soln).soln.P.img = nansum(h.inv_soln(h.current_inv_soln).Rm(vx_idx,:),1)';
    case 'normalized PSF'
        h.inv_soln(h.current_inv_soln).soln.P.img = nansum(h.inv_soln(h.current_inv_soln).nRm(:,vx_idx),2);
    case 'normalized CTF'
        h.inv_soln(h.current_inv_soln).soln.P.img = nansum(h.inv_soln(h.current_inv_soln).nRm(vx_idx,:),1)';
    case 'Spatial Dispersion PSF'
        h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).SD_PSF';
    case 'Spatial Dispersion CTF'
        h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).SD_CTF';
    case 'Overall Amplitude PSF'
        oa = nansum(abs(h.inv_soln(h.current_inv_soln).Rm(:,:)),2); 
        oa = oa/max(abs(oa)); % normalizing max(abs) so yield 0 to 1 scale
         h.inv_soln(h.current_inv_soln).soln.P.img = oa;   % sum of overall PSFs across whole brain   
    case 'Overall Amplitude CTF'
        oa = nansum(abs(h.inv_soln(h.current_inv_soln).Rm(:,:)),1)';
        oa = oa/max(abs(oa)); % normalizing max(abs) so yield 0 to 1 scale
         h.inv_soln(h.current_inv_soln).soln.P.img = oa;   % sum of overall CTFs across whole brain   
end

h.inv_soln(h.current_inv_soln).soln.P.img_org = h.inv_soln(h.current_inv_soln).soln.P.img; 




