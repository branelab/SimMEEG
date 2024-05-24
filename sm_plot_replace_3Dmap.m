function sm_plot_replace_3Dmap(varargin)

global h

if ~isempty(h.listbox_peaks_found.String)
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
                h.panel_3D_image_plot_msg_txt.String = sprintf('%s\n%s mapping',h.inv_soln(h.current_inv_soln).Type,h.menu_invsoln_map_type.String{h.menu_invsoln_map_type.Value}); drawnow;
    switch  h.menu_invsoln_map_type.String{h.menu_invsoln_map_type.Value}
        case 'image'
            if h.radio_find_spatiotemp_peaks.Value == 0 && h.btn_3D_plot_peak_waves.Value == 0  % Do not revert to org_img when spatiotemp plotting
                h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).org_img; 
            end
            
        case 'normalized image to max(abs)'
            h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).org_img./ max(abs(h.inv_soln(h.current_inv_soln).org_img));
            
        case 'Point Spread Function (PSF)'
            if ~isempty(h.inv_soln(h.current_inv_soln).Rm)
                h.inv_soln(h.current_inv_soln).soln.P.img = abs(nansum(h.inv_soln(h.current_inv_soln).Rm(:,vx_idx),2));
            else
                h.inv_soln(h.current_inv_soln).soln.P.img = zeros(size(h.inv_soln(h.current_inv_soln).soln.P.img));
            end

        case 'Cross Talk Function (CTF)'
            if ~isempty(h.inv_soln(h.current_inv_soln).Rm)
                h.inv_soln(h.current_inv_soln).soln.P.img = abs(nansum(h.inv_soln(h.current_inv_soln).Rm(vx_idx,:),1)');
            end
        case 'normalized PSF'
            if ~isempty(h.inv_soln(h.current_inv_soln).nRm)
                h.inv_soln(h.current_inv_soln).soln.P.img = nansum(h.inv_soln(h.current_inv_soln).nRm(:,vx_idx),2);
            end
        case 'normalized CTF'
            if ~isempty(h.inv_soln(h.current_inv_soln).nRm )
                h.inv_soln(h.current_inv_soln).soln.P.img = nansum(h.inv_soln(h.current_inv_soln).nRm(vx_idx,:),1)';
            end
        case 'Spatial Dispersion PSF'
            if ~isempty(h.inv_soln(h.current_inv_soln).SD_PSF )
                h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).SD_PSF';
            end
        case 'Spatial Dispersion CTF'
            if ~isempty(h.inv_soln(h.current_inv_soln).SD_CTF )
                h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).SD_CTF';
            end
        case 'Overall Amplitude PSF'
            if ~isempty(h.inv_soln(h.current_inv_soln).Rm )
                oa = nansum(abs(h.inv_soln(h.current_inv_soln).Rm(:,:)),2);
                oa = oa/max(abs(oa)); % normalizing max(abs) so yield 0 to 1 scale
                h.inv_soln(h.current_inv_soln).soln.P.img = oa;   % sum of overall PSFs across whole brain
            end
        case 'Overall Amplitude CTF'
            if ~isempty(h.inv_soln(h.current_inv_soln).Rm )
                oa = nansum(abs(h.inv_soln(h.current_inv_soln).Rm(:,:)),1)';
                oa = oa/max(abs(oa)); % normalizing max(abs) so yield 0 to 1 scale
                h.inv_soln(h.current_inv_soln).soln.P.img = oa;   % sum of overall CTFs across whole brain
            end
    end
    
    h.inv_soln(h.current_inv_soln).soln.P.img_org = h.inv_soln(h.current_inv_soln).soln.P.img;
end

% min_max = [min(h.inv_soln(h.current_inv_soln).soln.P.img) max(h.inv_soln(h.current_inv_soln).soln.P.img)];
% h.edit_3D_min_max.String = min_max; h.slider_3D_image_thresh.Value=min_max(1);  set_3D_image_caxis();



