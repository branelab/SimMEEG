function bs_plot_map_time(varargin)
global h

if h.btn_3D_plot_peak_waves.Value == 1
    h.axes_source_waves.ButtonDownFcn = @bs_plot_map_time;
else
    h.axes_source_waves.ButtonDownFcn = [];
end

if length(h.listbox_inv_solns.Value)==1
    if h.btn_3D_plot_peak_waves.Value==0    % plot original map
        %     if isfield(h.inv_soln(h.current_inv_soln),'org_img')    % stored orignal inverse solution map
        %         h.inv_soln(h.current_inv_soln).org_img = h.inv_soln(h.current_inv_soln).soln.P.img;
        %     else
        h.inv_soln(h.current_inv_soln).soln.P.img = h.inv_soln(h.current_inv_soln).org_img;
        %     end
        bs_plot_inv_soln;
    elseif h.btn_3D_plot_peak_waves.Value==1    % plot spatiotemporal map
        
        [xy]=get(h.axes_source_waves,'CurrentPoint');
        h.current_swf_time = xy;
        if isfield(h,'current_swf_time_plot')
            if isvalid(h.current_swf_time_plot)
                h.current_swf_time_plot.XData = [xy(1) xy(1)];
                h.current_swf_time_plot.YData = h.axes_source_waves.YLim;
            else
                h.current_swf_time_plot = plot(h.axes_source_waves,[xy(1) xy(1)], h.axes_source_waves.YLim,'r--');
            end
        else
            h.current_swf_time_plot = plot(h.axes_source_waves,[xy(1) xy(1)], h.axes_source_waves.YLim,'r--');
        end
        lat = h.cfg.study.lat_sim;
        try
            ss = find(lat<=h.current_swf_time_plot.XData(1)); s = ss(end); if s<=0; s=1; elseif s>length(lat); s=length(lat); end
        catch
            ss = find(lat<=0); s = ss(end); if s<=0; s=1; elseif s>length(lat); s=length(lat); end
        end
        
        
        
        switch h.inv_soln(h.current_inv_soln).Type
            case {'SPA' 'SIA' 'MIA' 'LCMV (FT)' 'sLORETA (FT)' 'SAM (FT)' 'sMCMV' 'bRAPBeam' 'TrapMUSIC'}    % BRANE Lab beamformers
                swf = abs(h.inv_soln(h.current_inv_soln).soln.wts' * squeeze(nanmean(h.sim_data.sens_final(s,h.anatomy.sens.good_sensors,:),3))');
                if h.radio_normalize_swf.Value == 1; swf_base = abs(h.inv_soln(h.current_inv_soln).soln.wts' * squeeze(nanmean(h.sim_data.sens_final(h.sim_data.cfg.study.base_samps,h.anatomy.sens.good_sensors,:),3))'); end
            case {'eLORETA (FT)' 'MNE (FT)'}    % Field Trips inverse solutions
                swf=[];
                for ox = 1:size(h.inv_soln(h.current_inv_soln).soln.wts,3)
                    swf(ox,:)=squeeze(nanmean(h.sim_data.sens_final(s,h.anatomy.sens.good_sensors,:),3))*squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,:,ox));
                    if h.radio_normalize_swf.Value == 1; swf_base(:,ox,:)=squeeze(nanmean(h.sim_data.sens_final(h.sim_data.cfg.study.base_samps,h.anatomy.sens.good_sensors,:),3))*squeeze(h.inv_soln(h.current_inv_soln).soln.wts(:,:,ox)); end
                end
                swf = squeeze(rms(swf,1))'; % taking RMS of waveforms across orientations
                if h.radio_normalize_swf.Value == 1; swf_base = squeeze(rms(swf_base,1))'; end % taking RMS of waveforms across orientations
            case {'LCMV (BST)' 'MNE (BST)' 'sLORETA (BST)'}
                
                if isempty(h.inv_soln(h.current_inv_soln).soln.ImageGridAmp)
                    swf = h.inv_soln(h.current_inv_soln).soln.ImagingKernel * squeeze(nanmean(h.sim_data.sens_final(s,h.anatomy.sens.good_sensors,:),3))';
                else
                    swf = h.inv_soln(h.current_inv_soln).soln.ImageGridAmp(:,s);
                end
                
                if h.radio_normalize_swf.Value == 1; swf_base = h.inv_soln(h.current_inv_soln).soln.ImagingKernel * squeeze(nanmean(h.sim_data.sens_final(h.sim_data.cfg.study.base_samps,h.anatomy.sens.good_sensors,:),3))'; end
                switch h.inv_soln(h.current_inv_soln).maxvectorori_Type
                    case 'RMS' % transforming vector dipole orientation with rms img power between active and control interval into scalar image
                        iVertSource = 1:size(swf,1);
                        idx1 = 1:3:length(iVertSource); idx2 = 2:3:length(iVertSource); idx3 = 3:3:length(iVertSource);
                        swf = squeeze(sqrt(swf(idx1,:).^2 + swf(idx2,:).^2 + swf(idx3,:).^2));
                        if h.radio_normalize_swf.Value == 1; swf_base = squeeze(sqrt(swf_base(idx1,:).^2 + swf_base(idx2,:).^2 + swf_base(idx3,:).^2)); end
                    case 'Max'    % transforming vector dipole orientation with maximum img power between active and control interval into scalar image
                        iVertSource = 1:size(swf,1);
                        idx1 = 1:3:length(iVertSource); idx2 = 2:3:length(iVertSource); idx3 = 3:3:length(iVertSource);
                        swf = max(abs(squeeze(cat(2, swf(idx1,:), swf(idx2,:), swf(idx3,:)))), [], 2);
                        if h.radio_normalize_swf.Value == 1; swf_base = max(abs(squeeze(cat(3, swf_base(idx1,:), swf_base(idx2,:), swf_base(idx3,:)))), [], 3); end
                end
        end
        
        if h.radio_normalize_swf.Value == 1 % plot normalized waves
            %% normalizing z-transform relative to standard deviation of the baseline - This is not correct for some inv solutions because they have already normalized to baseline
            mu_base = squeeze(nanmean(swf_base,2));
            sigma_base = squeeze(nanstd(swf_base,[],2));
            
            norm_swf = (swf-mu_base) ./ sigma_base;    % Z-score relative to baseline
            h.inv_soln(h.current_inv_soln).soln.P.img = abs(norm_swf);
        else
            h.inv_soln(h.current_inv_soln).soln.P.img = abs(swf);
        end
        bs_plot_inv_soln;
        
        %     axes(h.axes_source_waves);
        title(h.axes_source_waves, sprintf('%s (%.f ms)',h.inv_soln(h.current_inv_soln).Type,h.current_swf_time(1)*1000));
    end
else
    msgbox('please select only 1 Inveres Solution in List to plot');
end
