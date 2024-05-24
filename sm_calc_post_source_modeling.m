function sm_calc_post_source_modeling(varargin)
global h


if isempty(h.inv_soln(h.current_inv_soln).soln)
    warndlg(sprintf('%s solution not computed\nCheck matlab command line for errors',h.inv_soln(h.current_inv_soln).Type),'Warning!');
else

    % sMCMV can return "inf" so setting these to zero
    h.inv_soln(h.current_inv_soln).soln.P.img(isinf(h.inv_soln(h.current_inv_soln).soln.P.img)) = 0;

    h.inv_soln(h.next_inv_soln).soln.compute_time = toc(h.tt0); % time to compute inv_soln

    %% Listbox Name

    %% original map that gets changed with spatiotemporal mapping
    h.inv_soln(h.current_inv_soln).org_img = h.inv_soln(h.current_inv_soln).soln.P.img;

    %% calculating voxel-by-voxel spread dispersion functions for point-source spread functions (PSF) and cross-talk functions (CTF)
    if h.radio_3D_calc_SPF_CTF.Value==1
        fprintf('Calculating voxel-by-voxel spread dispersion functions for point-source spread functions (PSF) and cross-talk functions (CTF). This may take some time ...\n');
        [h.inv_soln(h.current_inv_soln).SD_PSF, h.inv_soln(h.current_inv_soln).SD_CTF, h.inv_soln(h.current_inv_soln).Rm, h.inv_soln(h.current_inv_soln).nRm] = sm_calc_SPF_CTF_inv_soln(h.inv_soln(h.current_inv_soln).leadfield.H,h.inv_soln(h.current_inv_soln).soln.ori,...
            h.inv_soln(h.current_inv_soln).soln.wts,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos);
    else
        h.inv_soln(h.current_inv_soln).SD_PSF=[];
        h.inv_soln(h.current_inv_soln).SD_CTF=[];
        h.inv_soln(h.current_inv_soln).Rm=[];
        h.inv_soln(h.current_inv_soln).nRm=[];
    end

    % min_max = round([min(h.inv_soln(h.current_inv_soln).soln.P.img) max(h.inv_soln(h.current_inv_soln).soln.P.img)],3,'significant');
    min_max = [0 max(abs(h.inv_soln(h.current_inv_soln).soln.P.img))];

    %% setting min to null_thresh
    % min_max(1) = h.inv_soln(h.current_inv_soln).soln.null_thresh;
    % if min_max(1)==min_max(2); min_max(2)=min_max(1)+1;
    % elseif min_max(1)>min_max(2); min_max(2)=min_max(1)+1;
    % end

    h.inv_soln(h.current_inv_soln).soln.plot_min_max = real(min_max);
    % h.inv_soln(h.current_inv_soln).soln.plot_thresh = round(h.inv_soln(h.current_inv_soln).soln.plot_min_max(2)*.5,3,'significant'); %h.inv_soln(h.current_inv_soln).soln.null_thresh;
    h.inv_soln(h.current_inv_soln).soln.plot_thresh = h.inv_soln(h.current_inv_soln).soln.null_thresh(end);
    update_image_thresh_txt;

    % updating threshold slider
    h.slider_3D_image_thresh.Max = h.inv_soln(h.current_inv_soln).soln.plot_min_max(2);
    h.inv_soln(h.current_inv_soln).soln.plot_min_max(1) =0; % setting to zero to find all possible then using slider to threshold image <--- this speeds up slider processing so as not to call in bs_plot_inv_soln everytime slider is moved
    h.slider_3D_image_thresh.Min = h.inv_soln(h.current_inv_soln).soln.plot_min_max(1);
    h.slider_3D_image_thresh.Value = 0;
    h.inv_soln(h.current_inv_soln).soln.plot_thresh = 0; % slider value will be updated to this when calling bs_plot_peak_waves.m


    bs_update_3D_listbox();
    bs_plot_inv_soln;

    %% setting h.inv_soln(h.current_inv_soln).peak_voxels once then updating current_3d_peak_voxels using slider threshold
    % for this inv_soln result these remain unchanged moving forward
    % h.inv_soln(h.current_inv_soln).peak_idx = h.current_3D_peak_idx;
    % h.inv_soln(h.current_inv_soln).peak_voxels = h.current_3D_peak_voxels;

    %% updating slider threshold to null_thresh
    h.inv_soln(h.current_inv_soln).soln.plot_min_max(2) = max(abs(h.inv_soln(h.current_inv_soln).soln.P.img));
    h.inv_soln(h.current_inv_soln).soln.plot_min_max(1) = h.inv_soln(h.current_inv_soln).soln.null_thresh(end);
    h.slider_3D_image_thresh.Min = h.inv_soln(h.current_inv_soln).soln.plot_min_max(1);
    h.inv_soln(h.current_inv_soln).soln.plot_min_max(1) = h.inv_soln(h.current_inv_soln).soln.null_thresh(end); % setting to zero to find all possible then using slider to threshold image <--- this speeds up slider processing so as not to call in bs_plot_inv_soln everytime slider is moved
    h.slider_3D_image_thresh.Value = h.inv_soln(h.current_inv_soln).soln.null_thresh(end);
    set_3D_image_thresh; % re-thresholding data based on null_thresh

    % setting slider values to be within range of uicontrol min max
    h.slider_3D_image_thresh_text_max.String = num2str(h.slider_3D_image_thresh.Max);
    if h.slider_3D_image_thresh.Value>h.slider_3D_image_thresh.Max
        h.slider_3D_image_thresh.Value = h.slider_3D_image_thresh.Max;
    elseif h.slider_3D_image_thresh.Value < h.slider_3D_image_thresh.Min
        h.slider_3D_image_thresh.Value = h.slider_3D_image_thresh.Min;
    end

    % hm = msgbox(sprintf('Calculating Performance Results\n\nThis may take time depending on how many initial peaks found\n'),'Performance Results');
    % sm_calc_results_loc_performance(); % calculates thresholded MCC and FPR/TPR results
    % close(hm);

    fprintf('Inverse Modeling Completed for %s\n',h.inv_soln(h.next_inv_soln).Type);

    %% converting to single precision for reduce memory storage
    h.inv_soln(h.next_inv_soln).soln.wts = single(h.inv_soln(h.next_inv_soln).soln.wts);
    h.inv_soln(h.next_inv_soln).soln.ori = single(h.inv_soln(h.next_inv_soln).soln.ori);
    h.inv_soln(h.next_inv_soln).soln.ori = single(h.inv_soln(h.next_inv_soln).soln.ori);
    h.inv_soln(h.next_inv_soln).soln.P.img = single(h.inv_soln(h.next_inv_soln).soln.P.img);
    % reducing redundancy to save storage space
    if isfield(h.inv_soln(h.next_inv_soln).leadfield,'cfg')
        h.inv_soln(h.next_inv_soln).leadfield = rmfield(h.inv_soln(h.next_inv_soln).leadfield,{'cfg'});
    end

    h.inv_soln(h.next_inv_soln).leadfield = double2single(h.inv_soln(h.next_inv_soln).leadfield);
    h.inv_soln(h.current_inv_soln).TFR_results =[];

    h.run_inv_soln_flag = 0; % user pressed "run modeling" button
    set_3D_image_thresh();
    bs_calc_errors_inv_soln(); % final update to plot bar graphs

    %% cleaning up
    axis(h.axes_3D_images,'off');
    % parfevalOnAll(@clearvars, 0) % clears parfor memory
    if h.monte_carlo_flag ~= 1
        h.waitfor_panel.Visible='off'; h.waitfor_txt.String = sprintf('Default Message');
    end

end
