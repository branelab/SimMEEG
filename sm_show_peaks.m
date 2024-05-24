function sm_show_peaks(p_idx,vis_on)
% This function hides scattter locs and waveforms of the inv_soln peaks indexed as hide_idx = [1 2 3 ...] within the h.inv_soln(h.current_inv_soln).peak_voxels
global h

%% need to find peaks that are visible
p_idx2 = find(isprop(h.map3D_peak_locs,'Marker')); 
if ~isempty(p_idx2)
    p_idx = p_idx2(p_idx);
else
    p_idx = [];
end

%  set(h.func3D(:),'Visible','off');
% h.map3D_peak_ori = handle(h.map3D_peak_ori);
if ~isempty(h.map3D_peak_locs) && ~isempty(h.map3D_peak_ori) && ~isempty(p_idx)
%     try
        set(h.map3D_peak_locs(p_idx),'Visible',vis_on);
        set(h.map3D_peak_ori(p_idx),'Visible',vis_on);
        %% slices correspnding to peak locations
        % dimensions for func3D are surfaces in the order of [xslice yslice zslice] for each location
        % dims = [3 size(h.func3D,1)/3]
        % func3d_idx =
        xlocs = [h.map3D_peak_locs(p_idx).XData];   % getting plane locations from p_idx scatter plots
        ylocs = [h.map3D_peak_locs(p_idx).YData];
        zlocs = [h.map3D_peak_locs(p_idx).ZData];

        % finding slices above threshold
        xslice_idx = []; yslice_idx = []; zslice_idx = [];
        
        for s = 1:length(h.func3D)
            if ismember(unique(h.func3D(s).XData),xlocs) % found slice at all thresholded peaks in xlocs
                xslice_idx = [xslice_idx s];
            end
             if ismember(unique(h.func3D(s).YData),ylocs) % found slice at all thresholded peaks in xlocs
                yslice_idx = [yslice_idx s];
            end
            if ismember(unique(h.func3D(s).ZData),zlocs) % found slice at all thresholded peaks in xlocs
                zslice_idx = [zslice_idx s];
            end
       end
        
%% setting visibility for slices 
        set(h.func3D(xslice_idx),'Visible',vis_on);
        set(h.func3D(yslice_idx),'Visible',vis_on);
        set(h.func3D(zslice_idx),'Visible',vis_on);
try
        %% waveforms and FFT
        set(h.current_inv_swf_plots(p_idx),'Visible',vis_on);
        set(h.current_inv_swf_plots_fft(p_idx),'Visible',vis_on);
end     

    
end


