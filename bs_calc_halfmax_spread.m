function [half_width,xslices,yslices,zslices]=bs_calc_halfmax_spread(img,voxel_pos,peak_idx)

% img = h.inv_soln(c).soln.P.img;
% voxel_pos = h.inv_soln(c).leadfield.voxel_pos;
% peak_locs = h.inv_soln(c).leadfield.voxel_pos(h.current_3D_peak_idx,:);
voxel_pos = double(voxel_pos);
peak_locs = double(voxel_pos(peak_idx,:));
img = double(img); 

% clear xslices yslices zslices
for v=1:length(peak_idx)
    xslices(v).idx = find( voxel_pos(:,2)== peak_locs(v,2) & voxel_pos(:,3)== peak_locs(v,3) )';
    yslices(v).idx = find( voxel_pos(:,1)== peak_locs(v,1) & voxel_pos(:,3)== peak_locs(v,3) )';
    zslices(v).idx = find( voxel_pos(:,1)== peak_locs(v,1) & voxel_pos(:,2)== peak_locs(v,2) )';
    xslices(v).vals = img(xslices(v).idx);
    yslices(v).vals = img(yslices(v).idx);
    zslices(v).vals = img(zslices(v).idx);
    
    %% find index in vals of maximum at the selected peak
    peak_xidx = find(xslices(v).idx == peak_idx(v));
    peak_yidx = find(yslices(v).idx == peak_idx(v));
    peak_zidx = find(zslices(v).idx == peak_idx(v));
    
    %% upsampling data using fit functions to get half-max width
%     up_samp=100;
     wx=[]; wy=[]; wz=[];
    
    % X-Slice
    try
%         xvals=interp(xslices(v).vals,up_samp); xdata=interp(voxel_pos(xslices(v).idx,1),up_samp);
        xdata = min(voxel_pos(xslices(v).idx,1)):.1:max(voxel_pos(xslices(v).idx,1));
        ft = fit(voxel_pos(xslices(v).idx,1), xslices(v).vals,'smoothingspline','SmoothingParam',1);
        xvals = feval(ft, xdata);

        % finding half-max edges for multiple peaks on same plane
        % X-Slice
        hmax_vals = xslices(v).vals(peak_xidx)/2;
        tvals = xvals; tvals(tvals<=hmax_vals) = hmax_vals*4;  tvals([1 end]) = hmax_vals*4;
        [~,w]=findpeaks(-tvals);
        
        if peak_xidx==1 % peak is on the outer edge of the head model
            w = [1; w]; half_factor=2; % need to double width because only half of the spread is acounted for
        elseif peak_xidx==size(xslices(v).vals,1) % peak is on the outer edge of the head model
             w = [w; size(xdata,1)]; half_factor=2; % need to double width because only half of the spread is acounted for
        else
            half_factor=1;
        end
        for s=1:2:length(w)     % finding peak that is sandwiched between half-max edges
            if xdata(w(s)) <= voxel_pos(xslices(v).idx(peak_xidx),1) && xdata(w(s+1)) >= voxel_pos(xslices(v).idx(peak_xidx),1)
                wx = s;
            end
        end
        half_width(v).edges_x = [xdata(w(wx)) xdata(w(wx+1))];
        half_width(v).val_x = range(half_width(v).edges_x)*half_factor;
        
    catch me
        half_width(v).edges_x = [nan nan];
        half_width(v).val_x = nan;
    end
    
    % Y-Slice
    try
%         yvals=interp(yslices(v).vals,up_samp); ydata=interp(voxel_pos(yslices(v).idx,2),up_samp);
        ydata = min(voxel_pos(yslices(v).idx,2)):.1:max(voxel_pos(yslices(v).idx,2));
        ft = fit(voxel_pos(yslices(v).idx,2), yslices(v).vals,'smoothingspline','SmoothingParam',1);
        yvals = feval(ft, ydata);
        
        % Y-Slice
        hmax_vals = yslices(v).vals(peak_yidx)/2;
        tvals = yvals; tvals(tvals<=hmax_vals) = hmax_vals*4;   tvals([1 end]) = hmax_vals*4;
        [~,w]=findpeaks(-tvals);
        
         if peak_yidx==1 % peak is on the outer edge of the head model
            w = [1; w]; half_factor=2; % need to double width because only half of the spread is acounted for
        elseif peak_yidx==size(yslices(v).vals,1) % peak is on the outer edge of the head model
             w = [w; size(ydata,1)]; half_factor=2; % need to double width because only half of the spread is acounted for
        else
            half_factor=1;
         end
         for s=1:2:length(w)     % finding peak that is sandwiched between half-max edges
            if ydata(w(s)) <= voxel_pos(yslices(v).idx(peak_yidx),2) && ydata(w(s+1)) >= voxel_pos(yslices(v).idx(peak_yidx),2)
                wy = s;
            end
        end
        half_width(v).edges_y = [ydata(w(wy)) ydata(w(wy+1))];
        half_width(v).val_y = range(half_width(v).edges_y)*half_factor;
    catch me
        half_width(v).edges_y = [nan nan];
        half_width(v).val_y = nan;
    end
    
    % Z-Slice
    try
%         zvals=interp(zslices(v).vals,up_samp); zdata=interp(voxel_pos(zslices(v).idx,3),up_samp);
        zdata = min(voxel_pos(zslices(v).idx,3)):.1:max(voxel_pos(zslices(v).idx,3));
        ft = fit(voxel_pos(zslices(v).idx,3), zslices(v).vals,'smoothingspline','SmoothingParam',1);
        zvals = feval(ft, zdata);

        % Z-Slice
        hmax_vals = zslices(v).vals(peak_zidx)/2;
        tvals = zvals; tvals(tvals<=hmax_vals) = hmax_vals*4;   tvals([1 end]) = hmax_vals*4;
        [~,w]=findpeaks(-tvals);
        if peak_zidx==1 % peak is on the outer edge of the head model
            w = [1; w]; half_factor=2; % need to double width because only half of the spread is acounted for
        elseif peak_zidx==size(zslices(v).vals,1) % peak is on the outer edge of the head model
            w = [w; size(zdata,1)]; half_factor=2; % need to double width because only half of the spread is acounted for
        else
            half_factor=1;
        end
         
         for s=1:2:length(w)     % finding peak that is sandwiched between half-max edges
            if zdata(w(s)) <= voxel_pos(zslices(v).idx(peak_zidx),3) && zdata(w(s+1)) >= voxel_pos(zslices(v).idx(peak_zidx),3)
                wz = s;
            end
        end
        half_width(v).edges_z = [zdata(w(wz)) zdata(w(wz+1))];
        half_width(v).val_z = range(half_width(v).edges_z);
    catch me
        half_width(v).edges_z = [nan nan];
        half_width(v).val_z = nan;
    end
    %% average source leakage
    half_width(v).average = nanmean([half_width(v).val_x half_width(v).val_y half_width(v).val_z]);
end
