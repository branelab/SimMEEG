function [half_width,xslices,yslices,zslices]=bs_calc_halfmax_spread(img,voxel_pos,peak_idx)

% img = h.inv_soln(c).soln.P.img;
% voxel_pos = h.inv_soln(c).leadfield.voxel_pos;
% peak_locs = h.inv_soln(c).leadfield.voxel_pos(h.current_3D_peak_idx,:);

peak_locs = voxel_pos(peak_idx,:);

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
    up_samp=100; 
    % X-Slice
%     ft = fit(voxel_pos(xslices(v).idx,1), xslices(v).vals,'smoothingspline','SmoothingParam',1);
%     xdata = min(voxel_pos(xslices(v).idx,1)):.1:max(voxel_pos(xslices(v).idx,1));
%     xvals = feval(ft, xdata);
    xvals=interp(xslices(v).vals,up_samp); xdata=interp(voxel_pos(xslices(v).idx,1),up_samp);
    % Y-Slice
%     ft = fit(voxel_pos(yslices(v).idx,2), yslices(v).vals,'smoothingspline','SmoothingParam',1);
%     ydata = min(voxel_pos(yslices(v).idx,2)):.1:max(voxel_pos(yslices(v).idx,2));
%     yvals = feval(ft, ydata);
    yvals=interp(yslices(v).vals,up_samp); ydata=interp(voxel_pos(yslices(v).idx,2),up_samp);
    % Z-Slice
%     ft = fit(voxel_pos(zslices(v).idx,3), zslices(v).vals,'smoothingspline','SmoothingParam',1);
%     zdata = min(voxel_pos(zslices(v).idx,3)):.1:max(voxel_pos(zslices(v).idx,3));
%     zvals = feval(ft, zdata);
     zvals=interp(zslices(v).vals,up_samp); zdata=interp(voxel_pos(zslices(v).idx,3),up_samp);
   
    %% finding half-max edges for multiple peaks on same plane
    % X-Slice
    hmax_vals = xslices(v).vals(peak_xidx)/2;
    tvals = xvals; tvals(tvals<=hmax_vals) = hmax_vals*4;   % 
    [~,w]=findpeaks(-tvals); 
    for s=1:2:length(w)     % finding peak that is sandwiched between half-max edges
        if xdata(w(s)) < voxel_pos(xslices(v).idx(peak_xidx),1) && xdata(w(s+1)) > voxel_pos(xslices(v).idx(peak_xidx),1)
            wx = s;
        end
    end
    half_width(v).edges_x = [xdata(w(wx)) xdata(w(wx+1))]; 
    half_width(v).val_x = range(half_width(v).edges_x);
    
    % Y-Slice
    hmax_vals = yslices(v).vals(peak_zidx)/2;
    tvals = yvals; tvals(tvals<=hmax_vals) = hmax_vals*4;   % 
    [~,w]=findpeaks(-tvals); 
    for s=1:2:length(w)     % finding peak that is sandwiched between half-max edges
        if ydata(w(s)) < voxel_pos(yslices(v).idx(peak_yidx),2) && ydata(w(s+1)) > voxel_pos(yslices(v).idx(peak_yidx),2)
            wx = s;
        end
    end
    half_width(v).edges_y = [ydata(w(wx)) ydata(w(wx+1))]; 
    half_width(v).val_y = range(half_width(v).edges_y);
    
    % Z-Slice
      hmax_vals = zslices(v).vals(peak_zidx)/2;
    tvals = zvals; tvals(tvals<=hmax_vals) = hmax_vals*4;   % 
    [~,w]=findpeaks(-tvals); 
    for s=1:2:length(w)     % finding peak that is sandwiched between half-max edges
        if zdata(w(s)) < voxel_pos(zslices(v).idx(peak_zidx),3) && zdata(w(s+1)) > voxel_pos(zslices(v).idx(peak_zidx),3)
            wx = s;
        end
    end
    half_width(v).edges_z = [zdata(w(wx)) zdata(w(wx+1))]; 
    half_width(v).val_z = range(half_width(v).edges_z);
    
    %% average source leakage
    half_width(v).average = nanmean([half_width(v).val_x half_width(v).val_y half_width(v).val_z]);
end
