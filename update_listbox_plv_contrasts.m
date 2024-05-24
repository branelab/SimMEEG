function update_listbox_plv_contrasts(varargin)
global h
peak_clr = [1 .6 .0]; % true color
true_clr = [0 0 0]; % true color
seed_clr = [0 .6 0]; % seed colors for plv
peak_idx = [];
true_idx = [];

switch h.menu_inv_tfr_type.String{h.menu_inv_tfr_type.Value}
    case {'Total Power' 'Evoked Power' 'Induced Power'}
        if ~isempty(h.current_3D_plv_contrasts_listbox_order)
            h.peak_clr = bsxfun(@mtimes,ones(length(h.current_3D_plv_contrasts_listbox_order),3),peak_clr);  %lines(length(h4.current_3D_peak_idx));
            h.true_clr = bsxfun(@mtimes,ones(length(h.current_3D_plv_contrasts_listbox_order),3),true_clr);  %lines(length(h4.current_3D_peak_idx));
            h.peak_clr(1:size(h.cfg.source.src_clr,1),:) = h.cfg.source.src_clr;
            h.true_clr(1:size(h.cfg.source.src_clr,1),:) = h.cfg.source.src_clr;
            
            peak_idx = h.inv_soln(h.current_inv_soln).plv_comp_idx;
            true_idx = h.sim_data.cfg.source.vx_idx;
        end
    case {'PLV' 'PLI' 'dPLI'}
        h.peak_clr = bsxfun(@mtimes,ones(length(h.current_3D_plv_contrasts_listbox_order),3),peak_clr);  %lines(length(h4.current_3D_peak_idx));
        h.true_clr = bsxfun(@mtimes,ones(length(h.current_3D_plv_contrasts_listbox_order),3),true_clr);  %lines(length(h4.current_3D_peak_idx));
        
        %         peak_clr = bsxfun(@mtimes,ones(length(h.current_3D_plv_contrasts_listbox_order),3),[.9 .6 .3]*.5);  %lines(length(h4.current_3D_peak_idx));
        % %         peak_clr(1:3,:) = h.plv_clr;
        %         peak_clr(1:size(h.cfg.source.src_clr,1),:) = h.cfg.source.src_clr;
        %      try
        peak_idx = h.inv_soln(h.current_inv_soln).plv_contrasts(h.current_3D_plv_contrasts_listbox_order,:);
        %          catch; peak_idx = []; end
        true_idx = h.sim_data.cfg.source.plv_contrast_idx;
        true_idx2 = h.sim_data.cfg.source.plv_contrasts;
        h.peak_clr(1:size(true_idx,1),:) = repmat(seed_clr,size(true_idx,1),1);   % setting hits with color of true sources
        
end
h.listbox_plv_contrasts_txt.ForegroundColor = peak_clr(end,:);
h.listbox_true_plv_contrasts_txt.ForegroundColor = true_clr(end,:);



if max(h.listbox_true_plv_contrasts.Value) > length(true_idx); h.listbox_true_plv_contrasts.Value=1; end
if max(h.listbox_plv_contrasts.Value) > length(peak_idx); h.listbox_plv_contrasts.Value=1; end

%% Peak FC contrasts
h.listbox_plv_contrasts.UserData.clr_str=''; h.listbox_plv_contrasts.String=''; %h.listbox_plv_contrasts.Value = 1;
pre = '<HTML><FONT color="'; post = '</FONT></HTML>';
for v=1:length(peak_idx)
    if any(isnan(peak_idx(v)))
        if size(peak_idx,1)==1  % TFR plots
            peak_name = 'Miss';
        else
            peak_name = sprintf( 'Miss (%.f vs %.f)', peak_idx(v,1), peak_idx(v,2) );
        end
%         clr_str = reshape( dec2hex( round([0 0 0]*255), 2 )',1, 6);
        clr_str = reshape( dec2hex( round(h.peak_clr(v,:)*255), 2 )',1, 6);
    else
        if size(peak_idx,2)==2
            if v<=length(true_idx)  % print out corresponding true source # in parentheses to match up with hit voxel #
                peak_name = sprintf( '(%.f) %.f vs (%.f) %.f', true_idx(v,1), peak_idx(v,1), true_idx(v,2), peak_idx(v,2) );
            else
                peak_name = sprintf('%.f vs %.f',peak_idx(v,1),peak_idx(v,2));
            end
            
        else
            if v<=length(true_idx)  % print out corresponding true source # in parentheses to match up with hit voxel #
                peak_name = sprintf('(%.f) %.f ', v, peak_idx(v));
            else
                peak_name = sprintf('%.f',peak_idx(v));
            end
        end
        
        clr_str = reshape( dec2hex( round(h.peak_clr(v,:)*255), 2 )',1, 6);
    end
    h.listbox_plv_contrasts.UserData.clr_str{v} = [pre clr_str '">' peak_name post];
end
h.listbox_plv_contrasts.String = h.listbox_plv_contrasts.UserData.clr_str;


%% True FC contrasts
h.listbox_true_plv_contrasts.UserData.clr_str=''; h.listbox_true_plv_contrasts.String=''; %h.listbox_true_plv_contrasts.Value = 1;
pre = '<HTML><FONT color="'; post = '</FONT></HTML>';
for v=1:length(true_idx)
    if isnan(true_idx(v))
        peak_name = 'Miss';
    else
        if size(true_idx,2)==2
            peak_name = sprintf( '(%.f) %.f vs (%.f) %.f', true_idx(v,1), true_idx2(v,1), true_idx(v,2), true_idx2(v,2) );
            peak_name = sprintf('%.f vs %.f',true_idx(v,1),true_idx(v,2));
        else
            peak_name = sprintf('(%.f) %.f',v,true_idx(v));
        end
        
    end
    clr_str = reshape( dec2hex( round(h.true_clr(v,:)*255), 2 )',1, 6);
    h.listbox_true_plv_contrasts.UserData.clr_str{v} = [pre clr_str '">' peak_name post];
end
h.listbox_true_plv_contrasts.String = h.listbox_true_plv_contrasts.UserData.clr_str;
if h.listbox_plv_contrasts.Value > length(h.current_3D_plv_contrasts); h.listbox_plv_contrasts.Value = 1; end
if h.listbox_true_plv_contrasts.Value > length(h.listbox_true_plv_contrasts.String);  h.listbox_true_plv_contrasts.Value = 1; end





