function [fit_error, dip, plot_params] = bl_dipfit_error_ga(dip_idx_ori, data, params, plot_params, plot_flag)
% function [dip, fit_error] = bl_dipfit_error(data, leadfield, dip_idx)
% INPUT
%   dip_idx = indices of leadfield to fit
%   data = [samps x chans];
%   params = struct('dip_pos',[-50 20 60; 0 10 57; 0 -10 57],'spatiotemp',true,'fit_on',[1 1 1],'fit_order',[1 1 1],'params.fit_idx',[1 1 1],'fit_ori',[1 1 1],'sym_idx',[0 0 0; 0 0 2; 0 0 0], 'vect_idx',[1 1 1], 'optimType','fminsearch','maxIter',1000,'verbose','off');

%% init vars
fit_error = 1;
dip = '';
dip_on = find(logical(params.fit_on)==1);
switch params.fitType
    case 'LeadField'
        %% get optimizer dip_idx
        dims = size(dip_idx_ori);
        dip_idx_ori = reshape(dip_idx_ori,dims(2)/3,3);
        didx = round(dip_idx_ori(:,1));
        %% get positions for didx
        params.dip_pos(params.fit_idx,:) = params.leadfield.pos(didx,:);
        %% set all dip_pos for all dips on
        params.dip_pos = params.dip_pos(dip_on,:);
        %% get dip_idx (lf2_idx) for all dips on
        dip_idx = find_nearest_voxel(params.dip_pos(dip_on,:), params.leadfield.pos);
        %% set dip_position to leadfield position
        params.dip_pos = params.leadfield.pos(dip_idx,:);
        params.dip_ori(params.fit_idx,:) = dip_idx_ori(:,2:3);

    case 'Continuous'
        num_dips = size(dip_idx_ori,2)/5;
        dip_idx_ori = reshape(dip_idx_ori,num_dips,5);
        params.dip_pos(params.fit_idx,:) = dip_idx_ori(:,1:3);
        params.dip_pos = params.dip_pos(dip_on,:);
        params.dip_ori(params.fit_idx,:) = dip_idx_ori(:,4:5);

if all( all( (params.dip_pos(:,1)>=min(params.hdm_tripos.pos(:,1)) & params.dip_pos(:,1)<=max(params.hdm_tripos.pos(:,1))) & ...
        (params.dip_pos(:,2)>=min(params.hdm_tripos.pos(:,2)) & params.dip_pos(:,2)<=max(params.hdm_tripos.pos(:,2))) &  ...
        (params.dip_pos(:,3)>=min(params.hdm_tripos.pos(:,3)) & params.dip_pos(:,3)<=max(params.hdm_tripos.pos(:,3))) & ...
        (params.dip_ori>=0 & params.dip_ori<=360 ) ))
   else
    return
end



end


%% set orientations
params.dip_ori = params.dip_ori(dip_on,:);
dip_ori = deg2rad(params.dip_ori);    % az, el
[x,y,z] = sph2cart(dip_ori(:,1), dip_ori(:,2), ones(size(dip_ori,1),1));
dip_ori = [x y z];

%% Next
num_dips = size(params.dip_pos,1);

%% Symmetrical constraints
sym_flag = true;
params.sym_idx = params.sym_idx(dip_on);
for v=1:size(params.sym_idx,2)   % rows=XYZ; clms=dips;
    if any(params.sym_idx(v)) && params.fit_idx(v)==1
        if params.sym_idx(v)~=0
            params.dip_pos(v,2) = -params.dip_pos(params.sym_idx(v),2);  % hemi symmetry
            params.dip_pos(v,[1 3]) = params.dip_pos(params.sym_idx(v),[1 3]);
            %% sym_dist criterion
            if params.sym_dist(v)>0 && ( abs(params.dip_pos(v,2) - params.dip_pos(params.sym_idx(v),2)) < params.sym_dist(v))
                sym_flag = false;
            end
        end
    end
end

%% check search radius
search_flag=true;
dip_pos = params.dip_pos(params.fit_idx,:);
search_pos = params.search_pos(params.fit_idx,:);
diff_pos = dip_pos-search_pos;
if any(any(abs(diff_pos) > params.search_radius(params.fit_idx)'))
    search_flag=false;
end

%% Run dipfit
if ~search_flag || ~sym_flag    % dip_pos outside of search_radius --> do not perform dipfit
   if ~strcmp(params.verbose,'off'); fprintf('Outside search radius\n'); end
else
    switch params.fitType
        case 'LeadField' %% Use pre-computed leadfields
            lf2 = params.leadfield.H(:,:,round(dip_idx));
            lf2 = lf2-mean(lf2); % average referencing leadfield;
        case 'Continuous'
            %% FieldTrip's leadfield calculator - but need to convert by moments*1e-3 to be close to BESA but they return diffreent lf compared to "bst_eeg_sph.m" - not sure why?
            % lf = ft_compute_leadfield(params.dip_pos, params.sens, params.ft_hdm, params.ft_opts{:});
            % These also work --> Brainstorm leadfield calculator - but need to convert by moments to be consistent with BESA
            % cidx = [1 2 3];
            % lf = bst_eeg_sph(params.dip_pos*1e-3, params.sens.chanpos*1e-3, params.ft_hdm.o*1e-3, params.ft_hdm.r(cidx)*1e-3, params.ft_hdm.cond(cidx));
            % center = [0.0030   -0.0009    0.0509]; R = [0.0855    0.0904    0.0972]; sigma = [0.3300    0.0042    0.3300];
            lf2 = bst_eeg_sph(params.dip_pos*1e-3, params.sens.chanpos*1e-3, params.bst_hdm.o*1e-3, params.bst_hdm.r*1e-3, params.bst_hdm.c);
            if isempty(lf2)
                if ~strcmp(params.verbose,'off')
                    fprintf('Dipole(s) outside of head model\n');
                end
                return
            end
            lf2 = lf2-mean(lf2); % average referencing leadfield;
            lf2 = lf2*1e-3;
            dims = size(lf2);
            lf2 = reshape(lf2, [dims(1) 3 dims(2)/3]);
    end
    %% set lf based on ori
    % lf2 = reshape(lf,size(lf,1),3,num_dips);  % struct = [chan x vect x dipole]
    lf = []; lf_idx = [];
    for vx=1:length(dip_on)
        v=dip_on(vx);
        if params.vect_idx(v)==1  % vector dipole so keep lfs
            lf = cat(2, lf, squeeze(lf2(:,:,v)) );
            lf_idx = [lf_idx v v v];
        else
            lf3 = squeeze(lf2(:,:,v))*dip_ori(v,:)';
            lf = cat(2, lf, lf3);
            lf_idx = [lf_idx v];
        end
    end
    % lfx = lf; load('lf2'); lf = lf(:,2); 
    % lf = lf/.9974; 
    % lf = lf-mean(lf); % average referencing leadfield;

    moments = pinv(lf)*data'; % getting source wave moments;
    fit_data = lf*moments;  % fitted sensor data
    dif_data = data-fit_data';

    %% Residual Variance
    fit_error = double(sum(dif_data(:).^2) / sum(data(:).^2));

    %% Output dip.
    dip.dip_pos = params.dip_pos;
    dip.leadfield = lf;
    dip.wts = pinv(lf)';
    % dip.wts = lf;
    dip.moments = moments;
    dip.fit_data = fit_data;
    dip.resid_var = fit_error;
    dip.lf_idx = lf_idx; 


    if plot_flag==1

        if ~isfield(plot_params,'dip_scatter'); plot_params.dip_scatter =[]; end
        if isempty(plot_params.dip_scatter)
            plot_params.hfig = figure;
            vw_angle = [-90 90; 180 0];
            for vw=1:2
                plot_params.haxes_anat(vw) = subplot(2,2,vw); cla; hold on;
                plot_params.lf_scatter(vw) = scatter3(params.leadfield.pos(:,1), params.leadfield.pos(:,2), params.leadfield.pos(:,3), 'ko','MarkerEdgeColor','none', 'MarkerFaceColor','k','MarkerFaceAlpha',.05,'SizeData', 2);
                plot_params.dip_scatter(vw) = scatter3(params.dip_pos(:,1), params.dip_pos(:,2), params.dip_pos(:,3), 'ro','filled');
                % plot_params.hdm_patch(vw) =patch(plot_params.haxes_anat(vw), 'Faces', params.hdm_tripos.tri,'Vertices', params.leadfield.pos,'FaceColor','k','FaceAlpha',.1,'EdgeColor','k','EdgeAlpha',.15);

                %% Plot Orientations
                ori_length = 25;
                for v=1:num_dips
                    if params.vect_idx(v)==0  % scalar dipole so plot 1-dipole plane
                         plot_params.dip_ori(v) = plot3(plot_params.haxes_anat(vw), params.dip_pos(v,1)+[0 dip_ori(v,1)*ori_length], params.dip_pos(v,2)+[0 dip_ori(v,2)*ori_length], params.dip_pos(v,3)+[0 dip_ori(v,3)*ori_length],'r-', 'LineWidth',2);
                    else % vector dipole so plot 3-dipole vectors
                        plot_params.dip_ori3(v,:) = plot3(plot_params.haxes_anat(vw), params.dip_pos(v,1)+[0 1; 0 0; 0 0]'*ori_length, params.dip_pos(v,2)+[0 0; 0 1; 0 0]'*ori_length, params.dip_pos(v,3)+[0 0; 0 0; 0 1]'*ori_length,'r-', 'LineWidth',2);
                    end
                end
                axis off; axis tight;
                view(plot_params.haxes_anat(vw),vw_angle(vw,:));
            end

            plot_params.haxes_waves = subplot(2,1,2); cla; hold on;
            plot_params.rec_data = plot(squeeze(data),'k');        % recorded data
            plot_params.fit_data = plot(squeeze(fit_data'),'g');   % fitted data
            plot_params.dif_data = plot(squeeze(dif_data),'r');    % residual data
            title(plot_params.haxes_waves, sprintf('RV = %.1f%%\n',fit_error*100))

        else    % Update dip locations and waveforms
            for vw=1:2
                set(plot_params.dip_scatter, 'XData',params.dip_pos(:,1))
                set(plot_params.dip_scatter, 'YData',params.dip_pos(:,2))
                set(plot_params.dip_scatter, 'ZData',params.dip_pos(:,3))
               %% Plot Orientations
                ori_length = 25;
                for v=1:num_dips
                    if params.vect_idx(v)==0  % scalar dipole so plot 1-dipole plane
                         plot_params.dip_ori(v).XData = params.dip_pos(v,1)+[0 dip_ori(v,1)*ori_length];
                         plot_params.dip_ori(v).YData = params.dip_pos(v,2)+[0 dip_ori(v,2)*ori_length];
                         plot_params.dip_ori(v).ZData = params.dip_pos(v,3)+[0 dip_ori(v,3)*ori_length];
                    else % vector dipole so plot 3-dipole vectors
                        plot_params.dip_ori3(v,:) = plot3(plot_params.haxes_anat(vw), params.dip_pos(v,1)+[0 1; 0 0; 0 0]'*ori_length, params.dip_pos(v,2)+[0 0; 0 1; 0 0]'*ori_length, params.dip_pos(v,3)+[0 0; 0 0; 0 1]'*ori_length,'r-', 'LineWidth',2);
                    end
                end
            end

            plot_params.haxes_waves = subplot(2,1,2); cla; hold on;
            plot_params.rec_data = plot(squeeze(data),'k');        % recorded data
            plot_params.fit_data = plot(squeeze(fit_data'),'g');   % fitted data
            plot_params.dif_data = plot(squeeze(dif_data),'r');    % residual data
            title(plot_params.haxes_waves, sprintf('RV = %.1f%%\n',fit_error*100))
        end
        drawnow

    end
end
if ~strcmp(params.verbose,'off')
    fprintf('RV = %.2f%%\n',fit_error*100);
end
