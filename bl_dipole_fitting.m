function [fit_error, dips, fit_params, plot_params] = bl_dipole_fitting(data, dips, leadfield, hdm_tripos, sens, bst_hdm, fit_params, plot_params, plot_flag, update_flag)
% function [fit_error, dip, fit_params, plot_fit_params] = bl_dipole_fitting(data, leadfield, fit_params, plot_fit_params, plot_flag)
%
% INPUT
%   data         = data to be fitted [samps x chans]
%   hdm_tripos   = headmodel.tri and headmodel.pos used to create "leadfield.pos"
%   sens         = sensor struct as defined for FieldTrip but also same as for BRANELab format
%   bst_hdm       = spherical headmodel with struct defined in Brainstorm
%   dips         = "dipole" class created using "dipoles.m"
%   fit_params   = "dipfit_params" class created using "dipfit_params.m"
%   plot_params  = handles to axes after plotting returned from "bl_dipdit_error.m", set to plot_fit_params=[]; to start new plot
%   plot_flag    = (0) no plots (1) plot iterations (2) plot final
%   update_flag  = (1) update leadfields but don't fit (0) perform nonlinear optimization fitting
%
% Based on equations from Scherg (1990) in Auditory Evoked Magnetic Fields and Electric Potentials vol 6 pp40-69
%
% U = C*S;  where U = [chans x samps]; C = forwards model [chans x sources]; S=source waves [sources x samps]
% S = C-1*S;  where C-1=inverse(C);
%

% xt = {'fminsearch' 'fmincon' 'fminmax' 'Genetic' 'GeneticMulti' 'ParticleSwarm' 'Pareto' 'Pattern' 'Surrogate' 'SimulatedAnnealing'};

%% set params for fitting dips
% clear params
%% dip params

params.dip_pos = reshape([dips.dip_pos], 3, size(dips,2))';
params.dip_ori = reshape([dips.dip_ori], 2, size(dips,2))';
params.dip_idx = find_nearest_voxel(params.dip_pos,leadfield.pos);
params.fit_on = [dips.dip_on];
params.fit_order = [dips.fit_order];
params.fit_ori = [dips.fit_ori];
params.fit_idx = [dips.fit_on];
params.sym_idx = [dips.sym_idx];
params.sym_dist = [dips.sym_dist];
params.vect_idx = [dips.vect_idx];
params.search_pos = reshape([dips.search_pos], 3, size(dips,2))';
params.search_radius = [dips.search_radius]; params.search_radius(params.search_radius<0) = nan;
params.fit_order = params.fit_order.*params.fit_idx;

%% fit_params
params.optimType = fit_params.optimType;
params.maxIter = fit_params.maxIter;
params.maxTime = fit_params.maxTime;
params.verbose = fit_params.verbose;
params.fitType = fit_params.fitType; 
params.parallel_compute = logical(fit_params.parallel_compute); % for parallel computing
params.residual_var = 1;
%% Field trip params needed for finding leadfields on the fly - currently hard-coding opts.
params.leadfield = leadfield;
params.hdm_tripos = hdm_tripos;
if ~isfield(bst_hdm,'o') % not a sphereical head model so creating one
    for s=1:3
        [bst_hdm.o(s,:), bst_hdm.r(s)] = fitsphere(bst_hdm.bnd(s).pos);
        % hold on; circle(bst_hdm.o(s,:), bst_hdm.r(s),1000);
    end
    [bst_hdm.r,sidx] = sort(bst_hdm.r);
    bst_hdm.o = bst_hdm.o(sidx(1),:); % selecting innermost sphere
    bst_hdm.c = [0.3300    0.0042    0.3300]; % default 3-sphere conductivities brain, skull, scalp
end
% disp(bst_hdm)
params.bst_hdm = bst_hdm;
params.sens = sens;
params.ft_opts = {'reducerank', [3], 'backproject', [], 'normalize', [], 'normalizeparam', [], 'weight', []};


%% fit only those selected in fit_params.fit_idx
% Perform fitting
if update_flag==1  %% Update not fitting
    %% set dipoles for fitting
    fit_idx = params.fit_order>0;
    params.fit_on(fit_idx) = true;  % turns on if was off
    params.fit_idx = fit_idx;

    switch params.fitType
        case 'LeadField'
            dip_idx_ori = [params.dip_idx(fit_idx)' params.dip_ori(fit_idx,:)];
            rshape_fact = 3; % for reshaping for ParticleSwarm, Genetic, ... 
        case 'Continuous'
            dip_idx_ori = [params.dip_pos(fit_idx,:) params.dip_ori(fit_idx,:)];
            rshape_fact = 5; % for reshaping for ParticleSwarm, Genetic, ... 
    end

else      %% fit
    num_fits = max([dips.fit_order]);
    for t=1:num_fits

        %% set dipoles for fitting
        fit_idx = params.fit_order==t;
        params.fit_on(fit_idx) = true;  % turns on if was off
        params.fit_idx = fit_idx;

        switch params.fitType
            case 'LeadField'
                dip_idx_ori = [params.dip_idx(fit_idx)' params.dip_ori(fit_idx,:)];
                num_dip = size(dip_idx_ori,1);
                lowbnd = repmat([1 0 0],num_dip,1); % lower bound for dip_pos(1:3) and ori(ax,el) in degrees
                upbnd = repmat([size(params.leadfield.pos,1) 360 360],num_dip,1); % lower bound for dip_pos(1:3) and ori(ax,el) in degrees
                rshape_fact = 3; % for reshaping for ParticleSwarm, Genetic, ...
            case 'Continuous'
                dip_idx_ori = [params.dip_pos(fit_idx,:) params.dip_ori(fit_idx,:)];
                num_dip = size(dip_idx_ori,1);
                lowbnd = repmat([min(params.hdm_tripos.pos) 0 0],num_dip,1); % lower bound for dip_pos(1:3) and ori(ax,el) in degrees
                upbnd = repmat([max(params.hdm_tripos.pos) 360 360],num_dip,1); % lower bound for dip_pos(1:3) and ori(ax,el) in degrees
                rshape_fact = 5; % for reshaping for ParticleSwarm, Genetic, ...
        end


        %% Perform initial fitting
        % Plot initial dip
        % plot_params = []; plot_flag=1;
        [fit_error, dip, plot_params] = bl_dipfit_error(dip_idx_ori, data, params, plot_params, plot_flag);
        % plot_flag=0;

        %% Optimization Fcn and params
        objcFcn = @(dip_idx_ori)bl_dipfit_error(dip_idx_ori, data, params, plot_params, 0);
        tic;
        output = '';

        %% optimization
        switch fit_params.optimType
            case 'fminsearch'
                % [x,fval] = fminsearch(f,dip_idx_ori);
                options = optimset(...
                    'MaxIter',fit_params.maxIter,...
                    'MaxFunEvals',2*fit_params.maxIter*size(params.dip_pos,1),...
                    'Display', fit_params.verbose,"MaxTime",params.maxTime);
                [dip_idx_ori, fval, exitflag, output] = fminsearch(@bl_dipfit_error, dip_idx_ori, options, data, params, plot_params, plot_flag);
            case 'fmincon'
                options = optimoptions("fmincon","MaxFunctionEvaluations",2*fit_params.maxIter,"MaxIterations",fit_params.maxIter,"Display",fit_params.verbose,"UseParallel",params.parallel_compute);
                [dip_idx_ori, fval, xflag, lambda,grad,hessian] = fmincon(objcFcn, dip_idx_ori,[],[],[],[], lowbnd, upbnd, [],options);
            case 'fminmax'
                options = optimoptions("fminimax","MaxFunctionEvaluations",2*fit_params.maxIter,"Display",fit_params.verbose,"UseParallel",params.parallel_compute);
                [dip_idx_ori, fval, exitflag, output] = fminimax(objcFcn,dip_idx_ori,[],[],[],[], lowbnd, upbnd,[],options);
            case 'Genetic'
                objcFcn = @(dip_idx_ori)bl_dipfit_error_ga(dip_idx_ori, data, params, plot_params, plot_flag);
                nvars = size(dip_idx_ori(:)',2);
                options = optimoptions("ga","Display",fit_params.verbose,"UseParallel",params.parallel_compute,"MaxTime",params.maxTime);
                [dip_idx_ori, fval, xflag,output,pop,scores] = ga(objcFcn, nvars,[],[],[],[], lowbnd(:)', upbnd(:)', [],[],options);
                dip_idx_ori = reshape(dip_idx_ori,num_dip, rshape_fact);
            case 'GeneticMulti'
                objcFcn = @(dip_idx_ori)bl_dipfit_error_ga(dip_idx_ori, data, params, plot_params, plot_flag);
                nvars = size(dip_idx_ori(:)',2);
                options = optimoptions("gamultiobj","Display",fit_params.verbose,"UseParallel",params.parallel_compute,"MaxTime",params.maxTime);
                [dip_idx_ori, fval, xflag,output,population,scores] = gamultiobj(objcFcn, nvars,[],[],[],[], lowbnd(:)', upbnd(:)', [],[],options);
                dip_idx_ori = reshape(dip_idx_ori,num_dip, rshape_fact);
            case 'ParticleSwarm'
                objcFcn = @(dip_idx_ori)bl_dipfit_error_ga(dip_idx_ori, data, params, plot_params, plot_flag);
                nvars = size(dip_idx_ori(:)',2);
                options = optimoptions("particleswarm","MaxIterations",fit_params.maxIter,"Display",fit_params.verbose,"UseParallel",params.parallel_compute,"MaxTime",params.maxTime,"SwarmSize",20);
                [dip_idx_ori, fval, xflag, output, pts] = particleswarm(objcFcn, nvars, lowbnd(:)', upbnd(:)', options);
                dip_idx_ori = reshape(dip_idx_ori,num_dip, rshape_fact);
            case 'Pareto'
                objcFcn = @(dip_idx_ori)bl_dipfit_error_ga(dip_idx_ori, data, params, plot_params, plot_flag);
                nvars = size(dip_idx_ori(:)',2);
                options = optimoptions("paretosearch","MaxFunctionEvaluations",2*fit_params.maxIter,"MaxIterations",fit_params.maxIter,"Display",fit_params.verbose,"UseParallel",params.parallel_compute,"MaxTime",params.maxTime);
                [dip_idx_ori, fval, xflag, output] = paretosearch(objcFcn,nvars,[],[],[],[],lowbnd(:)', upbnd(:)', [],options);
                dip_idx_ori = reshape(dip_idx_ori,num_dip, rshape_fact);
            case 'Pattern'
                options  = optimoptions("patternsearch",...
                    "MaxIterations",fit_params.maxIter,"MaxFunctionEvaluations",2*fit_params.maxIter,...
                    "Display",fit_params.verbose,"UseParallel",params.parallel_compute,"MaxTime",params.maxTime,...
                    "FunctionTolerance",fit_params.fcn_tolerance,"StepTolerance",fit_params.step_tolerance,"MeshTolerance",1e-10);
                [dip_idx_ori, fval, xflag, output] = patternsearch(objcFcn,dip_idx_ori,[],[],[],[], lowbnd, upbnd,options);
            case 'Surrogate'
                fprintf('Not Available: %s\n',fit_params.optimType);
                % options = optimoptions("surrogateopt","MaxFunctionEvaluations",200,"Display",fit_params.verbose,"PlotFcn",[],"UseParallel",params.parallel_compute);
                % [dip_idx_ori, fval, xflag, output] = surrogateopt(objcFcn, lowbnd, upbnd,[],[],[],[],[],options);
            case 'SimulatedAnnealing'
                options = optimoptions("simulannealbnd","MaxFunctionEvaluations",fit_params.maxIter,"Display",fit_params.verbose,"MaxTime",params.maxTime);
                [dip_idx_ori, fval, xflag, output] = simulannealbnd(objcFcn,dip_idx_ori, lowbnd, upbnd,options);
        end

        %% setting fitted dips back into params for next order of fits
        switch params.fitType
            case 'LeadField'
                params.dip_pos(fit_idx,:) = leadfield.pos(round(dip_idx_ori(:,1)),:);
                params.dip_ori(fit_idx,:) = dip_idx_ori(:,2:3);
            case 'Continuous'
                params.dip_pos(fit_idx,:) = dip_idx_ori(:,1:3);
                params.dip_ori(fit_idx,:) = dip_idx_ori(:,4:5);
        end


        %% update fit
        [fit_error, dip, plot_params] = bl_dipfit_error(dip_idx_ori, data, params, plot_params, plot_flag);
        if ~isempty(dip)
            params.dip_pos = dip.dip_pos;
        else
            switch params.fitType
                case 'LeadField'
                    dip_idx_ori = [params.dip_idx(fit_idx)' params.dip_ori(fit_idx,:)];
                case 'Continuous'
                    dip_idx_ori = [params.dip_pos(fit_idx,:) params.dip_ori(fit_idx,:)];
            end
            [fit_error, dip, plot_params] = bl_dipfit_error(dip_idx_ori, data, params, plot_params, plot_flag);
        end
        %% Verbose output
        if strcmp(fit_params.verbose,"on")
            disp(output);
            fprintf('Fit Time = %.1f sec\n',toc);
        end

    end
end
%% Final Update Fit
[fit_error, dip, plot_params] = bl_dipfit_error(dip_idx_ori, data, params, plot_params, plot_flag);
if ~isempty(dip)
    params.dip_pos = dip.dip_pos;
     switch params.fitType
            case 'LeadField'
               params.dip_ori(logical(params.fit_idx),:) = dip_idx_ori(:,2:3);
            case 'Continuous'
               params.dip_ori(logical(params.fit_idx),:) = dip_idx_ori(:,4:5);
     end
    params.residual_var = fit_error;

    %% update "dips"
    vidx = find(logical(params.fit_on));
    dims = size(dip.leadfield);

    % lf = reshape(dip.leadfield, [dims(1) 3 dims(2)/3]);
    xidx = find(~logical(params.fit_on));
    for nx=1:length(xidx);   dips(xidx(nx)).leadfield = []; end
    for v=1:length(vidx)
        dips(vidx(v)).dip_pos = params.dip_pos(vidx(v),:);
        dips(vidx(v)).dip_ori = params.dip_ori(vidx(v),:);
        dips(vidx(v)).leadfield = dip.leadfield(:,dip.lf_idx==vidx(v));
        dips(vidx(v)).wts = dip.wts(:,dip.lf_idx==vidx(v));
    end

else

    params.dip_pos = leadfield.pos(params.dip_idx,:);
    switch params.fitType
        case 'LeadField'
            dip_idx_ori = [params.dip_idx(fit_idx)' params.dip_ori(fit_idx,:)];
        case 'Continuous'
            dip_idx_ori = [params.dip_pos(fit_idx,:) params.dip_ori(fit_idx,:)];
    end
    [fit_error, dip, plot_params] = bl_dipfit_error(dip_idx_ori, data, params, plot_params, plot_flag);
    vidx = find(logical(params.fit_on));
    for v=vidx
        dips(v).dip_pos = params.search_pos(v,:);
        dips(v).leadfield = [];
    end

end


