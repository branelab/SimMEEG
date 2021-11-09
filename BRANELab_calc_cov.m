function [R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(data,act_samps, ctrl_samps)
%function [R,N,Rbar,Nbar,Rinv,Ninv]=BRANELab_calc_cov(data,act_samps, ctrl_samps)
% INPUT:
%   data = post-stimulus singla + noise matrix (samples x channels x trials)
%   act_samps = samples in active interval used for calculating signal-to-noise ratio.
%   ctrl_samps = samples in control interval used for calculating signal-to-noise ratio.
%
% OUTPUT:
%   R = covariance matrix for active interval to be used in "BRANElab_LCMV_Scalar_beamformer.m"
%   N = covariance matrix for control interval to be used in "BRANElab_LCMV_Scalar_beamformer.m"
%   Rinv = inverted covariance matrix for active interval to be used in "BRANElab_LCMV_Scalar_beamformer.m"
%   Ninv = inverted covariance matrix for control interval to be used in "BRANElab_LCMV_Scalar_beamformer.m"
%   Rbar = covariance matrix for active interval for multi-source event-related beamformer covariance Rbar =((b-mean(b))*(b-mean(b))'). see Moiseev et al. 2011 NeuroImage.
%   Nbar = covariance matrix for control interval for multi-source event-related beamformer covariance Rbar =((b-mean(b))*(b-mean(b))'). see Moiseev et al. 2011 NeuroImage.
%
% written by A. Herdman, UBC June 13, 2015

%% NOTES
%   tikhonov = Tikhonov regularization as a percent (default = 10%). NOTE: This is only used if Rinv or Ninv are badly scaled (not invertible). ;
%   Using pinv(R) instead of svd(R) and regularization

%%
%if nargin<4; tikhonov=10; end
RN_type = {'Active Interval' 'Control Interval'};
RN_name={'Rinv' 'Ninv'};
act_samps(act_samps==0)=1;
for c=1:2;  % loop trhough for active and control intervals
% first baselining data to mean for each trial;
    if c==1;
        xdata = data(act_samps,:,:) - repmat( nanmean(data(act_samps,:,:),1), [length(act_samps) 1 1]);
    elseif c==2;
        xdata = data(ctrl_samps,:,:) - repmat( nanmean(data(ctrl_samps,:,:),1), [length(ctrl_samps) 1 1]);
    end
    [num_samps,num_chans,num_trials]=size(xdata);
    
%    fprintf('calculate covariance then average across trials:  %s\n',RN_type{c});
    cov_dat=nan(num_chans,num_chans,num_trials);
    for t=1:num_trials;
        F=squeeze(xdata(:,:,t))';
        RN_data(:,:,t)=cov(F');
        %RN_data(:,:,t)=(F*F')./ ((ones(size(F,1),size(F,1))*num_samps)-1);     %This is equivalent to cov(F') in line above.
    end
    cov_dat=squeeze(nanmean(cov_dat,3));    % averaging across trials
    RN_data=squeeze(nanmean(RN_data,3));    % averaging across trials
    
    % calculating Rbar for multi-source event-related beamformer covariance Rbar =((b-mean(b))*(b-mean(b))'). see Moiseev et al. 2011 NeuroImage.
%    fprintf('average across trials then calculate covariance of average:  %s\n',RN_type{c});
    F=squeeze(nanmean(xdata,3));    % averaging across trials to get evoked data then calculating covariance
    Rbar_data=cov(F,1);
% Code below gives almost exactly same result as above Rbar_data=cov(F,1). difference is rounding errors close to the floating-poinjt precision of the system. 
%     for s=1:size(F,1); Rbar_data2(s,:,:) = ( F(s,:) -mean(F))'*(F(s,:)-mean(F)); end;
%     Rbar_data2=squeeze(nanmean(Rbar_data2));
    
    %    if c==1; Rinv=cov_dat; R=RN_data; elseif c==2; Ninv=cov_dat; N=RN_data; end
    if c==1;       % for act_int
        R=RN_data;
        Rinv=inv(RN_data); % covariance matrix should never be near singular for BRANELAb's LCMV beamformers - if it is then there is a problem.
        Rbar = Rbar_data;
    if rcond(Rinv)<=eps; 
        fprintf('\nWARNING: Inverted Matrix "%s" is badly scaled and is near singular to working precision!\nPlease view the covariance matrices using surf(%s)\n',RN_name{c},RN_name{c});
        Rinv=pinv(R);   % pseudo-inverse of covariance (similar to what is done in NutMEG
%         fprintf('Inverting R matrix using Tikhonov regularization: %.f%%\n',tikhonov);
%         % code from Brainstorm3 for inverting a matrix with regularization. NOT IMPLEMENTED because matrix should be invertible.
%         [U,S] = svd(R); % Covariance = 1/(m-1) * F * F'
%         S = diag(S); % Covariance = Cm = U*S*U'
%         Si = diag(1 ./ (S + S(1) * tikhonov / 100)); % 1/(S^2 + lambda I)
%         Rinv = U*Si*U'; % U * 1/(S^2 + lambda I) * U' = Cm^(-1)
      end
    elseif c==2;   % for ctrl_int
        N=RN_data;
        Ninv=inv(RN_data);  % covariance matrix should never be near singular for BRANELAb's LCMV beamformers - if it is then there is a problem.
        Nbar = Rbar_data;
   if rcond(Ninv)<=eps; 
        fprintf('\nWARNING: Inverted Matrix "%s" is badly scaled and is near singular to working precision!\nPlease view the covariance matrices using surf(%s)\n',RN_name{c},RN_name{c});
        Ninv=pinv(N);   % pseudo-inverse of covariance (similar to what is done in NutMEG
%         fprintf('Inverting N matrix using Tikhonov regularization: %.f%%\n',tikhonov);
%         % code from Brainstorm3 for inverting a matrix with regularization. NOT IMPLEMENTED because matrix should be invertible.
%         [U,S] = svd(N); % Covariance = 1/(m-1) * F * F'
%         S = diag(S); % Covariance = Cm = U*S*U'
%         Si = diag(1 ./ (S + S(1) * tikhonov / 100)); % 1/(S^2 + lambda I)
%         Ninv = U*Si*U'; % U * 1/(S^2 + lambda I) * U' = Cm^(-1)
    end
     end
end


