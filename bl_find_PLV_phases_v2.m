function [phase_dist,PLV_trials_est, PLV_evoked_est,PLI_trials_est,dPLI_trials_est]=bl_find_PLV_phases_v2(PLV_targets,PLV_evoked,PLI_targets,phase_lag,num_trials,plv_thresh,max_perm)
%function [phase_dist,PLV_trials_est, PLV_evoked_est,PLI_trials_est,dPLI_trials_est]=bl_find_PLV_phases_v2(PLV_targets,PLV_evoked,PLI_targets,num_trials,plv_thresh,max_perm)
%
% This function will simulate starting phases for each 3 signals with the Phase-locking Values (PLVs) for
% across trials within a signal (PLV_evoked) and across trials among signal comparisons sig1-sig2, sig1-sig3, and sig2-sig3 (PLV_targets).
% Found PLV will be stopped when PLVs are within the plv_thresh from the PLV_trial and PLV_evoked targets.
%
%
% INPUT:
%   PLV_targets = PLV target values for nchoosek(N,2) signals (WARNING! N==3).
%   PLV_evoked = PLV target values time-locked/consistent phase within a single signal. (WARNING! N==3).
%   PLI_targets = PLI target values for nchoosek(N,2) signals (WARNING! N==3). negative values = signal 1 leads signal 2; positive values = signal 2 leads signal 1.
%   phase_lag  = phase-lag (radians; relative to sample(1)) of each 3 signals within the signal interval --> (phase_lag/360)*2*pi);    cos(phase_lag) = correlation of signal relative to zero-phase onset
%                   Warning: large phase lag differences among sources can obliterate PLI and PLI_targets will not be found. The larger the PLI targets, the smaller phase_lag differences are needed among sources..
%   num_trials = number of trials in simulation
%   plv_thresh = .01; % minimum difference between target PLV and found PLV for all three comparisons
%   max_perm   = 5:9 % maximum number of permutations to search for PLV_targets (must be <10 or memory will fail on most computers) .
%
% OUTPUT:
%   phase_dist = [N x num_trials] distribution of phase of each signal N in order so estimated PLV will equal PLV_targs.
%   PLV_trials_est = estimated PLV_targets based on "phase_dist" differences across 3 signals .
%   PLV_evoked_est = estimated PLV_evoked based on "phase_dist" within each 3 signals .
%   PLI_trials_est = estimated PLI_targets based on "phase_dist" within each 3 signals .
%   dPLI_trials_est = estimated directed PLI_targets based on "phase_dist" within each 3 signals .
%           Note: dPLI is simply 1/2 of PLI with directionality of lead/lag --> negative values = signal 1 leads signal 2; positive values = signal 2 leads signal 1. .
%
%   written by Anthony Herdman July 24, 2018
%
%   NOTE: The mimization procedures is a time-consuming and is not
%   a mathematically elegant solution. This could be drastically improved
%   but it works for now.

% for testing
% PLV_targets=[0.5 0.5 0.5]; PLV_evoked=[0.3 0.3 0.3]; num_trials=100; plv_thresh=.05; num_iter=1e3;


%% Initializing outputs
phase_dist=[];PLV_trials_est=[]; PLV_evoked_est=[];PLI_trials_est=[];dPLI_trials_est=[];



% Note: if PLV_evoked are equal then error occurs because q1, q2, q3 all
% have equal phases after sorting, thus no PLV or PLI
fprintf('Adjusting PLV_evoked\n');
PLV_evoked_org=PLV_evoked;
while length(unique(PLV_evoked))==1 || sum(find(PLV_evoked>1))>0 % all values of PLV_evoked are the same, but need to adjust to all different to allow program to run; Also PLV can not > 1
    fprintf('.');
    PLV_evoked=PLV_evoked_org+(randn(3,1)*.01);
    fprintf('\b');
end
fprintf('Adjusted PLV_evoked to run program: PLV_evoked = [%.4f %.4f %.4f]\n',PLV_evoked)

if size(PLV_targets,1)>size(PLV_targets,2); PLV_targets=PLV_targets'; end  % reshaping so that PLV_targets are in a row not down a column;
if size(PLV_evoked,1)>size(PLV_evoked,2); PLV_evoked=PLV_evoked'; end  % reshaping so that PLV_evoked are in a row not down a column;
PLI_trials_est=[];
dPLI_trials_est=[];

%% PolyFit of PLV_fit data
sigma_fit=[0:1/1e2:1]; fit_num_trials=1e4;
clear PLV_fit
for t=1:length(sigma_fit)
    pdiff=(randi([-100 100],fit_num_trials,1)/100*sigma_fit(t))*pi;    % distribution of random numbers from -pi to pi radians multiplied by the standard deviation sigma
    PLV_fit(t)=abs((sum(exp(1i*(squeeze( pdiff' ))),2))')./fit_num_trials;
end
[p,s]=polyfit(PLV_fit,sigma_fit,9);
% plv=[0:1/1e4:1]; for t=1:length(plv); sm(t)=polyval(p,plv(t),s); end;
% clf; plot(PLV_fit,sigma_fit,'k.'); hold on; plot(plv,sm,'r'); xlabel('PLV'); ylabel('Std Dev of Phase Diff distribution'); axis([0 1 0 1]);

%% This function converts PLV_targs to stadard deviation of randomized phase differences (radians) based on aove PLV_fit. There is probably an exact mathematical solution but for now this works within a reasonable amount of error in the polyfit .
equal_rand=repmat([0:1/num_trials:1-(1/num_trials)]-0.5,[3 1]);
clear sigma_trial sigma_evk;
for cx=1:length(PLV_targets)
    sigma_trials(cx) = 2*pi*(polyval(p,PLV_targets(cx),s)); % standard deviation pf phases in radians 0:2*pi
    sigma_evk(cx) = 2*pi*(polyval(p,PLV_evoked(cx),s));
end

phase_trials=bsxfun(@times,equal_rand,sigma_trials');
phase_evk=bsxfun(@times,equal_rand,sigma_evk');
PLV_evoked_est=abs((sum(exp(1i*(squeeze( phase_evk ))),2))')./num_trials;
PLV_trials_est=abs((sum(exp(1i*(squeeze( phase_trials ))),2))')./num_trials;

%% sort phase_evk to yield largest PLVs among 3sigs, then reshuffle until PLV_trial.
% Sort phases
[sort_phases,sort_idx]=sort(phase_evk,2);   % if using equal_rand, then phase_evk are already sorted
% q1=sort_phases(1,:);
% q2=sort_phases(2,:);
% q3=sort_phases(3,:);
q1=sort_phases(1,:)+phase_lag(1);
q2=sort_phases(2,:)+phase_lag(2);
q3=sort_phases(3,:)+phase_lag(3);
% qdiff(1,:)=q1-q2; qdiff(2,:)=q1-q3; qdiff(3,:)=q2-q3;

%% Permutations for resuffling after down-sampled num_trials by number of permutations but indices are extended across the size of q1, q2, q3
% max_perm=10; % ~ 300 MBytes of data per permutation

if floor(num_trials/max_perm)~=num_trials/max_perm  % This means that max_perms can NOT divide equally into num_trials
    fprintf('WARNING! For more accurate PLV estimation, Set num_trials(%.f) to be equally divisible by max_perms(%.f). \n',num_trials,max_perm);
end

%% the following code needs work because it is too time consuming to iterate throught permutations.
% Use the old way of approximating
% % new way of randomly selecting possible permutations
% q2_perms=rand_perms(num_trials,max_perm,1);
% q3_perms=rand_perms(num_trials,max_perm,1);
% q1_perms=repmat(1:num_trials,max_perm,1);
%
% PLVtest=ones(1,3); kk=0;
%
% plv12=abs((sum(exp(1i*(squeeze( q1(q1_perms)-q2(q2_perms) ))),2))')./num_trials;
% plv12=abs((sum(exp(1i*(squeeze( q1(q1_perms)-q2(fliplr(q2_perms)) ))),2))')./num_trials;

%% old way of getting PLVs
num_perms=floor(num_trials/max_perm);
q2_perms=perms(1:max_perm)*num_perms; % last set of permutations
q3_perms=perms(1:max_perm)*num_perms; % last set of permutations
q1_perms=repmat(fliplr(1:max_perm),size(q2_perms,1),1)*num_perms; % same order for all permutations

%% search for PLVs
PLVtest=ones(1,3); kk=0;
q1_order=nan(num_perms,max_perm);
q2_order=q1_order;
q3_order=q1_order;
warn_str=[];
% for np=1:num_perms
np=1;
q1idx=q1_perms-(num_perms-np);
q2idx=q2_perms-(num_perms-np);
q3idx=q3_perms-(num_perms-np);

%     while max(abs(PLVtest-PLV_targets))>plv_thresh && kk<num_iter
%         kk=kk+1;
% find plv12, plv13, and plv23 with closest match to PLV_targets
% PLV q1 vs q2

% if sum(PLV_evoked)==0   % added small variance in the phases so that exact subtraction q1-q2, q1-q3, q2-q3 does not occur for the evoked phase matrix
% q1=q1-1e-8;
% end

plv12=abs((sum(exp(1i*(squeeze( q1(q1idx)-q2(q2idx) ))),2))')./max_perm;
plv13=abs((sum(exp(1i*(squeeze( q1(q1idx)-q3(q3idx) ))),2))')./max_perm;
xidx12=plv12<PLV_targets(1)+plv_thresh & plv12>PLV_targets(1)-plv_thresh; % these are all pssible permutations that will yield a plv12 within the plv_tresh limits
xidx13=plv13<PLV_targets(2)+plv_thresh & plv13>PLV_targets(2)-plv_thresh; % these are all pssible permutations that will yield a plv12 within the plv_tresh limits
xidx=find(xidx12.*xidx13==1);
% find indices that exist in both xidx12 and xidx13
%
% if isempty(xidx)
%     warndlg('ERROR! PLV_targets could not be found because they are not compatible across signals.\nPLV_trials_est based on random phases. Please revise PLV_targets.\n');
%     rn1=randperm(length(q1)); rn2=randperm(length(q1)); rn3=randperm(length(q1));
%     phase_dist(1,:)=q1(rn1); phase_dist(2,:)=q2(rn2); phase_dist(3,:)=q3(rn3);
%     PLV_trials_est(1)= abs((sum(exp(1i*(squeeze( phase_dist(1,:)-phase_dist(2,:) ))),2))')./num_trials;
%     PLV_trials_est(2)= abs((sum(exp(1i*(squeeze( phase_dist(1,:)-phase_dist(3,:) ))),2))')./num_trials;
%     PLV_trials_est(3)= abs((sum(exp(1i*(squeeze( phase_dist(2,:)-phase_dist(3,:) ))),2))')./num_trials;
%     return
% end

q23idx2=q2idx(xidx,:);    % only using indices for "good" plv12 but keeping same order as that found for plv12
r23=randperm(length(xidx));
q23idx3=q3idx(xidx(r23),:);   % randomizing permutation order across q3idx for those with plv13 meeting plv_thresh criteria

plv23=abs((sum(exp(1i*(squeeze( q2(q23idx2)-q3(q23idx3) ))),2))')./max_perm;
xidx2=xidx(plv23<PLV_targets(3)+plv_thresh & plv23>PLV_targets(3)-plv_thresh); % these are all pssible permutations that will yield a plv12 within the plv_tresh limits
xidx3=xidx(r23(plv23<PLV_targets(3)+plv_thresh & plv23>PLV_targets(3)-plv_thresh)); % these are all pssible permutations that will yield a plv12 within the plv_tresh limits
%         xidx23=xidx(r23(q23idx3(xidx3)));
%         xidx12=xidx(q23idx2(xidx3));

%% if size of remaining num_perms > 100,000 % processing take along time so limiting to <100,000.
num_idx=1e4;
max_iter=max_perm;    % maximum number of iterations before stopping the while loop.
kk=1; kstop=0;
while kk<max_iter
    fprintf('Starting loop #%.f\n',kk);
    if size(xidx2,2)>num_idx
        rn_idx=randperm(size(xidx2,2)); % randomly selecting the indices
        rn_idx=rn_idx(1:num_idx);
        xidx2=xidx2(rn_idx);
        xidx3=xidx3(rn_idx);
    end
    
    clear plv
    % recalculating PLV to only include those values that meet criteria for PLV_targets,
    plv(1,:)=abs((sum(exp(1i*(squeeze( q1(q1idx(xidx2,:))-q2(q2idx(xidx2,:)) ))),2))')./max_perm;   % plv 1 vs 2
    plv(2,:)=abs((sum(exp(1i*(squeeze( q1(q1idx(xidx3,:))-q3(q3idx(xidx3,:)) ))),2))')./max_perm;   % plv 1 vs 3
    plv(3,:)=abs((sum(exp(1i*(squeeze( q2(q2idx(xidx2,:))-q3(q3idx(xidx3,:)) ))),2))')./max_perm;   % plv 1 vs 2
    %% simply filling in next order
    rn1=randperm(num_perms)-1;
    rn2=randperm(num_perms)-1;
    rn3=randperm(num_perms)-1;
    x1=nan(length(xidx2),num_perms,size(q1idx,2)); x2=x1; x3=x1;
    for t=1:num_perms
        x1(:,t,:)=q1idx(xidx2,:)+rn1(t);
        x2(:,t,:)=q2idx(xidx2,:)+rn2(t);
        x3(:,t,:)=q3idx(xidx3,:)+rn3(t);
    end
    x1=reshape(x1,[size(x1,1) size(x1,2)*size(x1,3)]);
    x2=reshape(x2,[size(x1,1) size(x1,2)*size(x1,3)]);
    x3=reshape(x3,[size(x1,1) size(x1,2)*size(x1,3)]);
    clear qdiff;
    qdiff(:,:,1)=[q1(x1)-q2(x2)]'; qdiff(:,:,2)=[q1(x1)-q3(x3)]'; qdiff(:,:,3)=[q2(x2)-q3(x3)]';
    
    % pd_add=repmat([-pi:.01*pi:pi],100,1); % adding phases from -pi to pi to find PLI that is similar to PLV_targets.
    pd_add=[-pi:.05*pi:pi];
    try
        %% Finding PLI values in perms
        fprintf('Finding PLI values. This may take some time...\n')
        pqdiff=single(nan(3,size(x1,2),size(x1,1)));
        dPLI=single(nan(size(pd_add,2),3,size(x1,1)));
        for p=1:size(pd_add,2)
            pqdiff(1,:,:)=[q1(x1)-(q2(x2)+pd_add(p))]';
            pqdiff(2,:,:)=[q1(x1)-(q3(x3)+pd_add(p))]';
            pqdiff(3,:,:)=[(q2(x2)+pd_add(p))-(q3(x3)+pd_add(p))]';
            % finding PLIs that match
            [~,dPLI(p,:,:)]=calc_pli_phases(permute(squeeze(pqdiff(:,:,:)),[1 2 3]));
        end
        %dPLI=dPLI;
        clear xidx;
        xidx(:,1,:)=(dPLI(:,1,:)>PLI_targets(1)-plv_thresh & dPLI(:,1,:)<PLI_targets(1)+plv_thresh);  % finding PLI within limits of PLI_targets
        xidx(:,2,:)=(dPLI(:,2,:)>PLI_targets(2)-plv_thresh & dPLI(:,2,:)<PLI_targets(2)+plv_thresh);
        xidx(:,3,:)=(dPLI(:,3,:)>PLI_targets(3)-plv_thresh & dPLI(:,3,:)<PLI_targets(3)+plv_thresh);
        
        x=squeeze(sum(xidx,2)); % finding perms that have dPLI witihin targets for all three contrasts 1-2,1-3,2-3
        y=zeros(size(x)); y(x==3)=1;  % finding perms that have dPLI witihin targets for all three contrasts 1-2,1-3,2-3
        y=sum(y); % if any phases added to pli amount to getting perm plis within limits to PLI_targets, then y>1; else y=0;
        y(y>0)=1; pli_idx=logical(y);
        pli=nan(size(xidx)); pli_err2=pli;
        for cx=1:3
            pli(:,cx,:)=squeeze(dPLI(:,cx,:)).*pli_idx;
            pli_err2(:,cx,:)=squeeze(pli(:,cx,:))-PLI_targets(cx);
        end
        pli(:,:,pli_idx==0)=nan;
        pli_err2(:,:,pli_idx==0)=nan;   % setting pli not within PLI_targets limits to nan
        pli_err=squeeze(nansum(abs(pli_err2),2)); % perm pli difference to PLI_targets
        pli_err(:,pli_idx==0)=nan;  % setting pli not within PLI_targets limits to nan
        
        
        %% finding best PLV and PLI matches
        plv_err=nansum(abs(plv-PLV_targets'));   % perm plv differenceto PLV_targets
        plv_err(pli_idx==0)=nan;
        
        % multiplying pli_err by plv_err will yield composite score smallest values will reflect best overall fits between pli and plv.
        comp_err=bsxfun(@mtimes,pli_err,plv_err);
        
        [~,perm_idx]=min(min(comp_err));    % best-fit permutation
        [~,phase_idx]=min(comp_err(:,perm_idx));    % best-fit added phase
        
        %     pqdiff(1,:,:)=[q1(x1)-(q2(x2)+pd_add(phase_idx))]';
        %     pqdiff(2,:,:)=[q1(x1)-(q3(x3)+pd_add(phase_idx))]';
        %     pqdiff(3,:,:)=[(q2(x2)+pd_add(phase_idx))-(q3(x3)+pd_add(phase_idx))]';
        %     [~,x]=calc_pli_phases(permute(squeeze(pqdiff(:,:,:)),[1 2 3]));
        
        
        q1_ord=x1(perm_idx,:);
        q2_ord=x2(perm_idx,:);
        q3_ord=x3(perm_idx,:);
        
    catch
        warndlg('PLV_targets or PLI_targets are likely to different across sources to find possible solution');
        return
    end
    % adding in extra trials missed in the permutation if unequal num_trial/max_perms; the phase_diff=0 among these remaining sources 1-2, 1-3, 2-3.
    if length(q1_ord)<num_trials
        x=setdiff(1:num_trials,q1_ord)';
        q1_ord=[q1_ord; x]; q2_ord=[q2_ord; x]; q3_ord=[q3_ord; x];
    end
    
    %     pqdiff(1,:,:)=[q1(x1)-(q2(x2)+pd_add(phase_idx))]';
    %     pqdiff(2,:,:)=[q1(x1)-(q3(x3)+pd_add(phase_idx))]';
    %     pqdiff(3,:,:)=[(q2(x2)+pd_add(phase_idx))-(q3(x3)+pd_add(phase_idx))]';
    %     [~,x]=calc_pli_phases(permute(squeeze(pqdiff(:,:,:)),[1 2 3]));
    
    q2a=q2+pd_add(phase_idx); % adding phase delay to get proper PLI_targets
    q3a=q3+pd_add(phase_idx); % adding phase delay to get proper PLI_targets
    
    PLV_trial_est(1)=abs((sum(exp(1i*(squeeze( q1(q1_ord)-q2a(q2_ord) ))),2))')./length(q1_ord);
    PLV_trial_est(2)=abs((sum(exp(1i*(squeeze( q1(q1_ord)-q3a(q3_ord) ))),2))')./length(q1_ord);
    PLV_trial_est(3)=abs((sum(exp(1i*(squeeze( q2a(q2_ord)-q3a(q3_ord) ))),2))')./length(q1_ord);
    clear pqdiff;
    pqdiff(1,:)=q1(q1_ord)-q2a(q2_ord);
    pqdiff(2,:)=q1(q1_ord)-q3a(q3_ord);
    pqdiff(3,:)=q2a(q2_ord)-q3a(q3_ord);
    
    [PLI_trials_est,dPLI_trials_est]=calc_pli_phases(squeeze(pqdiff));
    PLI_trials_est=PLI_trials_est'; dPLI_trials_est=dPLI_trials_est';
    %% Stopping criteria
    plv_diff_est=abs(PLV_trial_est-PLV_targets);
    if sum(plv_diff_est>plv_thresh)>0
        %         warndlg(sprintf('PLV Estimates = %.3f %.3f %.3f\n   are not within limits of\nPLV Targets = %.3f %.3f %.3f', PLV_trial_est,PLV_targets));
        fprintf('PLV Estimates = %.3f %.3f %.3f are not within limits of PLV Targets = %.3f %.3f %.3f\n', PLV_trial_est,PLV_targets);
        kk=kk+1;
        kstop=0;
    else
        kstop=1;
    end
    pli_diff_est=abs(dPLI_trials_est-PLI_targets');
    if sum(pli_diff_est>=plv_thresh)>0
        %         warndlg(sprintf('PLI Estimates = %.3f %.3f %.3f\n   are not within limits of\nPLI Targets = %.3f %.3f %.3f', dPLI_trials_est,PLI_targets));
        fprintf('PLI Estimates = %.3f %.3f %.3f are not within limits of PLI Targets = %.3f %.3f %.3f\n', dPLI_trials_est,PLI_targets);
        kk=kk+1;
        kstop=0;
    else
        kstop=kstop+1;
    end
    if kstop==2     % stopping loop once PLV and PLI are within targets
        kk=max_iter+1;
    end
end

pli_diff_est=abs(dPLI_trials_est-PLI_targets');
if sum(pli_diff_est>=plv_thresh)>0
    warndlg(sprintf('\nPLI Estimates = %.3f %.3f %.3f are not within +/- %.2f of \nPLI Targets = %.3f %.3f %.3f\nPlease Try Again.', dPLI_trials_est,plv_thresh,PLI_targets));
end
plv_diff_est=abs(PLV_trial_est-PLV_targets);
if sum(plv_diff_est>plv_thresh)>0
    warndlg(sprintf('PLV Estimates = %.3f %.3f %.3f are not within +/- %.2f of \nPLV Targets = %.3f %.3f %.3f\nPlease Try Again.', PLV_trial_est,plv_thresh,PLV_targets));
end

q1sort=q1(q1_ord); q2sort=q2a(q2_ord); q3sort=q3a(q3_ord);

%%
% pd=permute(repmat(cat(1,q1sort,q2sort,q3sort),1,1,256),[3 1 2]);
% [PLI]=calc_PLI_ath(pd,256,[0:1/256:1-1/256],.05,0,nchoosek(1:3,2),0,0);
% [PLV]=calc_PLV_ath(pd,nchoosek(1:3,2),0,0);
% fprintf('PLV = %.2f %.2f %.2f\n',PLV.PLV(:,1));
% fprintf('PLI = %.2f %.2f %.2f\n',PLI.PLI(:,1));
% fprintf('dPLI = %.2f %.2f %.2f\n',(PLI.dPLI(:,1)-.5)*2);

%% Finalizing output
phase_dist(1,:)=q1sort;
phase_dist(2,:)=q2sort;
phase_dist(3,:)=q3sort;
PLV_trials_est(1)= abs((sum(exp(1i*(squeeze( q1sort-q2sort ))),2))')./num_trials;
PLV_trials_est(2)= abs((sum(exp(1i*(squeeze( q1sort-q3sort ))),2))')./num_trials;
PLV_trials_est(3)= abs((sum(exp(1i*(squeeze( q2sort-q3sort ))),2))')./num_trials;

clear rn1 rn2 x;
for t=1:10; rn1(t,:)=randperm(length(q1sort)); rn2(t,:)=randperm(length(q1sort));end
x(1,:)= abs((sum(exp(1i*(squeeze( q1sort(rn1)-q1sort(rn2) ))),2))')./num_trials;
x(2,:)= abs((sum(exp(1i*(squeeze( q2sort(rn1)-q2sort(rn2) ))),2))')./num_trials;
x(3,:)= abs((sum(exp(1i*(squeeze( q3sort(rn1)-q3sort(rn2) ))),2))')./num_trials;
PLV_surrogate_est=nanmean(x,2);


fprintf('Final PLV Surrogate [Target (1-2)=%.3f (1-3)=%.3f (2-3)=%.3f Found (1-2)=%.3f (1-3)=%.3f (2-3)=%.3f]\n',PLV_evoked,PLV_surrogate_est);
fprintf('Final PLV Trial [Target (1-2)=%.3f (1-3)=%.3f (2-3)=%.3f Found (1-2)=%.3f (1-3)=%.3f (2-3)=%.3f]\n',PLV_targets,PLV_trials_est);
fprintf('Final PLI Trial [Target (1-2)=%.3f (1-3)=%.3f (2-3)=%.3f Found (1-2)=%.3f (1-3)=%.3f (2-3)=%.3f]\n',PLI_targets,dPLI_trials_est);

% randomizing phases across trials but keep signal order.
rn=randperm(size(phase_dist,2));
phase_dist=phase_dist(:,rn);

function [PLI,dPLI]=calc_pli_phases(pd_diff)
% PLI
my_sine = round(sin(pd_diff)*(10^6))/(10^6); % For MATLAB, 0 can be just an approximate.  Make approximates real zeros.
sign_test = sign(my_sine);PLI = squeeze(abs(mean(squeeze(nanmean(sign_test, 2)), 3)));
% dPLI
Y = zeros(size(my_sine)); Y(my_sine > 0) = 1; Y(my_sine == 0) = .5;
dPLI = 2*(squeeze(mean(mean(Y, 2), 4))-.5);



