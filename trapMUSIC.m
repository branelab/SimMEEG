function sOut = trapMUSIC(sIn)
%
% SYNTAX:
%
% Given full and noise covariances and a set of lead fields, finds
% specified number of sources using TRAP-MUSIC algorithm
%
% Input:
%   sIn     - input parameters structure with fields:
%       R       - (M x M) full covariance of the data
%       arrN    - (M x M) noise covariance of the data
%       *** either group of fields specifies forward solutions **********
%       *** if both are present then arrH/lstFlag/dims group is used ****
%           arrH    - (nV x 3 x M) lead fields array
%           lstFlag - (nV x 1) if flag = 1 - this voxel and lead field is used
%           dims    - (3 x 1) [nX,nY,nZ]' = number of grid nodes along each
%                      direction
%               OR
%           lfFile  - full path to lead fields data file (typically
%                     leadFields.dat) with format described in readLeadFields()
%           lstVox  - (nV x 3) list of voxel locations in meters; if empty - no
%                     check if voxels really correspond to FS in lfFile is done
%                     !!! IMPORTANT !!! When lstVox was generated from rect
%                     grid, it is assumed that C++ order is used - that is the
%                     Z coordinate of the voxels is the one changing fastest
%       *****************************************************************
%       nSrc    - number of sources to find
%
%       --- Optional fields ---------------------
%       pVal    - (default 1) the p-value threshold for the peaks. Only
%                 peaks higher than (1 - pVal) quantile boundary will be
%                 taken into account. Set to 1 to get all peaks. NOTE:
%                 can't use "alpha" as a field name because it clashes with
%                 built in func in latest MATLABs
%       gap     - (default 2) minimal distance in voxels between maxima. Min allowed value is 2
%       bPlotLambda     - (default FALSE) if !0 - the labmda values (see
%                         below) will be plotted
%       bVerbose    - (default FALSE) if !0 - print out iterations results
%
%   sOut  - output structure with fields:
%       arrH    - (nV x 3 x M) lead fields array
%       lstFlag - (nV x 1) if flag = 1 - this voxel and lead field is used
%       rtNm1   - (M x M) arrN^(-1/2)
%       lambda  - (M x 1) eigenvalues of the problem R e = lambda N e
%       EV      - (M x M) eigenvectors for lambda in whitened basis
%       lstMu   - (nSrc x 1) localizer values for found sources
%       idxSrc  - (nSrc x 1) source voxels
%       U3D     - (nSrc x 3) source orientations
%       Wr      - (M x nSrc) The traditional MCMV beamformer weigths for
%                 found sources based on the full covariance
%       Wr      - (M x nSrc) MCMV beamformer weigths based on the noise covariance
%       fImg    - (nSrc x nV) localizer spatial distribution, with NaNs set
%                 outside the head bounds
%
% A.Moiseev, BCNI, Aug 2019

lstMandatoryArgs = {'R','arrN','nSrc'};
bReadFSFromFile = false;    % Flag to read lead fields from file

% Check the inputs
if ~all(isfield(sIn, lstMandatoryArgs))
    error('Some of the mandatory args are missing');
end

if ~isfield(sIn, 'bVerbose')
    sIn.bVerbose = false;
end % if ~isfield(sIn, 'bVerbose')

if ~isfield(sIn, 'arrH')
    if ~isfield(sIn, 'lfFile')
        error('Either arrH or lfFile should be specified');
    else
        if ~isfield(sIn,'lstVox')
            error('lstVox field is missing');
        end % if ~isfield(sIn,'lstVox')
        
        bReadFSFromFile = true;
    end % if ~isfield(sIn, 'lfFile')
else    % arrH is specified
    if ~isfield(sIn, 'lstFlag')
        error('lstFlag field is missing');
    end % if ~isfield(sIn, 'lstFlag')
    
    if ~isfield(sIn, 'dims')
        error('dims field is missing');
    end % if ~isfield(sIn, 'lstFlag')
end % if ~isfield(sIn, 'arrH')

if ~isfield(sIn,'pVal')
    sIn.pVal = 1;
end % if ~isfield(sIn,'pVal')

if ~isfield(sIn,'gap')
    sIn.gap = 2;
end % if ~isfield(sIn,'gap')

if ~isfield(sIn,'bPlotLambda')
    sIn.bPlotLambda = false;
end % if ~isfield(sIn,'bPlotLambda')

% Generate variables corresponding to the field names of the input
% structure
args = fieldnames(sIn);

for iArg = 1:length(args)
    eval([args{iArg},' = sIn.',args{iArg},';']);
end % for iArg = 1:length(args)

if bReadFSFromFile
    [arrH, lstFlag] = readLeadFields(lfFile, lstVox); %#ok<NODEF>
    % Needed to calc dims for svlPeak
    ROI = [min(lstVox(:,1)), max(lstVox(:,1)); min(lstVox(:,2)), max(lstVox(:,2)); ...
                    min(lstVox(:,3)), max(lstVox(:,3))];
    step = abs(lstVox(1,1,2) - lstVox(1,1,1));   % Grid step
    dims = round((ROI(:,2)-ROI(:,1))/step+1);
end % if bReadFSFromFile

M = size(R,1);      %#ok<NODEF> % Number of channels
sOut.arrH = arrH;
sOut.lstFlag = lstFlag;

% --------------------------------------------------------------
% Transform to whitened coordinates. Reuse the same variables
% --------------------------------------------------------------
rtNm1 = sqrtm(invSPD(arrN));
R = rtNm1 * R * rtNm1;

% Prewhiten arrH
% There seem to be no tensor product in matlab. So do it in a loop
nV = size(arrH,1);
for iV = 1:nV
    tmp = squeeze(arrH(iV,:,:))*rtNm1;
    arrH(iV,:,:) = tmp;
end % for iV = 1:nV

sOut.rtNm1 = rtNm1;

% --------------------------------------------------------------
% Get the source subspace from the EV problem R u = lmbda u
% --------------------------------------------------------------
[EV, D] = eig((R+R')/2);    % Enforce EXACT symmetry of R to always have real-valued EV

% Sort EVs properly
[lambda,idx] = sort(diag(D),'descend');
EV = EV(:,idx);

if bPlotLambda
    semilogy((1:M)', abs(lambda));  % abs() to avoid warning because of tiny negative lambda
end % if bPlotLambda

sOut.lambda = lambda;
sOut.EV = EV;

% Take the 1st nSrc dimenvectors:
% lambda = lambda(1:nSrc); - unused
EV = EV(:,1:nSrc);

% --------------------------------------------------------------
% Do the iterations
% --------------------------------------------------------------
lstMu = NaN(nSrc,1);    % Localizer values
idxSrc = NaN(nSrc,1);   % Idx of the source voxels
U3D = NaN(3,nSrc);      % Source orientations
Hsrc = NaN(M,nSrc);     % Forward solutions of the found sources
U = EV;                 % Current source subspace
P = U*U';               % Current projector on the source subspace
Q = eye(M);             % Current out-projector of the already found sources
fImg = NaN(nSrc,nV);    % Functional images for each iteration step
bMatlab = false;        % This is because lstVox are listed in C++ order

if bVerbose
    fprintf('TRAP-MUSIC: starting iterations...\n');
end % if bVerbose

for iS = 1:nSrc
%       arrH    - (nV x 3 x M) lead fields array
%       lstFlag - (nV x 1) if flag = 1 - this voxel and lead field is used
    f = fImg(iS,:)';
    u = NaN(3,nV);      % Tmp storage for orientation of each voxel

    for iV = 1:nV
        if ~lstFlag(iV); continue; end
        
        Ht = squeeze(arrH(iV,:,:)); % Transposed vector lead field of a voxel
        Ht = Ht*Q;                  % Project out found sources
        H = Ht';
        HtH = Ht * H;
        
        if rank(HtH) < 3 % Exclude voxels with degenerate H, specifically found source voxels
            f(iV) = 0;
            u(:,iV) = [0,0,0]';  % Make sure we get an error if trying to use this u
            continue;
        end % if rank(HtH) < 3
        
        % The localizer value is max eigenvalue of the EV problem
        HtPH = Ht * P * H; HtPH = (HtPH + HtPH')/2;
        HtH = (HtH + HtH')/2;    % Enforce EXACT symmetry of eig() args
        [Usrc, D] = eig(HtPH, HtH);               
        [D,idx] = sort(diag(D),'descend');
        
        if any(abs(D)) > 1  % Just in case - this should never happen
            error('Found an eigenvalue > 1 - something is wrong');
        end % if any(D) > 1
        
        f(iV) = D(1);               % Localizer value
        u(:,iV) = Usrc(:,idx(1));   % Unnormalized orientation
    end % for iV = 1:nV
    
    % Search for peaks
    % We usd sIn.gap instead of just "gap" due to problems running in environments with all toolboxes connected
    % because sometimes variable gap clashes with function gap() defined in some of them 
    [peaks, arrIdx3D, idxMax] = svlPeak(f, dims, pVal, sIn.gap, bMatlab); %#ok<ASGLU>
    
    if isempty(peaks)
        % A pathological case, but it does happen (rarely)
        % Just break the search cycle after cleaning NaNs
        fImg(iS:end,:) = [];
        lstMu(iS:end) = [];
        idxSrc(iS:end) = [];
        U3D(:,iS:end) = [];
        break;
    end % if isempty(peaks)
    
    % Source found. Save the results
    fImg(iS,:) = f;
    lstMu(iS) = peaks(1);
    idxSrc(iS) = idxMax(1);
    uMax = u(:,idxMax(1));
    U3D(:,iS) = uMax/norm(uMax);
    Hsrc(:,iS) = squeeze(arrH(idxMax(1),:,:))' * uMax;    % Source scalar whitened LF
    
    % Project out found sources from the source subspace and find reduced
    % subspace
    [Usrc,~] = svd(Hsrc(:,1:iS),0);     % ONB of the found sources
    Q = (eye(M) - Usrc*Usrc');          % New out-projector
    
    %==========================================================
    % Working version:
    U = Q*U;                            % New source subspace, not truncated
    % This variant outprojects from the source subspace EV, not U
    % (verified to yield the same result up to rounding errors)
    % U = Q*EV;                            % New source subspace, not truncated
    %==========================================================
    
    [U,D] = svd(U,0);                   % ONB and singular values of the new source subspace
    
    [~,idx] = sort(diag(D),'descend');
    U = U(:,idx);
    
    %==========================================================
    % Working version:
    % TRAP step: remove dimension with the smallest singular value
    U = U(:,1:(end-1));    
    
    % Variant with outprojecting from original source subspace EV:
    % - truncate iS smallest SVs (produces practically same res as current
    % working version
    % U = U(:,1:(end-iS));    
    %==========================================================
    
    % Update current projector on the source subspace
    P = U*U';
    
    if bVerbose
        fprintf('Source %d of %d found: mu = %g\n', iS, nSrc, lstMu(iS));
    end % if bVerbose
end % for iS = 1:nSrc

sOut.lstMu = lstMu;
sOut.idxSrc = idxSrc;
sOut.U3D = U3D';

% Calculate R- and N-based beamformer weights in the original (not pre-whitened) basis
% The expression is W = N^(-1/2)* Wwhitened
Rm1Hsrc = invSPD(R)*Hsrc;
sOut.Wr = rtNm1 * Rm1Hsrc * invSPD(Hsrc'* Rm1Hsrc);
sOut.Wn = rtNm1 * Hsrc * invSPD(Hsrc'* Hsrc);   % For N-weights, simply set Rm1 = I

sOut.fImg = fImg;

