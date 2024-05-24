function sOut = doSMCMV(sIn)
%
% SYNTAX:
%   sOut = doSMCMV(sIn)
%
% Implements source subspace based MCMV source localization with or without
% the TRAP trunctaion of the source supspace.
%
% Procedure:
% 1. Transfer to whitened FS and covariances
%
% 2. Given the number nSrc of sources to look for, construct "source covariance"
%    Rs and source subspace Usrc (M x nSrc) as sensor-level covariance
%    based on 1st nSrc principal components, and subspace spanned by those, respectively 
%
% 3. For each of nSrc iterations:
%
% 4.    For each voxel, do:
%       - Find this voxel's best orientation u by projecting on the full or
%         truncated source subspace by finding max e-vector of the problem
%           Ht Usrc Usrc' H u = mu Ht H u
%         (whith the max eigenvalue mu being actually equal to MUSIC localizer)
%       - Add voxel's scalar FS h = H(iv)*u to the set of found FS
%       - Create an M-dimensional out-projector Q(k,iv) of the set [h1,...,hk,h]
%       - Directly calculate corresponding MCMV localizer P as a trace
%         expresssion:
%           P = tr(Usrc'*Q*A*Q*Usrc*(Usrc'*Q*B*Q*Usrc)^(-1)) - (k+1)
%         with A,B being matrices corresponding to the localizer
%       - Use 1/((mu-1)(P+1)) as a localizer value because svlPeak looks
%         for maxima, not minima
%
% 5. After the voxels cycle:
%       - Largest max found by the svlPeak() is the new source
%         position/orientation
%       - Append scalar FS of this source to already found sources
%       - Calculate ONB of subspace spanned by the found sources: Usrc
%       - Calculate new out-projector operator for all found sources:
%           Q = I - Ufound Ufound'
%       - Apply this out-projector to the current source subspace U and find
%         ONB for this new source subspace:
%           tmp1 = Q*Usrc
%           tmp2 = ONB(tmp1)
%       - TRAP step: trancate 1 dim to get the new source subspace
%           Usrc = drop 1 smallest dim from tmp2 
%       - Calculate projector on the new source subspace:
%           Psrc = Usrc * Usrc'
%     This ends the iteration cycle over the sources (RAP-cycle)
%
% For the moment, we do not process references in the sense MCMV does, that
% is - for nulling known sources. This is because all sources (and more) are already
% included in the source subspace, and providing them explicitly in ideal
% world would change nothing, in real world - inflate the source subspace
% if provided references will not exadctly belong to Usrc constructed using
% covariance matrix. In the future "references" can be added as additional
% null constraints, by simply merging their space with source space.
%
% Input:
%   sIn     - input structure with the following fields:
%       --------- Mandatory ---------------------------------------------
%       beamType  (string) The beamformer type. One of: 'MPZ', 'MAI'
%                 ("power" beamformers), 'MER', 'rMER' ("evoked"
%                 beamformers)
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
%       nSrc    - max number of sources to find. This parameter determines max
%                 beamformer order used for SMCMV and is also the
%                 preset dimension of the source space.
%       bRAPBeam - type of MCMV localizer to combine with mu-localizer 
%               if TRUE (default), RAP beamformer-based localizer will be
%               calculated, which is
%                   P = RAP_Beam/(1-mu)
%               RAP_Beam being MXX "single source" RAP beamformer;
%               if FALSE, SMCMV beamformer-based localizer will be returned,
%               which is
%                   P = (((1-mu)(L_MXX + 1))^(-1)
%           
%       --- Optional fields ---------------------
%       Cavg    - (M x M) matrix of 2nd moments of the averaged field;
%                   mandatory if evoked beamformer (MER, rMER) is used
%       pVal    - (default 1) the p-value threshold for the peaks. Only
%                 peaks higher than (1 - pVal) quantile boundary will be
%                 taken into account. Set to 1 to get all peaks. NOTE:
%                 can't use "alpha" as a field name because it clashes with
%                 built in func in latest MATLABs
%       gap     - (default 2) minimal distance in voxels between maxima.
%                 Min allowed value is 2 
%       bDoTRAP (default FALSE)
%               - if set to FALSE, then the TRAP step will not be applied
%       bPlotLambda
%               - (default FALSE) if !0 - the labmda values of EV 
%                  R u = lambda N u  will be plotted
%       bVerbose - (default FALSE) if !0 - print out iterations results
%                   Optional ground truth data:
%       trueSrcVox  - (default - missing) (nSrc x 1) true src voxel numbers
%                     If this field is present, then the next one (U3D) is
%                     assumed to be present also. In this case, for each
%                     iteration localizer values for all the true source
%                     will be printed
%       trueSrcU3D  - (default - missing) (nSrc x 3) true src orientations
%       trueSrcPos  - (default - missing) (nSrc x 3) true src positions, m
%                     (not used so far - left here for future needs)
%
% Output:
%   sOut    - output structure with the following fields:
%       lstVox  - (nV x 3) either lstVox supplied on input or that read
%                 from the lfFile, or nothing if neither was done
%       samHeader - as is, if lead fields were read from file
%       arrH    - (nV x 3 x M) lead fields array
%       lstFlag - (nV x 1) if flag = 1 - this voxel and lead field is used
%       rtNm1   - (M x M) arrN^(-1/2)
%       lambda  - (M x 1) eigenvalues of the problem R e = lambda N e
%       EV      - (M x M) eigenvectors for lambda in whitened basis
%       lstMu   - (nSrc x 1) localizer values for found sources
%       idxSrc  - (nSrc x 1) source voxels
%       U3D     - (nSrc x 3) source orientations
%       Wr      - (M x nSrc) The traditional MCMV beamformer weigths based
%                 on the full covariance
%       Wr      - (M x nSrc) MCMV beamformer weigths based on the noise covariance
%       fImg    - (nSrc x nV) localizer spatial distribution, with NaNs set
%                 outside the head bounds
%       bRAPBeam
%               - the actual value of bRAPBeam flag used in calculations
%       bDoTRAP
%               - the actual value of bDoTRAP flag used in calculations
%       ----------- Optional outputs -------------------------------------
%       samHeader  SAM file header v.2; only returned if lfFile instead of
%                  arrH was specified on input
%
% A.Moiseev, BCNI, Dec 2019.

lstMandatoryArgs = {'beamType', 'R','arrN','nSrc'};

lstBeam = {'MPZ','MAI','MER','RMER'};

% Beamformer type enumerator definition. It will be passed to all internal
% functions explicitly, to avoid using function nesting. The latter is
% prohibited because it does not allow create/assign variables dynamically
BEAMTYPE.IMPZ = 1;
BEAMTYPE.IMAI = 2;
BEAMTYPE.IMER = 3;
BEAMTYPE.IRMER = 4;

bReadFSFromFile = false;    % Flag to read lead fields from file

% Check the inputs
% Mandatory:
if ~all(isfield(sIn, lstMandatoryArgs))
    error('Some of the mandatory args are missing');
end

tmp = 1:length(lstBeam);
iBeam = ismember(lstBeam,upper(sIn.beamType));
iBeam = tmp(iBeam); % iBeam is the numerical beamforme type

if isempty(iBeam)
    error('Invalid beamformer type specified: %s', sIn.beamType);
end % if isempty(iBeam)

bEvoked = false;

if (iBeam == BEAMTYPE.IMER) || (iBeam == BEAMTYPE.IRMER)
    bEvoked = true;
    
    if ~isfield(sIn, 'Cavg')
        error('Cavg must be supplied if evoked beamformer (MER or rMER) is used');
    end % if ~isfield(sIn, 'Cavg')
end % if (iBeam == IMER) || (iBeam == IRMER)

if ~isfield(sIn, 'bRAPBeam')
    sIn.bRAPBeam = true;
end % if ~isfield(sIn, 'bRAPBeam')

if ~isfield(sIn, 'bDoTRAP')
    sIn.bDoTRAP = false;
end % if ~isfield(sIn, 'bDoTRAP')

if ~isfield(sIn, 'bVerbose')
    sIn.bVerbose = false;
end % if ~isfield(sIn, 'bVerbose')

if ~isfield(sIn, 'arrH')    % if arrH is not supplied
    if ~isfield(sIn, 'lfFile')
        error('Either arrH or lfFile should be specified');
    else    % lead fields file is specified
        if ~isfield(sIn,'lstVox')
            lstVox = [];
        end % if ~isfield(sIn,'lstVox')
        
        bReadFSFromFile = true;
    end % if ~isfield(sIn, 'lfFile')
else    % arrH is specified
    if ~isfield(sIn, 'lstFlag')
        error('lstFlag field is missing');
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

% Check if we need to calculate localizer for the true sources
if isfield(sIn, 'trueSrcVox')
    bCalcTrue = true;
else
    bCalcTrue = false;
end % isfield(sIn, 'trueSrcVox')

% Generate variables corresponding to the field names of the input
% structure
args = fieldnames(sIn);

for iArg = 1:length(args)
    eval([args{iArg},' = sIn.',args{iArg},';']);
end % for iArg = 1:length(args)

if bReadFSFromFile
    [arrH, lstFlag, ~, samHeader, lstVox] = readLeadFields(lfFile, lstVox);
    % Needed to calc dims for svlPeak
    ROI = [min(lstVox(:,1)), max(lstVox(:,1)); min(lstVox(:,2)), max(lstVox(:,2)); ...
                    min(lstVox(:,3)), max(lstVox(:,3))];
    step = abs(lstVox(1,1,2) - lstVox(1,1,1));   % Grid step
    dims = round((ROI(:,2)-ROI(:,1))/step+1);
    
    sOut.lstVox = lstVox;
    sOut.samHeader = samHeader;
end % if bReadFSFromFile

sOut.arrH = arrH;
sOut.lstFlag = lstFlag;

% Number of channels
M = size(R,1); %#ok<NODEF>

% --------------------------------------------------------------
% Transform to whitened coordinates. Reuse the same variables
% --------------------------------------------------------------
if bEvoked
    Rm1 = invSPD(R);    % Inverse in original basis
    rtN = sqrtm(arrN);
    arrE = Rm1 * Cavg * Rm1;
    arrE = rtN * arrE * rtN;    % arrE in pre-whitened basis
    arrE = (arrE + arrE')/2;    % Enforce exact symmetry
else
    arrE = [];
end % if bEvoked

rtNm1 = sqrtm(invSPD(arrN));
R = rtNm1 * R * rtNm1;
Rm1 = invSPD(R);
Rm2 = Rm1 * Rm1;

% Ensure all args to eig() are exactly symmetrical
R = (R + R')/2;
Rm1 = (Rm1 + Rm1')/2;
Rm2 = (Rm2 + Rm2')/2;

% Pre-whiten arrH
% There seem to be no tensor product in matlab. So do it in a loop
nV = size(arrH,1);
for iV = 1:nV
    tmp = squeeze(arrH(iV,:,:))*rtNm1;
    arrH(iV,:,:) = tmp;
end % for iV = 1:nV

sOut.rtNm1 = rtNm1;

% -----------------------------------
% Get the source subspace and lambda
% -----------------------------------
[EV, D] = eig(R);

% Sort EVs properly
[lambda,idx] = sort(diag(D),'descend');
EV = EV(:,idx);

if bPlotLambda
    semilogy((1:M)', abs(lambda));  % abs() to avoid warning because of tiny negative lambda
end % if bPlotLambda

sOut.lambda = lambda;
sOut.EV = EV;

% Take the 1st nSrc eigenvectors:
% lambda = lambda(1:nSrc); - no need

% This is source subspace
Usrc0 = EV(:,1:nSrc);  % TODO: VERIFY EV IS ORTHOGONAL MATRIX!! - YES every time I checked
Usrc = Usrc0;

% Calculate various matrices for gen EV problems depending on the beamformer
% type
Usrc2 = Usrc * Usrc';   % Current projector on the source subspace
Q = eye(M);             % Current out-projector of the already found sources
lstMu = NaN(nSrc,1);    % Localizer values
idxSrc = NaN(nSrc,1);   % Idx of the source voxels
U3D = NaN(3,nSrc);      % Source orientations
Hsrc = NaN(M,nSrc);     % Forward solutions of the found sources
fImg = NaN(nSrc,nV);    % Functional images for each iteration step
bMatlab = false;        % This is because lstVox are listed in C++ order

if bVerbose
    if bRAPBeam
        sSubType = 'RAP-beam';
    else
        sSubType = 'Subspace-beam';
    end % if bRAPBeam
    
    fprintf('doSMCMV: starting %s (%s) search...\n', beamType, sSubType);
end % if bVerbose

% --------------------------------------------------------------
% Do the RAP iterations
% --------------------------------------------------------------
for iS = 1:nSrc
    f = fImg(iS,:)';
    uImg = NaN(3,nV);      % Tmp storage for orientation of each voxel
    
    for iV = 1:nV
        if ~lstFlag(iV); continue; end        

        Ht = squeeze(arrH(iV,:,:)); % Transposed vector lead field of a voxel
        Ht = Ht*Q;                  % Project out found sources
        H = Ht';    % M x 3 lead field for the voxel

        [f(iV), uImg(:,iV)] = calcSMSMV4Voxel(iBeam, Usrc, Usrc2, H, Ht, ...
            Hsrc(:,1:(iS - 1)), Rm1, Rm2, arrE, BEAMTYPE, bRAPBeam);
    end % for iV = 1:nV

    % We usd sIn.gap instead of just "gap" due to problems running in environments with all toolboxes connected
    % because sometimes variable gap clashes with function gap() defined in some of them 
    [peaks, ~, idxMax, lstBorderMax] = svlPeak(f, dims, pVal, sIn.gap, bMatlab);

    if isempty(peaks)
        % A pathological case, but it does happen (rarely)
        % Just break the search cycle after cleaning NaNs
        fImg(iS:end,:) = [];
        lstMu(iS:end) = [];
        idxSrc(iS:end) = [];
        U3D(:,iS:end) = [];
        break;
    end % if isempty(peaks)
    
    % Sources found. Save the results
    fImg(iS,:) = f;
    lstMu(iS) = peaks(1);
    idxSrc(iS) = idxMax(1);
    U3D(:,iS) = uImg(:,idxSrc(iS));
    Hsrc(:,iS) = squeeze(arrH(idxSrc(iS),:,:))' * U3D(:,iS);    % Source scalar whitened LF    
    
    % Optionally, display the results for current and true sources
    if bCalcTrue
        fprintf('Iteration %d: the localizer max found = %g. Max for true sources:\n', iS, lstMu(iS));
        
        for iTrue = 1:length(trueSrcVox)
            Ht = squeeze(arrH(trueSrcVox(iTrue),:,:)); % Transposed vector lead field of a voxel
            Ht = Ht*Q;                  % Project out found sources
            H = Ht';    % M x 3 lead field for the voxel
            
            % Here we calculate on true sources USING APPROXIMATELY FOUND
            % sources of the previous iterations. At least at first
            % iterations, one of the true max should be larger than found
            % max
            pTrue = calcSMSMV4Voxel(iBeam, Usrc, Usrc2, H, Ht, ...
                Hsrc(:,1:(iS - 1)), Rm1, Rm2, arrE, BEAMTYPE);
            fprintf('True %d: %g\n', iTrue, pTrue);            
        end % for iTrue = 1:length(trueSrcVox)
        
        % List discarded border maxima, if any
        for iB = 1:size(lstBorderMax,1)
            iV = idx3Dto1D(lstBorderMax(iB,:), dims, bMatlab);
            
            Ht = squeeze(arrH(iV,:,:)); % Transposed vector lead field of a voxel
            Ht = Ht*Q;                  % Project out found sources
            H = Ht';    % M x 3 lead field for the voxel
            
            % Calculate discarded border max
            pBdr = calcSMSMV4Voxel(iBeam, Usrc, Usrc2, H, Ht, ...
                Hsrc(:,1:(iS - 1)), Rm1, Rm2, arrE, BEAMTYPE);
            fprintf('Discarded border max %d: %g\n', iB, pBdr);            
        end % for iB = 1:size(lstBorderMax,1)
        
        fprintf('\n');
    end % if bCalcTrue
    
    % Do the out-projections and truncations for the Usrc
    [Ufnd,~] = svd(Hsrc(:,1:iS),0);     % ONB of the found sources
    Q = (eye(M) - Ufnd*Ufnd');          % New out-projector
    
%     Usrc = Q*Usrc;                      % New source subspace, not truncated, recursively calculated
    Usrc = Q*Usrc0;                     % New source subspace, not truncated - directly calculated
                                        % (same result - checked) 
    [Usrc,D] = svd(Usrc,0);             % ONB and singular values of the new source subspace
    % NOTE: we can't avoid the above SVD irrespective to bDoTRAP, because
    % it is used in the Usrc2 (the source space projector) calculation
    
    if bDoTRAP  % TRAP step for the source subspace: remove dimension with the smallest singular value
        [D,idx] = sort(diag(D),'descend'); %#ok<ASGLU>
        Usrc = Usrc(:,idx);

%         Usrc = Usrc(:,1:(end-1));    % Discard just 1 SV for recursive calculation of Usrc
        Usrc = Usrc(:,1:(end-iS));      % Discard iS singular values for direct calc of Usrc
                                        % (same result - checked) 
    end % if bDoTRAP
    
    % Update current projector on the source subspace
    Usrc2 = Usrc*Usrc';
        
    if bVerbose
        fprintf('Source %d of %d found: mu = %g\n', iS, nSrc, lstMu(iS));
    end % if bVerbose
end % for iS = 1:nSrc

sOut.lstMu = lstMu;
sOut.idxSrc = idxSrc;
sOut.U3D = U3D';

% Calculate R- and N-based beamformer weights in the original (not pre-whitened) basis
% The expression is W = N^(-1/2)* Wwhitened
Rm1Hsrc = Rm1*Hsrc;
sOut.Wr = rtNm1 * Rm1Hsrc * invSPD(Hsrc'* Rm1Hsrc);
sOut.Wn = rtNm1 * Hsrc * invSPD(Hsrc'* Hsrc);   % For N-weights, simply set Rm1 = I

sOut.fImg = fImg;
sOut.bRAPBeam = bRAPBeam;
sOut.bDoTRAP = bDoTRAP;

%---------------------------------------------
function [P, u] = calcSMSMV4Voxel(iBeam, Usrc, Usrc2, H, Ht, Hsrc, Rm1, Rm2, arrE, BEAMTYPE, bRAPBeam)
%---------------------------------------------
% Calculate the SMCMV localizer value and correpsonding optimal orientation
%
% Input:
%   iBeam   - (int) beamformer type, one of IMPZ = 1, IMAI, IMER, IRMER
%   Usrc    - (M x nSrc) Full (no-TRAP), or truncated (TRAP) source subspace
%             which will be subjecto to out-projection of the (found + test) sources
%   Usrc2   - (M x M) Src subspaces product: Usrc*Usrc'
%   H       - (M x 3) Vector FS for a voxel (after out-projecting already
%             found sources)
%   Ht      - (3 x M) H'
%   Hsrc    - (M x nFound) Matrix of already found source lead fields (in
%             the original sensor space - no out-projections, of course)
%   Rm1     - (M x M) inverse of covariance (pre-whitened)
%   Rm2     - (M x M) inverse of covariance (pre-whitened) squared: Rm1*Rm1
%   arrE    - (M x M) (whitened) Rm1*Cavg*Rm1 for MER beamformer, or empty 
%   BEAMTYPE - beamformer type enumerator
%   bRAPBeam -
%             if TRUE, RAP beamformer-based localizer will be calculated,
%             which is
%               P = RAP_Beam/(1-mu)
%             RAP_Beam being MXX "single source" RAP beamformer;
%             if FALSE, SMCMV beamformer-based localizer will be returned,
%             which is
%               P = (((1-mu)(L_MXX + 1))^(-1)
%
% Output:
%   P       - (double) the localizer value (in this case 1/trace(...))
%   u       - (3 x 1) optimal orientation

% Find optimal orientation as the one with largest projection on source subspace
HtUUH = Ht*Usrc2*H; HtUUH = (HtUUH + HtUUH')/2;
HtH = Ht*H; HtH = (HtH + HtH')/2;   % Ensure exact symmetry, as HtUUH, HtH are used by eig() 

if rank(HtH) < 3 % Exclude voxels with degenerate H, specifically found source voxels
    P = 0;
    u = [0,0,0]';  % Make sure we get an error if trying to use this u
    return;
end % rank(HtH) < 3
                
% The orientation is max EV of a problem Ht Usrc Usrc' H u = mu Ht H u
[U,D] = eig(HtUUH,HtH);
[mu,idx] = sort(diag(D),'descend');

if any(abs(mu)) > 1  % Just in case - this should never happen
    error('Found an eigenvalue > 1 - something is wrong');
end % if any(abs(mu)) > 1

U = U(:,idx);
u = U(:,1);

% This is the MUSIC version of source orientation:
u = u/norm(u);  % For generalized EV |u| is not always 1

if bRAPBeam
    % Find source orientation as EVec of a problem Du = lambda Fu, with
    % D=G, F=S for MAI, D = S, F = T for MPZ,
    % D = E, F = T for MER
    switch(iBeam)
        case BEAMTYPE.IMPZ            % PZ = S/T
            D = Ht * Rm1 * H;
            F = Ht * Rm2 * H;
        case BEAMTYPE.IMAI  %  MAI: G/S
            D = HtH;
            F = Ht * Rm1 * H;
        case BEAMTYPE.IMER  % MER = E/T
            D = Ht * arrE * H;
            F = Ht * Rm2 * H;
        case BEAMTYPE.IRMER  %  rMER: A = Rm1, B = Rm1 * Cavg * Rm1
            error('Not implemented');
    end % switch
    
    [U,D] = eig(D,F);   % Rm1, Rm2, arrE are enforced to be real symmetric
                        % so U, D are always real
    [~,idx] = sort(diag(D),'descend');

    U = U(:,idx);
    uBeam = U(:,1);
    % This is the scalar FS for the voxel:
    h = H*uBeam;

    % RAP-beamformer version of the localizer
    % NOTE that we already have h=Q_prev * h_org, where Q_prev is an
    % out-projector of previously found sources. Therefore we simply use a
    % straight single source expressions here
    switch(iBeam)
        case BEAMTYPE.IMPZ  % MPZ: A = Rm2, B = Rm1
            % PZ = S/T
            P = (h'*Rm1*h)/(h'*Rm2*h);
        case BEAMTYPE.IMAI  %  MAI
            % AI = G/S = 1/S
            P = (h'* h)/(h'*Rm1*h);        
        case BEAMTYPE.IMER  %  MER: A = Rm1 * Rm1, B = Rm1 * Cavg * Rm1
            % SER = E/T
            P = (h'*arrE*h)/(h'*Rm2*h);
        case BEAMTYPE.IRMER  %  rMER: A = Rm1, B = Rm1 * Cavg * Rm1
            error('Not implemented');
    end % switch(iBeam)
    % The final localizer is P/(1-mu)
    P=P/(1-mu(1));
else
    % This is the scalar FS for the voxel:
    h = H*u;

    % SMCMV version of the localizer
    % Form an out-projector of all previously found source FS with this one
    M = length(h);
    Hsrc1 = [Hsrc,h];   % Augment found source set with h
    Q = eye(M) - Hsrc1 * pinv(Hsrc1);   % This is an out-projector

    % Find the new out-projected source subspace Usrc1
    Usrc1 = Q*Usrc; % Now Usrc1 is not orhonormal

    % Calculate the localizer on the out-projected source subspace Usrc
    switch(iBeam)
        case BEAMTYPE.IMPZ  % MPZ: A = Rm2, B = Rm1
            % MPZ = tr(S*T^-1) = tr(Usrc' Rm1 Usrc * (Usrc' Rm2 Usrc)^(-1))
            P = trace(Usrc1'*Rm1*Usrc1*inv(Usrc1'*Rm2*Usrc1))-size(Usrc1,2); %#ok<MINV>
            % Q: Is the above expression valid if Usrc1 is not ONB?
        case BEAMTYPE.IMAI  %  MAI
            % MAI = tr(G S^-1) = tr(Usrc' Rm1 Usrc)^(-1)) for orthonormal Usrc
    %         P = trace(Usrc1'*Usrc1*invSPD(Usrc1'*Rm1*Usrc1));
            % This expression does not need ONB 
            P = trace(Usrc1'*Usrc1*inv(Usrc1'*Rm1*Usrc1))-size(Usrc1,2); %#ok<MINV>        
        case BEAMTYPE.IMER  %  MER: A = Rm1 * Rm1, B = Rm1 * Cavg * Rm1
            P = trace(Usrc1'*arrE*Usrc1*inv(Usrc1'*Rm2*Usrc1)); %#ok<MINV>
        case BEAMTYPE.IRMER  %  rMER: A = Rm1, B = Rm1 * Cavg * Rm1
            error('Not implemented');
    end % switch(iBeam)

    % The final localizer is (1-mu)^(-1)*(P+1)
    P=1/(1-mu(1))/(P+1);
end % if bRAPBeam

