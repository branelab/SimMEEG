function [peaks, arrIdx3D, idxMax, lstBorderMax] = svlPeak(data, dims, alpha, gap, bMatlab)
%
%   SYNTAX:
%       [peaks, arrIdx3D, idxMax, lstBorderMax] = svlPeak(data, dims, alpha, gap, bMatlab)
%
%   3D image peak finder analogous to the CTF's svlPeak routine.
%   Generalized to multi-component (vector) images by searching for the
%   maximum of the norm of the vector (or rather a sum of squares of the
%   components). The function only finds maxima. For minima - apply this routine
%   to (-data).
%   MIND TO SET VOXELS THAT SHOULD NOT BE INCLUDED INTO THRESHOLD
%   ESTIMATION TO NaN !!!
%
%   Input:
%       data:       (nVoxels x 1) vector of image values, or (nVoxels x d)
%                   list of vector values for each voxel. The data
%                   represents a 3D rectangular region, with C++ array
%                   oredring of values (that is, Z-index IS THE FASTEST)
%                   when bMatlab flag is set to 0, and with X being the
%                   fastest one when it is 1.
%       dims:       1 x 3 list of ROI dimensions in voxels: [nx ny nz]
%       alpha:      the p-value threshold for the peaks. Only peaks higher
%                   than (1 - alpha) quantile boundary will be taken into
%                   account.
%       gap:        minimal distance in voxels between maxima. Min gap is
%                   2.
%       bMatlab:    Set 1 if 3D data is stored in memory according to the
%                   Matlab conventions (the 1st index is the fastest one);
%                   set to 0 if it is stored in C/C++ convention (the 3d
%                   index is the fastest)
%                   !!! Mind that data read by readSvlImage.m is stored !!!
%                   !!! in MATLAB convention. Set bMatlab = TRUE !!!!!!!!!!
%                   (although in the SVL file C/C++ convension is used)
%                   !!! Mind also that if lstFlag is used to mask pixels
%                   outside the head - lstFlag uses C/C++ ordering. Means
%                   that in this case either data or lstFlag should be
%                   reordered !!!
%   Ouput:
%       peaks:      nPeaks x 1 - the list of the max values of the data. For
%                   vector data, L2-norm of the vector is returned
%       arrIdx3D:   nPeaks x 3 - integer vector of voxel coordinates [ix iy
%                   iz] for the maxima
%       idxMax:     nPeaks x 1 - the linear indecies of the maxima into the
%                   data array
%       lstBorderMax:
%                   nBorderMax x 3 - a list of voxel coordinates [ix iy iz]
%                   for the ROI border maxima (if any)
%
% A. Moiseev, DSRF, Oct 2008

% Have something to return even when no peaks are found:
peaks = zeros(0,1);
arrIdx3D = zeros(0,3);
idxMax = zeros(0,1);
lstBorderMax = zeros(0,3);

nx = dims(1); ny = dims(2); nz = dims(3);
nVoxels = size(data, 1);

if nx*ny*nz ~= nVoxels
    fprintf(2,'svlPeak ERROR: supplied dimensions (dims) are not consistent with the data\n');
    return;
end

if (alpha < 0) || (alpha > 1)
    fprintf(2,'svlPeak ERROR: alpha should belong to [0, 1] interval\n');
    return;
end

if gap < 2
    fprintf(2,'svlPeak ERROR: gap < 2 specified\n');
    return;
end

gp = gap - 1; % This determines the box boundaries around each voxel for peak finding
bVectorData = 0;    % Flag that the data is vector-valued

if size(data, 2) > 1    % vector data is supplied
    tmp = zeros(nVoxels,1);
    bVectorData = 1;
    
    for i=1:nVoxels
        tmp(i) = data(i,:)*data(i,:)';
    end
    
    data = tmp; % Now data here is a scalar image; original function argument
                % "data" should remain unchanged
end     % if size(data, 2) > 1

% Define the threshold. Mind that Analyze images may have NaNs as voxel
% values
notNaN = ~isnan(data);  % Mask for numeric voxels
th = quantile(data(notNaN), 1 - alpha);

aboveIdx = 1:nVoxels;
aboveIdx = aboveIdx(notNaN & (data > th)); % List of voxels above the threshold

isMax = logical(zeros(nVoxels,1));
isMax(aboveIdx) = true;

lstBorderMax = zeros(0,3);

for i  = 1:length(aboveIdx)
    iV = aboveIdx(i);
    
    if isMax(iV) == false
        continue;   % skip voxels already marked as not maxima
    end
    
    bSmallerFound = 0;  % Flag that a smaller voxel in the vicinity already found 
	[ix, iy, iz] = idx1Dto3D(iV, dims, bMatlab);
    
    % QQQQQQ debug
    if (ix < 1) || (iy < 1) || (iz < 1) || (ix > nx) || (iy > ny) || (iz > nz)
        fprintf(2,'ERROR: index out of bounds\n');
        return;
    end % QQQQQQQQQQQQ
	
    bNextI = 0; % Flag to break all internal loops and continue with the next i
    
    for jx = max(1, ix - gp) : min(nx, ix + gp)
  		for jy = max(1, iy - gp) : min(ny, iy + gp)
    		for jz = max(1, iz - gp) : min(nz, iz + gp)
                % Find a linear idx of current voxel in the box
                jV = idx3Dto1D([jx, jy, jz], dims, bMatlab);
                
                if jV == iV
                    continue;
                end
                
                if isnan(data(jV))
                    continue;
                end
                
                if data(jV) < data(iV)
                    isMax(jV) = false;
                    bSmallerFound = 1;
                else
                    if (data(jV) == data(iV)) && (bSmallerFound ~= 0)
                        isMax(jV) = false;
                    else
                        isMax(iV) = false;
                        bNextI = 1;
                        break;  % Get out of the 3-level loop
                    end % if (data(jV) == data(iV))                                    
                end % if data(jV) < data(iV)                             
            end % for jz
            
            if bNextI == 1, break, end
        end % for jy
        
        if bNextI == 1, break, end
    end % for jx
    
    if bNextI == 1, continue, end
    
    % We get here iff iV is a maximum. Check for boundaries
    if(ix == 1) || (iy == 1) || (iz == 1) || (ix == nx) || (iy == ny) || (iz == nz)
        lstBorderMax(size(lstBorderMax,1)+1,:)= [ix iy iz];
  		isMax(iV) = false;
    end
end     % i  = 1:length(aboveIdx)

idxMax = [1:nVoxels];
idxMax = idxMax(isMax);

nPeaks = length(idxMax);

if nPeaks == 0
    return;
end

peaks = data(idxMax);

[peaks, perm] = sort(peaks, 1, 'descend');
idxMax = idxMax(perm); % Now idxMax is also sorted

arrIdx3D = zeros(nPeaks, 3);

for i = 1:nPeaks
    [ix iy iz] = idx1Dto3D(idxMax(i), dims, bMatlab);
    arrIdx3D(i,:) = [ix iy iz];
end % i = 1:nPeaks

if bVectorData ~= 0
    peaks = sqrt(peaks);
end

if size(lstBorderMax,1) > 0
    fprintf('svlPeak[] WARNING: one or more maxima on the borders were discarded\n');
end

% All done
return;

% ========================================================

function [ix, iy, iz] = idx1Dto3D(iV, dims, bMatlab)
% Convert one-dimentional index to a vector of 3 coordinate indecies
nx = dims(1); ny = dims(2); nz = dims(3);

if bMatlab == 0
    ix = floor((iV - 1)/(ny*nz)) + 1;
    iy = floor((iV - 1 - (ix - 1)*ny*nz)/nz) + 1;
    iz = (iV - 1 - (ix - 1)*ny*nz) - (iy - 1)*nz + 1;
else
    iz = floor((iV - 1)/(ny*nx)) + 1;
    iy = floor((iV - 1 - (iz - 1)*ny*nx)/nx) + 1;
    ix = (iV - 1 - (iz - 1)*ny*nx) - (iy - 1)*nx + 1;
end

function iV = idx3Dto1D(idx3D, dims, bMatlab)
% Convert a vector of 3 coordinate indecies to a one-dimentional index 
nx = dims(1); ny = dims(2); nz = dims(3);
ix = idx3D(1); iy = idx3D(2); iz = idx3D(3);

if bMatlab == 0
    iV = (ix - 1)*ny*nz + (iy - 1)*nz + iz;
else
    iV = (iz - 1)*ny*nx + (iy - 1)*nx + ix;
end
