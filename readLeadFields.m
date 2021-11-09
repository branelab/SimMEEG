function [arrH, lstFlag, lstSensorNames, samHeader, lstVoxels] = readLeadFields(strFile, lstVoxels)
%   SYNTAX:
%   [arrH, lstFlag, lstSensorNames, samHeader, lstVoxels] = readLeadFields(strFile, lstVoxels)
%
%  * The lead fields data file (leadFields.dat) has the following structure:
%  * 	Name				Length			Comment
%  * 	----				------			-------
%  * 	ID  				8			Data identifier (string): LEADFLDS
%  *	SAM_HDR_v2  		768			Standard header common for SAM files
%  *	MEG sensor names	16*M		M - number of primary sensors used
%  *	[flag + pos + lead fields]	(4+8*(3+3*M))*V	V - total number of voxels
%  *
%  * NOTES:
%  * 	- data is written as LITTLE ENDIAN (not like other SAM files)
%  * 	- number of physical channels (M) may be found at offset 260 from the start of
%  * 	  SAM header, or 268 (0x10c) from the start of the file. This is a 4 byte integer
%  * 	- number of voxels (V) may be found at offset 264 from the start of the header,
%  * 	  or 272 (0x110) from the start of the file
%  * 	- for each voxel, [flag + pos + lead fields] record consists of inside head flag (4 bytes),
%  * 	  voxel coordinates in m (3 doubles) followed by 3 lead field vectors for X, Y, Z directions, respectively
%  * 	  Each lead field vector is an array of M doubles
%  *
%  * 	Thus, the size of the output data file should be equal to:
%  * 	(4+8*(3+3*M))*V + 16 M + 776 
%   Input:
%       strFile     full path to the data file
%       lstVoxels   list of voxel locations for the lead fields (METERS!)
%                   (is used to verify that fields DO correspond to the voxels)
%                   nVoxels x 3. The ordering of voxels should be in
%                   accordance with C/C++ array storage, that is the voxel
%                   at postion ix, iy, iz in the grid is stored at the
%                   index idx = (ix-1)*ny*nz + (iy-1)*nz + iz
%                   If EMPTY ([]) - then comparison of the grid voxels is not
%                   performed, and actual voxel positions stored in the
%                   file are returned
%   Output:
%       arrH:       nVoxels x 3 x nSensors array of LCMV lead fields
%       lstFlag:    (nVoxwls). If flag = 1, the voxel is inside the head, 0
%                   otherwise
%       lstSensorNames: nVoxels x 5 character array of names of the physical
%                   sensors being used
%   samHeader       SAM file header v.2
%   lstVoxels       (nVoxels x 3) list of voxel positions in (m). If this
%                   one was supplied on input - the same list is returned,
%                   otherwise list stored in the file is returned
%
% Update Aug 2011:
%   For voxels that have FS identically equal to 0 - such as those beyond
%   the sensor shell or too close to the center of the head if a spherical
%   model is used - the flag is set to 0, irrespective to the fact that the
%   voxel may be inside the head.
%
% Update Aug 2019:
%   If voxel list is not supplied for verification, the list stored inside
%   the file is returned
%
% A. Moiseev, DSRF, Sep 2008

fp = fopen(strFile,'rb');

% Was:
% fseek(fp, 268, 'bof');
% M = fread(fp, 1, 'int32');  % Number of channels
% V = fread(fp, 1, 'int32');  % Number of voxels
% Now, using readSAMHeaderV2():
fseek(fp, 8, 'bof');    % Skip file tag
samHeader = readSAMHeaderV2(fp);
M = samHeader.NumChans;
V = samHeader.NumWeights;

if ~isempty(lstVoxels)
    bVerifyGrid = true;
    nVoxels = size(lstVoxels,1);

    if nVoxels ~= V
        fprintf('ERROR: Number of voxels does not match the number in the leadFields file\n');
        fclose(fp);
        return;
    end
else
    bVerifyGrid = false;
    lstVoxels = NaN(V, 3);
end % if ~isempty(lstVoxels)

fseek(fp, 776, 'bof');  % Move to the channel names
lstSensorNames=char(zeros(M,5));

for i=1:M
    s = fread(fp, 16, 'schar');     % full name: MLT32-1234
    lstSensorNames(i,:)=s(1:5);     % short name: MLT32
end

lstFlag=zeros(V,1,'int32');
arrH=zeros(V,3,M);

% Read the flags and the forward solutions

for i=1:V
    lstFlag(i)=fread(fp, 1, 'int32');
    
    % Read voxel position and compare with lstVoxels
    dd = fread(fp, 3, 'float64');
    
    if bVerifyGrid
        if norm(dd-lstVoxels(i,:)') > 1e-4      % lstVoxels(i,:) is a row vector      
            fprintf('ERROR: Voxels locations do not match those in the leadFields file\n');
            fclose(fp);
            return;
        end     % if
    else
        lstVoxels(i,:) = dd';
    end % if bVerifyGrid
    
    for iH=1:3
        h = fread(fp, double(M), 'float64');
        arrH(i,iH,:) = h;
        
        % For "ivalid" voxels such as beyound the sensor or too close to
        % the center of the head for the spherical model, the lead fields
        % program sets FS to exact zero. The flag may not necessarily be
        % zero.        
        if h'*h == 0
            lstFlag(i) = 0;
        end
    end     % for iH
end     % for i

fclose(fp);
