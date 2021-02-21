function samHeader = readSAMHeaderV2(fp)
%
% SYNTAX:
%       samHeader = readSAMHeaderV2(fp)
%
% Reads the SAM header (v.2) structure from the binary file. It is assumed
% that file is already opened for reading in proper mode (either big-endian
% or little-endian), so no additional mode specification is done here. File
% stays open.
%
% Input:
%   fp              The file pointer
% Output:
%   samHeader       SAM svl image header v.2 (usually obtained from the lead fields
%                   file)
%
% A. Moiseev, DSRF, Apr 2009

samHeader.Version = cast(fread(fp, 1, 'int32'), 'int32');    % Version
samHeader.SetName = cast(fread(fp, 256, 'schar'),'uint8');    % Name of parent dataset, 256 bytes
samHeader.NumChans = cast(fread(fp, 1, 'int32'), 'int32');    % as is
samHeader.NumWeights = cast(fread(fp, 1, 'int32'), 'int32');    % 0 for svl files
samHeader.pad_bytes1 = cast(fread(fp, 1, 'int32'), 'int32');    % as is
samHeader.XStart = fread(fp, 1, 'double');    % as is, in m
samHeader.XEnd = fread(fp, 1, 'double');    % as is, in m
samHeader.YStart = fread(fp, 1, 'double');    % as is, in m
samHeader.YEnd = fread(fp, 1, 'double');    % as is, in m
samHeader.ZStart = fread(fp, 1, 'double');    % as is, in m
samHeader.ZEnd = fread(fp, 1, 'double');    % as is, in m
samHeader.StepSize = fread(fp, 1, 'double');    % (m)
samHeader.HPFreq = fread(fp, 1, 'double');    % (Hz)
samHeader.LPFreq = fread(fp, 1, 'double');    % (Hz)
samHeader.BWFreq = fread(fp, 1, 'double');    % (Hz)
samHeader.MeanNoise = fread(fp, 1, 'double');    % (T)
samHeader.MriName = cast(fread(fp, 256, 'schar'),'uint8');    % 256 bytes
samHeader.Nasion = cast(fread(fp, 3, 'int32'), 'int32');    % MRI voxel index for nasion, 3 values
samHeader.RightPA = cast(fread(fp, 3, 'int32'), 'int32');    % MRI voxel index for right, 3 values
samHeader.LeftPA = cast(fread(fp, 3, 'int32'), 'int32');    % MRI voxel index for left, 3 values
samHeader.SAMType = cast(fread(fp, 1, 'int32'), 'uint8');    % SAM_TYPE_IMAGE (0), SAM_TYPE_WT_ARRAY (1),  SAM_TYPE_WT_LIST (2)
samHeader.SAMUnit = cast(fread(fp, 1, 'int32'), 'int32');    % SAM_UNIT_COEFF (0), SAM_UNIT_MOMENT (1), ... see docs
samHeader.pad_bytes2 = cast(fread(fp, 1, 'int32'), 'int32');    % as is
samHeader.MEGNasion = fread(fp, 3, 'double');    % MEG dewar coordes for nasion, 3 values
samHeader.MEGRightPA = fread(fp, 3, 'double');    % MEG dewar coordes for right, 3 values
samHeader.MEGLeftPA = fread(fp, 3, 'double');    % MEG dewar coordes for left, 3 values
samHeader.SAMUnitName = cast(fread(fp, 32, 'schar'),'uint8');    % as is
