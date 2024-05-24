function h = sm_batch_sm_bst_noisecov(varargin)
% Anthony Herdman extracted much of the following code from "bst_noisecov.m" in Brainstorm3. It was altered to fit within SimMEEG software platform.
% Looks like Noise and Data Covariances calculated in Brainstorm are based on averaged (evoked) data - This is a different approach than for SPA/SIA/MIA pseudo-z covariance
%
% Here is the copyright information for the code 
%
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2009-2022
%
% Altered by Anthony Herdman (ath) for entering in variables to the calculations

h = varargin{1};  % ath added

data = varargin{end};
dims_data = size(data); % SimMEEG structure = [samps x sens x trial];
nChanAll = dims_data(2); 

%% ===== COMPUTE NOISE COVARIANCE =====
Ntotal   = zeros(nChanAll);
DataCov = zeros(nChanAll);
NoiseCov = zeros(nChanAll);
FourthMoment = zeros(nChanAll);

%% Create DataMat
DataMat = struct('ChannelFlag',[],'F',[],'Leff',[],'Time',[],'nAvg',[]);
%    [DataMat, iTimeBaseline, iTimeCov] = ReadRecordings(); 
%    iTimeBaseline = base_samps;
%    iTimeCov = act_samps;

DataMat.ChannelFlag = 1:dims_data(2); 
DataMat.F = squeeze(nanmean(data,3))';   % struct = [samps x sens]; Brainstorm uses averaged data for Cov calculations. 
DataMat.Leff = 1; DataMat.nAvg = 1; 
DataMat.Time = h.cfg.study.lat_sim;
iTimeBaseline = h.inv_soln(h.current_inv_soln).params.ctrl_samps;
iTimeCov = h.inv_soln(h.current_inv_soln).params.act_samps;
iGoodChan = 1:dims_data(2);
%%
N = length(iTimeCov);

% === Compute average ===
% Average baseline values
Favg = mean(DataMat.F(:,iTimeBaseline), 2);
% === Compute covariance ===
% Remove baseline average
DataMat.F(iGoodChan,iTimeBaseline) = bst_bsxfun(@minus, DataMat.F(iGoodChan,iTimeBaseline), Favg(iGoodChan,1));
%% DataCov: Compute covariance for this file
fileCov  = DataMat.Leff .* (DataMat.F(iGoodChan,iTimeCov)    * DataMat.F(iGoodChan,iTimeCov)'   );
fileCov2 = DataMat.Leff .* (DataMat.F(iGoodChan,iTimeCov).^2 * DataMat.F(iGoodChan,iTimeCov)'.^2);
% Add file covariance to accumulator
DataCov(iGoodChan,iGoodChan)     = DataCov(iGoodChan,iGoodChan)     + fileCov;
FourthMoment(iGoodChan,iGoodChan) = FourthMoment(iGoodChan,iGoodChan) + fileCov2;
Ntotal(iGoodChan,iGoodChan) = Ntotal(iGoodChan,iGoodChan) + N;

% Remove zeros from N matrix
Ntotal(Ntotal <= 1) = 2;
% Divide final matrix by number of samples
DataCov     = DataCov     ./ (Ntotal - 1);
FourthMoment = FourthMoment ./ (Ntotal - 1);
% Display result in the command window
nSamplesTotal = max(Ntotal(:));
% disp(['Number of time samples used for the noise covariance: ' num2str(nSamplesTotal)]);
% Check for NaN values
if (nnz(isnan(DataCov)) > 0)
    error('The output covariance contains NaN values. Please check your recordings and tag the bad channels correctly.');
end

%% add to SimMEEG structure
h.inv_soln(h.current_inv_soln).soln.bst_cov.DataCovMat.NoiseCov = DataCov;
h.inv_soln(h.current_inv_soln).soln.bst_cov.DataCovMat.Comment = 'Data covariance';
h.inv_soln(h.current_inv_soln).soln.bst_cov.DataCovMat.nSamples = [];
h.inv_soln(h.current_inv_soln).soln.bst_cov.DataCovMat.FourthMoment = FourthMoment;
h.inv_soln(h.current_inv_soln).soln.bst_cov.DataCovMat.History(1) = {datetime};
h.inv_soln(h.current_inv_soln).soln.bst_cov.DataCovMat.History(2) = {'compute'};
h.inv_soln(h.current_inv_soln).soln.bst_cov.DataCovMat.History(3) = {'Computed using SimMEEG'};

%% NoiseCov: Compute covariance for this file
fileCov  = DataMat.Leff .* (DataMat.F(iGoodChan,iTimeBaseline)    * DataMat.F(iGoodChan,iTimeBaseline)'   );
fileCov2 = DataMat.Leff .* (DataMat.F(iGoodChan,iTimeBaseline).^2 * DataMat.F(iGoodChan,iTimeBaseline)'.^2);
% Add file covariance to accumulator
NoiseCov(iGoodChan,iGoodChan)     = NoiseCov(iGoodChan,iGoodChan)     + fileCov;
FourthMoment(iGoodChan,iGoodChan) = FourthMoment(iGoodChan,iGoodChan) + fileCov2;
Ntotal(iGoodChan,iGoodChan) = Ntotal(iGoodChan,iGoodChan) + N;

% Remove zeros from N matrix
Ntotal(Ntotal <= 1) = 2;
% Divide final matrix by number of samples
NoiseCov     = NoiseCov     ./ (Ntotal - 1);
FourthMoment = FourthMoment ./ (Ntotal - 1);
% Display result in the command window
nSamplesTotal = max(Ntotal(:));
% disp(['Number of time samples used for the noise covariance: ' num2str(nSamplesTotal)]);
% Check for NaN values
if (nnz(isnan(NoiseCov)) > 0)
    error('The output covariance contains NaN values. Please check your recordings and tag the bad channels correctly.');
end

%% add to SimMEEG structure
h.inv_soln(h.current_inv_soln).soln.bst_cov.NoiseCovMat.NoiseCov = NoiseCov;
h.inv_soln(h.current_inv_soln).soln.bst_cov.NoiseCovMat.Comment = 'Data covariance';
h.inv_soln(h.current_inv_soln).soln.bst_cov.NoiseCovMat.nSamples = [];
h.inv_soln(h.current_inv_soln).soln.bst_cov.NoiseCovMat.FourthMoment = FourthMoment;
h.inv_soln(h.current_inv_soln).soln.bst_cov.NoiseCovMat.History(1) = {datetime};
h.inv_soln(h.current_inv_soln).soln.bst_cov.NoiseCovMat.History(2) = {'compute'};
h.inv_soln(h.current_inv_soln).soln.bst_cov.NoiseCovMat.History(3) = {'Computed using SimMEEG'};

