function sdata = genSurrogateData(data, seed)
%
% SYNTAX:
%
% Generate a surrogate MEG/EEG dataset with the following properties:
%   - it has the same spatial (sensor) covariance matrix as the original
%     real data 
%   - each sensor signal has the same power spectrum as original real data
%   - any synchronization/connectivity that possibly existed in the
%     original real data is destroyed
% This is done by 1) extracting time courses of principal components of the
% original data 2) generating new time courses by randomizing phases of the
% frequency components of each time course (therefore power spectrum of the
% PC is preserved). Technically, those new time courses will no longer be
% completely uncorrelated (one gets correlations ~0.01(, but we ignore
% those. 3) building a new dataset by transforming time courses from PC
% basis to the original (sensor) basis.
%
% NOTE: Continuous (not epoched) original data is used as input. One can
% always generate a single-epoch dataset from multi-epoch data by merging
% the epochs. However if those epochs do not form continuous data the sensor
% spectra of the resulting surrogate may differ from the original due to 
% ringing artifacts at discontinuities.
%
% OBSERVATIONS for 151 x 120000 data:
%   - To get a good match between original and surrogate covariances
%     (Frob err^2 < 0.01) one needs around 10000 samples of data
%   - After phase reshuffle, new PC time courses have correlations around
%     0.01 - 0.02
%   - The surrogate channel spectra generally fit the original ones but
%     some smaller details may be destroyed (likely because we partially
%     changed time correlations between the physical channels)
%
% Input:
%   data        - (M x T) Original data used to generate a surrogate. Here
%                 M is the number of sensor channels, T is the number of
%                 time points
% Optional:
%   seed        - if specified - a seed value for the random generator, to
%                 be used in a call: "rng(seed);"
%
% Output:
%   sdata       - (M x T) Surrogate data
%
% A.Moiseev, BCNI, July 2018
NMANDATORY = 1;

if nargin > NMANDATORY
    rng(seed);
end % if nargin > NMANDATORY

[M, T] = size(data);
dmean = mean(data,2);    % (1 x M) vector of mean values of the channels

% Subtract mean values. Not exactly necessary, but may help to avoid huge
% FFT amps due to discontinuity at the ends if the data is unfiltered
for iCh = 1:M
    data(iCh,:) = data(iCh,:) - dmean(iCh);
end % for iCh = 1:M

R = cov(data');     % (M x M) covariance matrix of the original data
[V,D] = eig(R);     % Columns of V are principal vectors, diag(D) are PC powers
D = diag(D);        %#ok<NASGU> % Now vector of PC powers

% Generate PC time courses (as columns)
sdata = data' * V;  % (T x M); i-th column is a time course of i-th PC
                    % Verified that mutual corr of those is around 1e-13

% Calculate spectrum
sdata = fft(sdata); % (T x M) complex data

% We now need to randomize phases but to ensure that FFT(k) = FFT(-k)*, to
% keep the signal real. The expression is
%                      N
%        X(k) =       sum  x(n)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.
%                     n=1
% The first component X(1) is DC and is always real
iMax = floor(T/2);
bOdd = iMax < T/2;  % Odd number of samples
iMax1 = iMax + 1;
iMax2 = iMax + 2;

for iS = 1:M
    for iT=2:iMax
        sdata(iT,iS) = randphase(sdata(iT,iS));
        sdata(T-iT+2,iS) = conj(sdata(iT,iS));
    end % for i=1:T
    
    % For even T (bOdd = false) element iMax + 1 should always be real, so
    % we have nothing to do. For odd T we need additionally process iMax+1,
    % iMax +2
    
    if bOdd
        sdata(iMax1,iS) = randphase(sdata(iMax1,iS));
        sdata(iMax2,iS) = conj(sdata(iMax1,iS));
    end
end % for iS = 1:M

% Now return back to the time domain
sdata = real(ifft(sdata));  % Still (T x M)

% % Scale the time courses to get exactly the same variances as original PCs
% % NO NEED: it is always perfectly so 
% for iCh = 1:M
%     sdata(:,iCh) = sdata(:,iCh) * sqrt(D(iCh)/var(sdata(:,iCh)));
% end % for iCh = 1:M

% Now back to the sensor base
% sdata = (sdata * V)';   % M x T - seems wrong
sdata = V * sdata';     % M x T

% And restore the original offsets (means)
for iCh = 1:M
    sdata(iCh,:) = sdata(iCh,:) + dmean(iCh);
end % for iCh = 1:M

% ------------------------------------------------------------
function y = randphase(x)
% ------------------------------------------------------------
% y has the same module as x but with random phase
y = sqrt(x*conj(x))*exp(1i*2*pi*rand);

%========================================================================
% TEST PROGRAM
%========================================================================
% clear;
% close all;
% dbstop if error;
% fs = filesep;
% 
% addpath('/opt/matlab/toolbox/ctfds');
% projectBase = '/data/noisedir';
% noiseDs = [projectBase, fs, 'am_20080728_BN_lp_300tr.ds'];
% 
% % Read the noise trial dataset. !!! THE DATASET IS IN FT !!!
% [data , chNames, SR, ds] = ctfDsRead(noiseDs, 'MEG');
% [nSamples, nChannels, nTrials] = size(data);
% data = data * 1e-15;    % Convert to Tesla
% 
% % Make a continuous dataset
% if(nTrials > 1)
%     data = reshape(                                 ...
%                     permute(data,[1,3,2]),          ... % samples x trials x channels
%                     nSamples*nTrials, nChannels   ...
%                     );  % Now data is continuous ntotalsamples x nchannels
%     nSamples = nSamples*nTrials;
% end
% 
% % Verify for uneven number of samples, if necessary
% % data = data(1000:11000,:);   % 10001 sample
% data = data(1001:11000,:);   % 10000 samples
% 
% % Calculate original covariance (data now is T x M)
% R0 = cov(data);
% 
% % Calculate original spectra
% for iCh = 1:nChannels
%     sp = pwelch(data(:,iCh));
%     
%     if iCh == 1
%         sp0 = zeros(length(sp), nChannels);
%     end % if iCh == 1
%     
%     sp0(:,iCh) = sp;
% end % for iCh = 1:nChannels
% 
% % seed = 12345;
% % sdata = genSurrogateData(data', seed);
% sdata = genSurrogateData(data');
% 
% % Surrogate covariance
% R1 = cov(sdata');
% 
% % Check the difference using Frobenious norm
% dR = R1 -R0;
% fprintf('Squared Frob norm of relative difference (R1-R0): %g \n',...
%     trace(dR*dR)/trace(R0*R0));
% 
% iCh = 123;
% 
% % Plot the original and new time courses
% plot([data(1:200,iCh),sdata(iCh,1:200)']);
% legend('original','surrogate');
% 
% % Plot the original and new spectra
% figure;
% loglog([pwelch(data(:,iCh)), pwelch(sdata(iCh,:)')]);
% legend('original','surrogate');
