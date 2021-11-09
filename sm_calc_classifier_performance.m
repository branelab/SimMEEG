function [perf]=sm_calc_classifier_performance(Hit,Miss,CR,FA)
%function [perf]=calc_classifier_performance(Hit,Miss,CR,FA);
%
% REturns the classifier performance for given Hit, Miss, CR, and FA
% values.
%  INPUT
%     Hit = number of hits
%     Miss = number of misses
%     CR = number of correct rejections
%     FA = number of false alarms
%  OUTPUT 
%   perf. 
%       prevalence = prevalence of # true positives
%       TPR = True Positive Rate = sensitivity of classifier 
%       TNR = True Negative Rate = specificity of classifier 
%       FPR = False Positive Rate  
%       ...
%
% created by A. Herdman, UBC  Sep 1, 2015
%
%   Note: based on equations from https:././en.wikipedia.org./wiki./Sensitivity_and_specificity

perf.prevalence = (Hit + Miss) ./ (Hit + Miss + CR + FA); % prevalence
perf.accuracy = (Hit + CR) ./ (Hit + Miss + CR + FA);    % Accuracy 
perf.TPR = Hit ./ (Hit+Miss);    % sensitivity aka True Positive Rate
perf.TNR = CR ./ (FA + CR);  % True-Negative Rate

perf.FNR = Miss ./ (Hit + Miss); % False Negative Rate
perf.PPV = Hit ./ (Hit + FA);   % Positive Predictive Value aka Precision
perf.FDR = FA ./ (Hit + FA);   % False Discovery Rate

%perf.FPR = FA ./ (FA + CR);  % False-Positive Rate
perf.FPR = 1-perf.TNR;  % False-Positive Rate

perf.FOR = Miss ./ (CR + Miss);  % False Omission Rate
perf.NPV = CR ./ (CR + Miss);    % Negative Predictive Value

perf.PLR = perf.TPR./perf.FPR;   % positive Likelihood Ratio (LR+)
perf.NLR = perf.FNR./perf.TNR;           % Negative Likelihood Ratio (LR-)
perf.DOR = perf.PLR./perf.NLR;           % Diagnostic Odds Ratio
perf.F1 = (2.*Hit) ./ ( (2.*Hit) + FA + Miss); % Harmonic mean of precision./sensitivity
perf.MCC = ( (Hit.*CR)-(FA.*Miss) )  ./ sqrt( (Hit+FA).*(Hit+Miss).*(CR+FA).*(CR+Miss) );     % Mathews correlation coefficient
perf.informedness = perf.TPR + perf.TNR -1; % Informedness
perf.markedness = perf.PPV + perf.NPV -1;                   % markedness


