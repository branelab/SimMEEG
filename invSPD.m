function [Am1, bException] = invSPD(A)
%
% SYNTAX:
%   Am1 = invSPD(A)
% 
% This function finds an inverse of a SYMMETRIC POSITIVE DEFINITE matrix.
% Supposed to be more accurate than using inv(), because the 
% latter does not exploit the special form of A
%
% Input:
%   A           (n x n) SPD matrix
% Ouput:
%   Am1         A^(-1)
%   bException  true if fall back to inv() was necessary
%
% A. Moiseev, DSRF, July 2011

bException = false;

try
    R=chol(A);      % Find upper triangular R such that R' * R = A;
    I=eye(size(A)); 
    Am1=R\(R'\I);    % NOTE: A\B = inv(A)*B
catch
    bException = true;
    Am1=inv(A);
    warning('Falled back to inv(A)');
%    display(A);
end

