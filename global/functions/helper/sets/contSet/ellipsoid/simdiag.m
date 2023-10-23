function [T,D] = simdiag(M1,M2,TOL)
% simdiag - Find T such that
%           T*M1*T' = I and T*M2*T' = D (diagonal)
%
% Syntax:
%    [T,D] = simdiag(M1,M2,TOL)
%
% Inputs:
%    M1 - numerical matrix
%    M2 - numerical matrix
%    TOL - tolerance
%
% Outputs:
%    T - transformation matrix
%    D - result of T*M2*T' (diagonal matrix)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Victor Gassmann
% Written:       06-June-2022 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~all(size(M1)==size(M2))
    throw(CORAerror('CORA:dimensionMismatch',M1,M2));
end
% check if both are symmetric, M1 is pd, and M2 psd (probably too strict of a condition)
if ~isApproxSymmetric(M1,TOL) || ~isApproxSymmetric(M2,TOL) || ...
        min(eig(M1))<TOL || min(eig(M2))<-TOL
    throw(CORAerror('CORA:specialError',...
        'Both matrices need to be symmetric, first matrix needs to be pd, second needs to be psd!'));
end

[U1,S1,~] = svd(M1);
S1_12inv = diag(1./sqrt(diag(S1)));
[U2,~,~] = svd(S1_12inv*U1'*M2*U1*S1_12inv);
T = U2'*S1_12inv*U1';

if nargout > 1
    D = T*M2*T';
end

% ------------------------------ END OF CODE ------------------------------
