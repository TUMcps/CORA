function res = isApproxSymmetric(Q,TOL)
% isApproxSymmetric - Checks if a shape matrix is symmetric (within
%    tolerance)
%
% Syntax:
%    res = isApproxSymmetric(Q,TOL)
%
% Inputs:
%    Q - square matrix
%    TOL - tolerance
%
% Outputs:
%    res - true/false indicating whether Q is symmetric (within tolerance)
%
% Example: 
%    E = ellipsoid([1,0;0,1/2],[1;1]);
%    res = isApproxSymmetric(E.Q,E.TOL);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   15-October-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% take default value for tolerance if none given
if ~exist('TOL','var')
    TOL = ellipsoid.empty(1).TOL;
end

res = all(all(withinTol(triu(Q),tril(Q)',TOL)));

% ------------------------------ END OF CODE ------------------------------
