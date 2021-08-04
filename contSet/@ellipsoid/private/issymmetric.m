function res = issymmetric(Q,TOL)
% issymmetric - Checks whether Q matrix is symmetric
%
% Syntax:  
%    res = issymmetric(Q,TOL) Computes whether Q is symmetric using TOL
% Inputs:
%    Q - matrix to check
%    TOL - used tolerance
%
% Outputs:
%    res - logical indicating whether Q is symmetric (within tolerance)
%
% Example: 
%    E = ellipsoid([1,0;0,1/2],[1;1]);
%    res = issymmetric(E.Q,E.TOL);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  15-October-2019
% Last revision:---
%------------- BEGIN CODE --------------
res = all(all(withinTol(triu(Q),tril(Q)',ellipsoid().TOL)));
%------------- END OF CODE --------------