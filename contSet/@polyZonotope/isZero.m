function res = isZero(pZ,TOL)
% isZero - Checks if a polynomial zonotope is all zero
%
% Syntax:  
%    res = isZero(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%    TOL- (optional) tolerance
%
% Outputs:
%    res - 1 if set is all zero, 0 otherwise
%
% Example: 
%    pZ = polyZonotope([0;0],[1,1],zeros(2,0),[1,2;3,4;5,6]);
%    res = isZero(pZ);    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope, isPolytope

% Author:       Victor Gassmann
% Written:      12-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if ~exist('TOL','var')
    TOL = eps;
end
res = all(abs([pZ.c,pZ.G,pZ.Grest])<=TOL,2);
%------------- END CODE --------------