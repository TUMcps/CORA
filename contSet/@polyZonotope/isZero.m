function res = isZero(pZ,varargin)
% isZero - Checks if a polynomial zonotope only represents the origin
%
% Syntax:  
%    res = isZero(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    pZ = polyZonotope([0;0],[1,2;0,0],[],[1,2;3,4;5,6]);
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

% default value
tol = setDefaultValues({0},varargin);

% parse input arguments
inputArgsCheck({{pZ,'att','polyZonotope'};
                {tol,'att','numeric',{'nonnegative','scalar','nonnan'}}});

tmp = abs([pZ.c,pZ.G,pZ.Grest]);
res = all(( withinTol(tmp,0,tol) ),2);

%------------- END OF CODE --------------