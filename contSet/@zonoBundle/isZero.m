function res = isZero(zB,varargin)
% isZero - Checks if a zonotope bundle only represents the origin; if a
%    tolerance is given, it is checked whether the capsule is contained in
%    the ball centered at the origin with radius tolerance
%
% Syntax:  
%    res = isZero(zB)
%    res = isZero(zB,tol)
%
% Inputs:
%    zB - zonoBundle object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    zB
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      17-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default value
tol = setDefaultValues({0},varargin);

% parse input arguments
inputArgsCheck({{zB,'att','zonoBundle'};
                {tol,'att','numeric',{'nonnegative','scalar','nonnan'}}});

% empty case
if isemptyobject(zB)
    res = false; return
end

% convert to constrained zonotope
res = isZero(conZonotope(zB),tol);

%------------- END OF CODE --------------