function res = isZero(cZ,varargin)
% isZero - Checks if a constrained zonotope only represents the origin; if
%    a tolerance is given, it is checked whether the capsule is contained
%    in the ball centered at the origin with radius tolerance
%
% Syntax:  
%    res = isZero(cZ)
%    res = isZero(cZ,tol)
%
% Inputs:
%    cZ - conZonotope object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    cZ1 = conZonotope([0;0]);
%    res = isZero(cZ1);
%
%    cZ2 = conZonotope([0.02;-0.01],[0.03 0; 0 0.02]);
%    tol = 0.05;
%    res = isZero(cZ2,tol);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      16-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default value
tol = setDefaultValues({0},varargin);

% parse input arguments
inputArgsCheck({{cZ,'att','conZonotope'};
                {tol,'att','numeric',{'nonnegative','scalar','nonnan'}}});

% empty case
if isemptyobject(cZ)
    res = false; return
end

% enclose by an interval and compute norm
res = norm(interval(cZ)) <= tol;

%------------- END OF CODE --------------
