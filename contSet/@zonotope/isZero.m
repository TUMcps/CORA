function res = isZero(Z,varargin)
% isZero - Checks if a zonotope only represents the origin; if a tolerance
%    is given, it is checked whether the capsule is contained in the ball
%    centered at the origin with radius tolerance
%
% Syntax:  
%    res = isZero(Z)
%    res = isZero(Z,tol)
%
% Inputs:
%    Z - zonotope object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    Z = zonotope([0;0.01],[0.01 0; 0.02 -0.01]);
%    tol = 0.05;
%    res = isZero(Z,tol);
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
inputArgsCheck({{Z,'att','zonotope'};
                {tol,'att','numeric',{'nonnegative','scalar','nonnan'}}});

% empty case
if isempty(Z)
    res = false; return
end

% enclose by an interval and compute norm
res = norm(interval(Z)) <= tol;

%------------- END OF CODE --------------
