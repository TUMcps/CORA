function res = isZero(P,varargin)
% isZero - Checks if an polytope only represents the origin; if a tolerance
%    is given, it is checked whether the capsule is contained in the ball
%    centered at the origin with radius tolerance
%
% Syntax:  
%    res = isZero(P)
%    res = isZero(P,tol)
%
% Inputs:
%    P - mptPolytope object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    P1 = interval(0,0);
%    res = isZero(P);
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
inputArgsCheck({{P,'att','mptPolytope'};
                {tol,'att','numeric',{'nonnegative','scalar','nonnan'}}});

% empty case
if isempty(P)
    res = false; return
end

% enclose by an interval and compute norm
res = norm(interval(P)) <= tol;

%------------- END OF CODE --------------