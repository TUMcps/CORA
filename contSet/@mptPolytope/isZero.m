function res = isZero(P,varargin)
% isZero - Checks if an polytope only represents the origin; if a tolerance
%    is given, it is checked whether the polytope is contained in the ball
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
%    P = mptPolytope([1 1; -2 1; 0 -1],[0.01;0.01;0.01]);
%    isZero(P)
%    isZero(P,0.05)
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
res = norm_(interval(P),2) <= tol;

%------------- END OF CODE --------------