function res = isZero(I,varargin)
% isZero - Checks if an interval only represents the origin; if a
%    tolerance is given, it is checked whether the capsule is contained in
%    the ball centered at the origin with radius tolerance
%
% Syntax:  
%    res = isZero(I)
%    res = isZero(I,tol)
%
% Inputs:
%    I - interval object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval(0,0);
%    res = isZero(I1);
%
%    I2 = interval([0.01;-0.01],[0.02;0.01]);
%    tol = 0.05;
%    res = isZero(I2,tol);
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
inputArgsCheck({{I,'att','interval'};
                {tol,'att','numeric',{'nonnegative','scalar','nonnan'}}});

% empty case
if isempty(I)
    res = false; return
end

% enclose by an interval and compute norm
res = norm(I) <= tol;

%------------- END OF CODE --------------