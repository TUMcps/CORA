function res = isZero(C,varargin)
% isZero - Checks if a capsule only represents the origin; if a tolerance
%    is given, it is checked whether the capsule is contained in the ball
%    centered at the origin with radius tolerance
%
% Syntax:  
%    res = isZero(C)
%    res = isZero(C,tol)
%
% Inputs:
%    C - capsule object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    C1 = capsule([0;0]);
%    res = isZero(C1)
%
%    C2 = capsule([0;0],[1;1],0);
%    res = isZero(C2)
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
inputArgsCheck({{C,'att','capsule'};
                {tol,'att','numeric',{'nonnegative','scalar','nonnan'}}});

% empty case
if isempty(C)
    res = false; return
end

if tol == 0
    % no tolerance -> quicker check
    res = all(center(C) == 0) && all(C.g == 0) && C.r == 0;

else
    % including tolerance...
    res = true;

    % two versions for generator addition
    if vecnorm(center(C) + C.g) + C.r > tol
        res = false;
    elseif vecnorm(center(C) - C.g) + C.r > tol
        res = false;
    end

end

%------------- END OF CODE --------------