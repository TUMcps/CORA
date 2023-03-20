function res = isZero(E,varargin)
% isZero - Checks if an ellipsoid only represents the origin; if a
%    tolerance is given, it is checked whether the capsule is contained in
%    the ball centered at the origin with radius tolerance
%
% Syntax:  
%    res = isZero(E)
%    res = isZero(E,tol)
%
% Inputs:
%    E - ellipsoid object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    E = ellipsoid(zeros(3),zeros(3,1));
%    res = isZero(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      17-March-2023
% Last update:  20-March-2023 (VG: allow class arrays)
% Last revision:---

%------------- BEGIN CODE --------------

% default value
tol = setDefaultValues({0},varargin);

% parse input arguments
inputArgsCheck({{E,'att','ellipsoid'};
                {tol,'att','numeric',{'nonnegative','scalar','nonnan'}}});

% for each ellipsoid in E
res = false(size(E));
for i=1:length(E)
    % empty case
    if ~isempty(E(i))
        % compute enclosing ball for ellipsoid and compare to tolerance
        res = vecnorm(E(i).q) + radius(E(i)) <= tol;
    end
end
%------------- END OF CODE --------------
