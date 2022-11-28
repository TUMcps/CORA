function res = eq(E1,E2)
% eq - Overloaded '==' operator for the comparison of ellipsoids
%
% Syntax:  
%    res = eq(E1,E2)
%
% Inputs:
%    E1 - ellipsoid object 
%    E2 - ellipsoid object 
%
% Outputs:
%    res - true/false
%
% Example: 
%    E1 = ellipsoid([1 0; 0 1]);
%    E2 = ellipsoid(E1.Q,E1.q);
%    res = E1==E2;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Victor Gassmann
% Written:      14-October-2019 
% Last update:  16-March-2021 (relative TOL)
%               04-July-2022 (VG: input checks)
% Last revision:---

%------------- BEGIN CODE --------------

% check input arguments
inputArgsCheck({{E1,'att','ellipsoid','scalar'};
                {E2,'att','ellipsoid','scalar'}});

% check if dimensions are equal
if dim(E1) ~= dim(E2)
    res = false;
    return;
end

% check for emptyness
if (isempty(E1) && ~isempty(E2)) || (~isempty(E1) && isempty(E2))
    res = false;
    return;
elseif isempty(E1) && isempty(E2)
    res = true;
    return;
end
    
% set tolerance for numerical comparsion
TOL = min(E1.TOL,E2.TOL);

% compare shape matrix and center numerically
res = all(all(withinTol(E1.Q,E2.Q,TOL))) && all(withinTol(E1.q,E2.q,TOL));

%------------- END OF CODE --------------