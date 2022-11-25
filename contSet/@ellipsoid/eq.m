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
%    res - result of comparison
%
% Example: 
%    E1=ellipsoid([1 0; 0 1]);
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
% Last revision:---

%------------- BEGIN CODE --------------
if ~isa(E1,'ellipsoid') || ~isa(E2,'ellipsoid')
    error('Both objects need to be ellipsoids!');
end
if E1.dim ~= E2.dim
    res = false;
    return;
end
if (isempty(E1) && ~isempty(E2)) || (~isempty(E1) && isempty(E2))
    res = false;
    return;
elseif isempty(E1) && isempty(E2)
    res = true;
    return;
end
    
TOL = min(E1.TOL,E2.TOL);

res = all(all(withinTol(E1.Q,E2.Q,TOL))) && all(withinTol(E1.q,E2.q,TOL));
%------------- END OF CODE --------------