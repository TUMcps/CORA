function res = test_ellipsoid_lift
% test_ellipsoid_lift - unit test function of lift
%
% Syntax:
%    res = test_ellipsoid_lift
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       19-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% init
E = ellipsoid([1 0;0 2],[0; 1]);

% should always fail
try
    E_lift = lift(E,5,[2,1]);
    % should fail
    resvec(end+1) = false;
catch ME
    resvec(end+1) = strcmp(ME.identifier,'CORA:notDefined');
end

% except no new dimensions are created
E_lift = lift(E,2,[2,1]);
E_proj = project(E,[2,1]);
resvec(end+1) = isequal(E_proj,E_lift);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
