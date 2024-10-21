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

% init
E = ellipsoid([1 0;0 2],[0; 1]);

% should always fail
try
    E_lift = lift(E,5,[2,1]);
    % should fail
    assert(false);
catch ME
    if ~strcmp(ME.identifier,'CORA:notDefined')
        rethrow(ME)
    end
end

% except no new dimensions are created
E_lift = lift(E,2,[2,1]);
E_proj = project(E,[2,1]);
assert(isequal(E_proj,E_lift));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
