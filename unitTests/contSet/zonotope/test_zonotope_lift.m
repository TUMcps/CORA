function res = test_zonotope_lift
% test_zonotope_lift - unit test function of lift
%
% Syntax:
%    res = test_zonotope_lift
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
Z = zonotope([1;1], [3 0; 0 2]);

% should always fail
try
    Z_lift = lift(Z,5,[2,1]);
    % should fail
    assert(false);
catch ME
    if ~strcmp(ME.identifier,'CORA:notDefined')
        rethrow(ME)
    end
end

% except no new dimensions are created
Z_lift = lift(Z,2,[2,1]);
Z_proj = project(Z,[2,1]);
assert(isequal(Z_proj,Z_lift));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
