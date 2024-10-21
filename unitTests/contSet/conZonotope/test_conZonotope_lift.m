function res = test_conZonotope_lift
% test_conZonotope_lift - unit test function of lift
%
% Syntax:
%    res = test_conZonotope_lift
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
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZ = conZonotope(Z,A,b);

% should always fail
try
    cZ_lift = lift(cZ,5,[2,1]);
    % should fail
    assert(false);
catch ME
    if ~strcmp(ME.identifier,'CORA:notDefined')
        rethrow(ME)
    end
end

% except no new dimensions are created
cZ_lift = lift(cZ,2,[2,1]);
cZ_proj = project(cZ,[2,1]);
assert(isequal(cZ_proj,cZ_lift));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
