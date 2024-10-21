function res = test_capsule_lift
% test_capsule_lift - unit test function of lift
%
% Syntax:
%    res = test_capsule_lift
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
% See also: none

% Authors:       Tobias Ladner
% Written:       19-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init
C = capsule([2; -1], [1; 2], 0.5);

% should always fail
try
    C_lift = lift(C,5,[2,1]);
    % should fail
    assert(false);
catch ME
    if ~strcmp(ME.identifier,'CORA:notDefined')
        rethrow(ME)
    end
end

% except no new dimensions are created
C_lift = lift(C,2,[2,1]);
C_proj = project(C,[2,1]);
assert(isequal(C_proj,C_lift));

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
