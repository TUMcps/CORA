function res = test_zonoBundle_lift
% test_zonoBundle_lift - unit test function of lift
%
% Syntax:
%    res = test_zonoBundle_lift
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
Z1 = zonotope([1;1], [3 0; 0 2]);
Z2 = zonotope([0;0], [2 2; 2 -2]);
zB = zonoBundle({Z1,Z2});

% should always fail
try
    zB_lift = lift(zB,5,[2,1]);
    % should fail
    resvec(end+1) = false;
catch ME
    resvec(end+1) = strcmp(ME.identifier,'CORA:notDefined');
end

% except no new dimensions are created
zB_lift = lift(zB,2,[2,1]);
zB_proj = project(zB,[2,1]);
resvec(end+1) = isequal(zB_proj,zB_lift);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
