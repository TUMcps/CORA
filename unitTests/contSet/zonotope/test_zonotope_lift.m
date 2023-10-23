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

resvec = [];

% init
Z = zonotope([1;1], [3 0; 0 2]);

% should always fail
try
    Z_lift = lift(Z,5,[2,1]);
    % should fail
    resvec(end+1) = false;
catch ME
    resvec(end+1) = strcmp(ME.identifier,'CORA:notDefined');
end

% except no new dimensions are created
Z_lift = lift(Z,2,[2,1]);
Z_proj = project(Z,[2,1]);
resvec(end+1) = isequal(Z_proj,Z_lift);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
