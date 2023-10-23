function res = test_polyZonotope_lift
% test_polyZonotope_lift - unit test function of lift
%
% Syntax:
%    res = test_polyZonotope_lift
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
pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0.5;0],[1 0 3;0 1 1]);

% should always fail
try
    pZ_lift = lift(pZ,5,[2,1]);
    % should fail
    resvec(end+1) = false;
catch ME
    resvec(end+1) = strcmp(ME.identifier,'CORA:notDefined');
end

% except no new dimensions are created
pZ_lift = lift(pZ,2,[2,1]);
pZ_proj = project(pZ,[2,1]);
resvec(end+1) = isequal(pZ_proj,pZ_lift);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
