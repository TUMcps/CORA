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

resvec = [];

% init
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1];
b = 1;
cZ = conZonotope(Z,A,b);

% should always fail
try
    cZ_lift = lift(cZ,5,[2,1]);
    % should fail
    resvec(end+1) = false;
catch ME
    resvec(end+1) = strcmp(ME.identifier,'CORA:notDefined');
end

% except no new dimensions are created
cZ_lift = lift(cZ,2,[2,1]);
cZ_proj = project(cZ,[2,1]);
resvec(end+1) = isequal(cZ_proj,cZ_lift);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
