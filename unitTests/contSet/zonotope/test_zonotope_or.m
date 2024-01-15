function res = test_zonotope_or
% test_zonotope_or - unit test function of or
%
% Syntax:
%    res = test_zonotope_or
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
% Written:       26-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% check empty zonotope
Z = zonotope.empty(2);
Zres = Z | Z;
resvec(end+1) = true;

% instantiate zonotope
zono1 = zonotope([4 2 2;1 2 0]);
zono2 = zonotope([3 1 -1 1;3 1 2 0]);
Zres = zono1 | zono2;
resvec(end+1) = true;

% axis aligned zonotopes
zono1 = zonotope(interval([2;3],[3;4]));
zono2 = zonotope(interval([6;3],[7;4]));
Zres = zono1 | zono2;
resvec(end+1) = true;

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
