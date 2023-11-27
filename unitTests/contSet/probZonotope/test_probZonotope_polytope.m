function res = test_probZonotope_polytope
% test_probZonotope_polytope - unit test to check whether the mean of a
%    probabilistic zonotope can be represented by a polytope
%
% Syntax:
%    res = test_probZonotope_polytope
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

% Authors:       Adrian Kulmburg
% Written:       16-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate probabilistic zonotope
Z1 = [0 1 0; 0 0 1];
Z2 = [0.6 1.2  ; 0.6 -1.2 ];
pZ3 = probZonotope(Z1,Z2,2);

% convert it to a polytope

P = polytope(pZ3);

res = 1;

% ------------------------------ END OF CODE ------------------------------
