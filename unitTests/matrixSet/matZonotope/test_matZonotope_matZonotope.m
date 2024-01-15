function res = test_matZonotope_matZonotope
% test_matZonotope_matZonotope - unit test function for constructor
% 
% Syntax:
%    res = test_matZonotope_matZonotope
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

% Authors:       Mark Wetzlinger
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% only center
C = 0;
matZ = matZonotope(C);
res(end+1,1) = all(withinTol(matZ.center,C));
C = [1; 1; 0];
matZ = matZonotope(C);
res(end+1,1) = all(withinTol(matZ.center,C));

% center and one generator
G{1} = [1; 0; 0];
matZ = matZonotope(C,G);
res(end+1,1) = all(withinTol(matZ.center,C)) && ...
    length(matZ.generator) == 1 && all(withinTol(matZ.generator{1},G{1}));

% center and multiple generators
C = [1 2 1; 3 2 0];
G{1} = [2 0 1; -1 1 -2];
G{2} = [3 1 0; -1 -1 4];
G{3} = [0 1 -1; 3 1 2];
matZ = matZonotope(C,G);
res(end+1,1) = all(all(withinTol(matZ.center,C))) && ...
    length(matZ.generator) == 3 && all(all(withinTol(matZ.generator{1},G{1}))) ...
    && all(all(withinTol(matZ.generator{2},G{2}))) ...
    && all(all(withinTol(matZ.generator{3},G{3})));

% copy constructor
matZ_ = matZonotope(matZ);

% conversion from zonotope: empty, only center, one/multiple generator(s)
Z = zonotope.empty(2);
matZ = matZonotope(Z);
Z = zonotope([1;2;1]);
matZ = matZonotope(Z);
Z = zonotope([1; 2],[1; 0]);
matZ = matZonotope(Z);
Z = zonotope([1; 2],[1 2 4 2; 0 2 -3 1]);
matZ = matZonotope(Z);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
