function res = test_matPolytope_matZonotope()
% test_matPolytope_matZonotope - unit test function for conversion to 
%    matZonotope
% 
% Syntax:
%    res = test_matPolytope_matZonotope
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

% Authors:       Lukas Schäfer
% Written:       13-April-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init matPolytope
A_C = [1.0, 0.1; 0.1, 1.0];
B_C = [0.09; 0.09];
B_V = [0.01, 0.0; 0.0, -0.04];
v = [4, 4; 4, -4; -4, 4; -4, -4]';
V = cat(3,[A_C, B_C+B_V*v(:,1)], ...
    [A_C, B_C+B_V*v(:,2)], ...
    [A_C, B_C+B_V*v(:,3)], ...
    [A_C, B_C+B_V*v(:,4)]);
matP = matPolytope(V);

% enclose by matZonotope
matZ = matZonotope(matP);

% check dimensions
assert(all(size(matZ.C)==size(matP.V(:,:,1))));

% check if vertices are contained
assert(all(contains(matZ,V)));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
