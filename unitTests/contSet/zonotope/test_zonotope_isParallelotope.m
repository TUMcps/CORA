function res = test_zonotope_isParallelotope
% test_zonotope_isParallelotope - unit test function of isParallelotope
%
% Syntax:  
%    res = test_zonotope_isParallelotope
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

% Author:       Mark Wetzlinger
% Written:      08-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check empty zonotope
Z = zonotope();
res(1) = ~isParallelotope(Z);

% instantiate parallelotope
c = [-2; 1];
G = [2 4; -2 3];
Z = zonotope(c,G);

% check result
res(2) = isParallelotope(Z);

% add zero-length generators
G = [G, zeros(2,2)];
Z = zonotope(c,G);

% still a parallelotope
res(3) = isParallelotope(Z);

% add generator
G = [G, [4; -2]];
Z = zonotope(c,G);

% not a parallelotope anymore
res(4) = ~isParallelotope(Z);

% no generator matrix
Z = zonotope(c);
res(5) = ~isParallelotope(Z);

% currently, no check for aligned generators...

% combine results
res = all(res);

%------------- END OF CODE --------------