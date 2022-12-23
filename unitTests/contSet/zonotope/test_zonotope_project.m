function res = test_zonotope_project
% test_zonotope_project - unit test function of project
%
% Syntax:  
%    res = test_zonotope_project
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      26-July-2016
% Last update:  09-August-2020 (MW, enhance randomness)
% Last revision:---

%------------- BEGIN CODE --------------

% create zonotope
Z = zonotope([-4, -3, -2, -1; 1, 2, 3, 4; 5, 5, 5, 5]);

% project zonotope
Z1 = project(Z,[1 3]);
c1 = center(Z1);
G1 = generators(Z1);

% logical indexing
Z2 = project(Z,[true false true]);
c2 = center(Z2);
G2 = generators(Z2);

% true result
true_c = [-4; 5];
true_G = [-3, -2, -1; ...
           5, 5, 5];

% check result
res = compareMatrices(c1,true_c) && compareMatrices(G1,true_G) ...
    && compareMatrices(c2,true_c) && compareMatrices(G2,true_G);

%------------- END OF CODE --------------
