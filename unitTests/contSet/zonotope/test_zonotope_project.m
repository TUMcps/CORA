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

% 1. Analytical Test ------------------------------------------------------

% create zonotope
Z = zonotope([-4, -3, -2, -1; 1, 2, 3, 4; 5, 5, 5, 5]);

% project zonotope
Zres = project(Z,[1 3]);
Zmat1 = Zres.Z;

% logical indexing
Zres = project(Z,[true false true]);
Zmat2 = Zres.Z;

% true result
true_mat = [-4, -3, -2, -1; ...
            5, 5, 5, 5];

% check result
res = all(all(Zmat1 == true_mat)) && all(all(Zmat2 == true_mat));

%------------- END OF CODE --------------
