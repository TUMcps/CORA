function res = test_zonotope_vertices
% test_zonotope_vertices - unit test function of vertices
%
% Syntax:  
%    res = test_zonotope_vertices
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

% Author:       Matthias Althoff
% Written:      26-July-2016
% Last update:  09-February-2017
% Last revision:---

%------------- BEGIN CODE --------------

% empty set
Z_e = zonotope();
res_e = isempty(vertices(Z_e)) && isnumeric(vertices(Z_e));

% create zonotope
Z1 = zonotope([-4, -3, -2, -1; 1, 2, 3, 4]);

% obtain result
Vmat = vertices(Z1);

% true result
true_Vmat = [-8, 0, 2, -4, -4, -10; ...
              2, 0, -8, -4, 6, 10];

% res
res = compareMatrices(Vmat,true_Vmat) && res_e;

%------------- END OF CODE --------------
