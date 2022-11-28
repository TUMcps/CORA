function res = eq(cPZ1,cPZ2)
% eq - overloaded '==' operator for exact comparison of two constrained polynomial zonotopes
%
% Syntax:  
%    res = eq(cPZ1,cPZ2)
%
% Inputs:
%    cPZ1 - conPolyZono object
%    cPZ2 - conPolyZono object
%
% Outputs:
%    res - true/false
%
% Example: 
%    cPZ1 = conPolyZono([0;0],[1 0 1;0 -1 1],[1 0 2;0 1 1],[1 -0.5], ...
%                       -1,[0 2;1 0],[0.4 0;0.1 1]);
%    cPZ2 = conPolyZono([0;0],[1 1 0;1 0 -1],[2 1 0;1 0 1],[0.5 -1], ...
%                       1,[2 0;0 1],[0 0.4;1 0.1]);
%    cPZ1 == cPZ2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mingrui Wang
% Written:      21-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = isequal(cPZ1,cPZ2);

%------------- END OF CODE --------------