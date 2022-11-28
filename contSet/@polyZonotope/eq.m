function res = eq(pZ1,pZ2)
% eq - overloaded '==' operator for exact comparison of two polynomial zonotopes
%
% Syntax:  
%    res = eq(pZ1,pZ2)
%
% Inputs:
%    pZ1 - polyZonotope object
%    pZ2 - polyZonotope object
%
% Outputs:
%    res - true/false
%
% Example: 
%    pZ1 = polyZonotope([0;0],[1 0 1;0 -1 1],[0.4 0;0.1 1],[1 0 2;0 1 1]);
%    pZ2 = polyZonotope([0;0],[1 1 0;1 0 -1],[0 0.4;1 0.1],[2 1 0;1 0 1]);
%    pZ1 == pZ2
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

res = isequal(pZ1,pZ2);

%------------- END OF CODE --------------