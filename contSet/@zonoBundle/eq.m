function res = eq(zB1,zB2)
% eq - overloaded '==' operator for exact comparison of two zonotope bundles
%
% Syntax:  
%    res = eq(zB1,zB2)
%
% Inputs:
%    zB1 - zonoBundle object
%    zB2 - zonoBundle object
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval([0;0],[2;2]);
%    I2 = interval([0;1],[2;3]);
%    I3 = interval([0;2],[2;4]);
%
%    zB1 = zonoBundle({zonotope(I1),zonotope(I2)});
%    zB2 = zonoBundle({zonotope(I1),zonotope(I3)});
%
%    zB1 = zB2
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

res = isequal(zB1,zB2);

%------------- END OF CODE --------------