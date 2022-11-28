function res = eq(Z1,Z2)
% eq - overloaded '==' operator for exact comparison of two zonotopes
%
% Syntax:  
%    res = eq(Z1,Z2)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
%
% Outputs:
%    res - true/false
%
% Example: 
%    Z1 = zonotope(zeros(3,1),rand(3,5));
%    Z2 = zonotope(zeros(3,1),rand(3,5));
%    Z1 == Z2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      20-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = isequal(Z1,Z2);

%------------- END OF CODE --------------
