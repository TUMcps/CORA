function res = eq(hs1,hs2)
% eq - overloaded '==' operator for exact comparison of two halfspaces
%
% Syntax:  
%    res = eq(hs1,hs2)
%
% Inputs:
%    hs1 - halfspace object
%    hs2 - halfspace object
%
% Outputs:
%    res - true/false
%
% Example: 
%    hs1 = halfspace([2;4], 4);
%    hs2 = halfspace([-3;5], 3);
%    hs1 == hs2
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

res = isequal(hs1,hs2);

%------------- END OF CODE --------------