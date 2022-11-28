function res = eq(hyp1,hyp2)
% eq - overloaded '==' operator for exact comparison of two hyperplanes
%
% Syntax:  
%    res = eq(hyp1,hyp2)
%
% Inputs:
%    hyp1 - hyperplane object
%    hyp2 - hyperplane object
%
% Outputs:
%    res - true/false
%
% Example: 
%    hyp1 = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2]);
%    hyp2 = conHyperplane(halfspace([1;-1],0),[1 0;-1 0],[2;2]);
%    hyp1 == hyp2
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

res = isequal(hyp1,hyp2);

%------------- END OF CODE --------------