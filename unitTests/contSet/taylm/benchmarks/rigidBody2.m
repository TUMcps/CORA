function res = rigidBody2(x1, x2, x3)
% rigidBody2 - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/control/RigidBody.scala
%
% Syntax:  
%    res = rigidBody2(x1, x2, x3)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:       Dmitry Grebenyuk
% Written:      10-October-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    res = 2*x1*x2*x3 + 3*x3*x3 - x2*x1*x2*x3 + 3*x3*x3 - x2;
end
%------------- END OF CODE --------------