function res = rigidBody1(x1, x2, x3)
% rigidBody1 - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/control/RigidBody.scala
%
% Syntax:  
%    res = rigidBody1(x1, x2, x3)
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
    res = -x1*x2 - 2*x2*x3 - x1 - x3;
end
%------------- END OF CODE --------------