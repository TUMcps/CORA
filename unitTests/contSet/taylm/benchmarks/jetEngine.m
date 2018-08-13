function res = jetEngine(x1, x2)
% jetEngine - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/control/JetEngine.scala
%
% Syntax:  
%    res = jetEngine(x1, x2) 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author:       Dmitry Grebenyuk
% Written:      09-October-2017
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    res = x1 + (...
      (2*x1*((3*x1*x1 + 2*x2 - x1)/(x1*x1 + 1))*((3*x1*x1 + 2*x2 - x1)/(x1*x1 + 1) - 3) +...
     x1*x1*(4*((3*x1*x1 + 2*x2 - x1)/(x1*x1 + 1))-6))...
      *(x1*x1 + 1) +...
    3*x1*x1*((3*x1*x1 + 2*x2 - x1)/(x1*x1 + 1)) + x1*x1*x1 + x1 + 3*((3*x1*x1 + 2*x2 -x1)/(x1*x1 + 1)));
end
%------------- END OF CODE --------------