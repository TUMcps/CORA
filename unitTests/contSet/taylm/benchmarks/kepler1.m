function res = kepler1(x1, x2, x3, x4)
% kepler1 - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/real2float/Kepler.scala
%
% Syntax:  
%    res = kepler1(x1, x2, x3, x4)
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
    res = x1 * x4 * (-x1 + x2 + x3 - x4) + x2 * (x1 - x2 + x3 + x4) + x3 * (x1 + x2 - x3 + x4) - ...
      x2 * x3 * x4 - x1 * x3 - x1 * x2 - x4;
end
%------------- END OF CODE --------------