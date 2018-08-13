function res = bspline2(u)
% bspline2 - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/rosa/Bsplines.scala
%
% Syntax:  
%    res = bspline2(u)
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
    res = (-3 * u*u*u  + 3*u*u + 3*u + 1) / 6.0;
end
%------------- END OF CODE --------------