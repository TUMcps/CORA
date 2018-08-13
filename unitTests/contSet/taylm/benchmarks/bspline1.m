function res = bspline1(u)
% bspline1 - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/rosa/Bsplines.scala
%
% Syntax:  
%    res = bspline1(u)
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
    res = (3 * u*u*u - 6 * u*u + 4) / 6.0;
end
%------------- END OF CODE --------------