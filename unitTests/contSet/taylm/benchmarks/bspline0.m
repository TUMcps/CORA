function res = bspline0(u)
% bspline0 - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/rosa/Bsplines.scala
%
% Syntax:  
%    res = bspline0(u)
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
    res = (1 - u) * (1 - u) * (1 - u) / 6.0;
end
%------------- END OF CODE --------------