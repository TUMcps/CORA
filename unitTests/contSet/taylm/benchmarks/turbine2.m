function res = turbine2(v, w, r)
% turbine2 - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/rosa/Turbine.scala
%
% Syntax:  
%    res = turbine2(v, w, r)
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
    res = 6*v - 0.5 * v * (w*w*r*r) / (1-v) - 2.5;
end
%------------- END OF CODE --------------