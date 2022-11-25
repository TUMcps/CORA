function res = turbine3(v, w, r)
% turbine3 - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/rosa/Turbine.scala
%
% Syntax:  
%    res = turbine3(v, w, r)
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
    res = 3 - 2/(r*r) - 0.125 * (1+2*v) * (w*w*r*r) / (1-v) - 0.5;
end
%------------- END OF CODE --------------