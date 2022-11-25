function res = turbine1(v, w, r)
% turbine1 - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/rosa/Turbine.scala
% v in [-4.5,-0.3] && w in [0.4, 0.9], r in [3.8, 7.8]
% Syntax:  
%    res = turbine1(v, w, r)
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

    res = 3 + 2/(r*r) - 0.125*(3-2*v)*(w*w*r*r)/(1-v) - 4.5;

end
%------------- END OF CODE --------------