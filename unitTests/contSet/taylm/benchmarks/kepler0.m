function res = kepler0(x1, x2, x3, x4, x5, x6)
% kepler0 - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/real2float/Kepler.scala
%
% Syntax:
%    res = kepler0(x1, x2, x3, x4, x5, x6)
%
% Inputs:
%    x1-x6 - see benchmark
%
% Outputs:
%    res - result 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Authors:       Dmitry Grebenyuk
% Written:       10-October-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = x2 * x5 + x3 * x6 - x2 * x3 - x5 * x6 + x1 * (-x1 + x2 + x3 - x4 + x5 + x6);
end

% ------------------------------ END OF CODE ------------------------------
