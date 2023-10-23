function res = bspline3(u)
% bspline3 - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/rosa/Bsplines.scala
%
% Syntax:
%    res = bspline3(u)
%
% Inputs:
%    u - see benchmark
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

res = -u*u*u / 6.0;
end

% ------------------------------ END OF CODE ------------------------------
