function res = himmilbeau(x1, x2)
% himmilbeau - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/real2float/Himmilbeau.scala
%
% Syntax:
%    res = himmilbeau(x1, x2)
%
% Inputs:
%    x1-x2 - see benchmark
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

res = (x1*x1 + x2 - 11)*(x1 * x1 + x2 - 11) + (x1 + x2*x2 - 7)*(x1 + x2*x2 - 7);
end

% ------------------------------ END OF CODE ------------------------------
