function res = doppler(u, v, T)
% doppler - a benchmark from https://github.com/malyzajko/daisy/blob/master/testcases/rosa/Doppler.scala
%
% Syntax:
%    res = doppler(u, v, T)
%
% Inputs:
%    u, v, T - see benchmark
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

t1 = 331.4 + 0.6 * T;
res = (- (t1) *v) / ((t1 + u)*(t1 + u));
end

% ------------------------------ END OF CODE ------------------------------
