function res = testLong_interval_eq
% testLong_interval_eq - unit test function of overloaded '==' operator
%
% Syntax:
%    res = testLong_interval_eq
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       04-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% number of tests
nrTests = 1000;

for i=1:nrTests

    % create interval
    n = floor(1 + 9*rand(1));
    lb1 = -10*rand(n,1);
    ub1 = 10*rand(n,1);
    I1 = interval(lb1, ub1);
    lb2 = -10*rand(n,1);
    ub2 = 10*rand(n,1);
    I2 = interval(lb2, ub2);
    
    % compute non-equality
    if ~(I1 == I1)
        res = false; break
    elseif I1 == I2
        res = false; break
    end

end

% ------------------------------ END OF CODE ------------------------------
