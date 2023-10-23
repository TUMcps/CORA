function res = testLong_interval_randPoint
% testLong_interval_randPoint - unit test function of randPoint
%
% Syntax:
%    res = testLong_interval_randPoint
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: @interval > containsPoint.m
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   23-March-2021 (MW, more thorough test, empty case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% random intervals
nrTests = 100;
nrRandPoints = 1;
res = true;
for i=1:nrTests
    
    % random dimension
    n = randi(25);
    % random bounds
    lb = -rand(n,1);
    ub = rand(n,1);
    
    % init interval
    I = interval(lb,ub);
    
    % check 'standard' method
    p = randPoint(I,nrRandPoints);
    % check 'extreme' method
    p = [p, randPoint(I,nrRandPoints,'extreme')];
    % check 'gaussian' method (only for run-time errors)
    pr = 0.8;
    randPoint(I,nrRandPoints,'gaussian',pr);
    
    if ~contains(I,p)
        res = false; break;
    end
    
end

% ------------------------------ END OF CODE ------------------------------
