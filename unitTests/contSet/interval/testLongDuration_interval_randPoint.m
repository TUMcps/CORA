function res = testLongDuration_interval_randPoint
% testLongDuration_interval_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = testLongDuration_interval_randPoint
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: @interval > containsPoint.m
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  23-March-2021 (MW, more thorough test, empty case)
% Last revision:---

%------------- BEGIN CODE --------------

% empty case
Int = interval();
res_empty = true;
try
    p = randPoint(Int); % <- should throw error here
    res_empty = false;
end

% random intervals
nrTests = 1000;
nrRandPoints = 1;
res_rand = true;
for i=1:nrTests
    
    % random dimension
    n = randi(25);
    % random bounds
    lb = -rand(n,1);
    ub = rand(n,1);
    
    % init interval
    Int = interval(lb,ub);
    
    % check 'standard' method
    p = randPoint(Int,nrRandPoints);
    % check 'extreme' method
    p = [p, randPoint(Int,nrRandPoints,'extreme')];
    % check 'gaussian' method (only for run-time errors)
    pr = 0.8;
    randPoint(Int,nrRandPoints,'gaussian',pr);
    
    if ~in(Int,p)
        res_rand = false; break;
    end
    
end


% combine results
res = res_rand && res_empty;

if res
    disp('test_randPoint successful');
else
    disp('test_randPoint failed');
end

%------------- END OF CODE --------------