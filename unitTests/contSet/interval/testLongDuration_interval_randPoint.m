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
I = interval();
res_empty = isempty(randPoint(I));

% random intervals
nrTests = 100;
nrRandPoints = 1;
res_rand = true;
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
        res_rand = false; break;
    end
    
end


% combine results
res = res_rand && res_empty;

if ~res
    path = pathFailedTests(mfilename());
    save(path,'lb','ub','n','p');
end

%------------- END OF CODE --------------