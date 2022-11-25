function res = testLongDuration_interval_dim
% testLongDuration_interval_dim - unit test function of dim
%
% Syntax:  
%    res = testLongDuration_interval_dim
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      12-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% 1. Empty case
I = interval();

% compute dimension
dimI = dim(I);
true_dim = 0;
res_empty = dimI == true_dim;


% 2. Random cases
res_rand = true;
nrOfTests = 1000;

for i=1:nrOfTests
    
    % random dimension
    n = randi(50);
    
    % init random interval
    lb = -rand(n,1);
    ub = rand(n,1);
    I = interval(lb,ub);

    % compute dimension
    dimI = dim(I);

    % check with correct solution
    if dimI ~= n
        res_rand = false; break;
    end

end

% combine results
res = res_empty && res_rand;

if res
    disp('test_dim successful');
else
    disp('test_dim failed');
end

%------------- END OF CODE --------------
