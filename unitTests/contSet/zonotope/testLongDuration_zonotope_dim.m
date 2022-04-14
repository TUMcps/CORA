function res = testLongDuration_zonotope_dim
% testLongDuration_zonotope_dim - unit test function of dim
%
% Syntax:  
%    res = testLongDuration_zonotope_dim
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
% Written:      11-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% number of tests
nrOfTests = 1000;

for i=1:nrOfTests
    % random dimension
    n = randi([2,50]);
    % random center
    c = randn(n,1);
    % random generator matrix (also empty)
    if rand(1) < 0.05
        G = [];
    else
        G = 5*randn(n,randi(10));
    end
    
    % instantiate zonotope
    Z = zonotope(c,G);
    
    % get dimension
    Zdim = dim(Z);
    
    % assert correctness
    if Zdim ~= n
        res = false; break
    end
end


if res
    disp('test_dim successful');
else
    disp('test_dim failed');
end

%------------- END OF CODE --------------