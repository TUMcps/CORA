function res = testLongDuration_zonotope_isFullDim
% testLongDuration_zonotope_isFullDim - unit test function of isFullDim
%
% Syntax:  
%    res = testLongDuration_zonotope_isFullDim
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
res_empty = true;

Z = zonotope();
if isFullDim(Z)
    res_empty = false;
end


% random tests
res_rand = true;
nrOfTests = 1000;

for i=1:nrOfTests
    % random dimension
    n = randi([2,50]);
    % random center
    c = randn(n,1);
    % random number of generators
    gamma = randi([n,100]);
    % random full-rank generator matrix
    G = randn(n,gamma);
    while rank(G) < n
        G = randn(n,gamma);
    end
    
    % instantiate zonotope
    Z = zonotope(c,G);
    
    % check if full-dimensional
    if ~isFullDim(Z)
        res_rand = false; break
    end
end


% combine results
res = res_empty && res_rand;

if res
    disp('test_isFullDim successful');
else
    disp('test_isFullDim failed');
end

%------------- END OF CODE --------------