function res = test_interval_dim
% test_interval_dim - unit test function of dim
%
% Syntax:
%    res = test_interval_dim
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
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. Empty case
I = interval();

% compute dimension
dimI = dim(I);
true_dim = 0;
res(1) = dimI == true_dim;

    
% init random interval
lb = [-3; -2; -5];
ub = [4; 2; 1];
I = interval(lb,ub);

% compute dimension
res(2) = dim(I) == 3;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
