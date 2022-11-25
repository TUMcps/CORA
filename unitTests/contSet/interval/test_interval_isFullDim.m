function res = test_interval_isFullDim
% test_interval_isFullDim - unit test function of isFullDim
%
% Syntax:  
%    res = test_interval_isFullDim
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
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% 1. Empty case
I = interval();

% compute dimension
res(1) = ~isFullDim(I);

% init full-dimensional interval
lb = [-2; -4; -7; -1; -2];
ub = [4; 2; 6; 4; 8];
I = interval(lb,ub);

% check with correct solution
res(2) = isFullDim(I);

% init random lower-dimensional interval
% ... by setting random dimension to 0
n0 = 3;
lb(n0) = 0;
ub(n0) = 0;
I = interval(lb,ub);

% check with correct solution
res(3) = ~isFullDim(I);

% combine results
res = all(res);

if res
    disp('test_isFullDim successful');
else
    disp('test_isFullDim failed');
end

%------------- END OF CODE --------------
