function res = test_interval_isequal
% test_interval_isequal - unit test function of isequal
%
% Syntax:  
%    res = test_interval_isequal
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

% Author:       Mark Wetzlinger
% Written:      17-September-2019
% Last update:  03-December-2022 (MW, add Inf case)
%               23-December-2022 (MW, add matrix case)
% Last revision:---

%------------- BEGIN CODE --------------

% assume true
res = true;

% instantiate interval
lb = [-3; -9; -4; -7; -1];
ub = [4;   2;  6;  3;  8];
I1 = interval(lb, ub);
ub = [4;   2;  6;  2;  8];
I2 = interval(lb, ub);

% check equality
if ~isequal(I1,I1) || isequal(I1,I2)
    res = false;
end

% interval with Inf values
lb = [-Inf; -9; -4; -7; -1];
ub = [4;   2;  6;  Inf;  8];
I1 = interval(lb, ub);
ub = [4;   2;  Inf;  2;  8];
I2 = interval(lb, ub);

if ~isequal(I1,I1) || isequal(I1,I2)
    res = false;
end

% instantiate matrix interval
lb = [-2 -3   -Inf; -4 -Inf -1];
ub = [ 2  Inf 5;    Inf 3    0];
I1 = interval(lb,ub);
ub = [ 2  Inf Inf;  Inf 3    0];
I2 = interval(lb,ub);
I3 = I1([1,2],[1,2]);

if ~isequal(I1,I1) || isequal(I1,I2) || isequal(I1,I3)
    res = false;
end

%------------- END OF CODE --------------