function res = test_interval_reduce
% test_interval_reduce - unit test function of reduce
%
% Syntax:
%    res = test_interval_reduce
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
I = interval.empty(2);
% reduce
I_ = reduce(I);
% should remain the same...
res = isequal(I,I_);

% init column interval
I = interval([-2;-1],[1;3]);
% reduce
I_ = reduce(I);
% check result
res(end+1,1) = isequal(I,I_);

% init matrix interval
I = interval([-2 -1; -4 -2],[1 3; 5 2]);
% reduce
I_ = reduce(I);
% check result
res(end+1,1) = isequal(I,I_);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
