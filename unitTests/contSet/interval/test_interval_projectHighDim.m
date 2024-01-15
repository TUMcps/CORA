function res = test_interval_projectHighDim
% test_interval_projectHighDim - unit test function of projectHighDim
%
% Syntax:
%    res = test_interval_projectHighDim
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

% Authors:       Tobias Ladner
% Written:       13-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% check default case
I = interval([2;3],[4;5]);
R = projectHighDim(I,5,[1,4]);
Rtrue = interval([2;0;0;3;0],[4;0;0;5;0]);
resvec(end+1) = isequal(R,Rtrue);

% gather results
res = all(resvec);

% empty cannot be projected
% check empty
I = interval.empty(1);
try
    R = projectHighDim(I,2,[]);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
