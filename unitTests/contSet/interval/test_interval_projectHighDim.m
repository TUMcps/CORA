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
% See also: -

% Authors:       Tobias Ladner
% Written:       13-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% check empty
I = interval();
R = projectHighDim(I,2,[]);
Rtrue = interval([0;0],[0;0]);
resvec(end+1) = isequal(R,Rtrue);

% check default case
I = interval([2;3],[4;5]);
R = projectHighDim(I,5,[1,4]);
Rtrue = interval([2;0;0;3;0],[4;0;0;5;0]);
resvec(end+1) = isequal(R,Rtrue);

% gather results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
