function res = test_interval_isZero
% test_interval_isZero - unit test function of isZero
%
% Syntax:  
%    res = test_interval_isZero
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

% Author:       Mark Wetzlinger
% Written:      17-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty case
res(1) = ~isZero(interval());

% only origin
I = interval(zeros(3,1),zeros(3,1));
res(2) = isZero(I);

% shifted center
I = interval([0.01;0.02],[0.03;0.025]);
res(3) = ~isZero(I);

% shifted center, contains origin within tolerance
I = interval([0.01;-0.01],[0.02;0.01]);
tol = 0.05;
res(4) = isZero(I,tol);
  
% combine results
res = all(res);

%------------- END OF CODE --------------
