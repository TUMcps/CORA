function res = test_mptPolytope_isZero
% test_mptPolytope_isZero - unit test function of isZero
%
% Syntax:  
%    res = test_mptPolytope_isZero
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
% Written:      17-March-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty case
res(1) = ~isZero(mptPolytope());

% only origin
P = mptPolytope([1 0; 0 1; -1 0; 0 -1],zeros(4,1));
res(end+1) = isZero(P);

% shifted center
P = P + [0.01; 0];
res(end+1) = ~isZero(P);
% add tolerance
tol = 0.02;
res(end+1) = isZero(P,tol);

% combine results
res = all(res);

%------------- END OF CODE --------------
