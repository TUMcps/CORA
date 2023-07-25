function res = test_mptPolytope_not
% test_mptPolytope_not - unit test function of complement operation
%
% Syntax:  
%    res = test_mptPolytope_not
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
% Written:      25-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% assume true
res = true;

% init polytope
P = mptPolytope([1 0],1);
% compute complement
P_ = ~P;

% true solution
P_true = mptPolytope([-1 0],-1);

% check
res(end+1,1) = P_ == P_true;


% combine results
res = all(res);

%------------- END OF CODE --------------