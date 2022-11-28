function res = test_mptPolytope_isemptyobject
% test_mptPolytope_isemptyobject - unit test function of isemptyobject
%
% Syntax:  
%    res = test_mptPolytope_isemptyobject
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      03-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate polytopes
P1 = mptPolytope();

A = [-1 0; 2 4; 1 -2];
b = [-1; 14; -1];

P2 = mptPolytope(A,b);

% check results
res = isemptyobject(P1) && ~isemptyobject(P2);

%------------- END OF CODE --------------
