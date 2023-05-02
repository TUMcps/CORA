function res = test_mptPolytope_uplus
% test_mptPolytope_uplus - unit test function of uplus
%
% Syntax:  
%    res = test_mptPolytope_uplus
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

% Author:       Tobias Ladner
% Written:      06-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

resvec = true(0);

% init
C = [1 0 -1 0 1; 0 1 0 -1 1]';
d = [3; 2; 3; 2; 1];
P = mptPolytope(C,d);

% plus
pP = +P;
resvec(end+1) = all(pP.P.A == C, 'all');
resvec(end+1) = all(pP.P.b == d, 'all');

% compare with P
resvec(end+1) = isequal(pP, P);

% test empty case
resvec(end+1) = isempty(+mptPolytope());

% add results
res = all(resvec);

%------------- END OF CODE --------------