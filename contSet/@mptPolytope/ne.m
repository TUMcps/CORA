function res = ne(P1,P2,varargin)
% ne - overloads the '~=' operator for exact comparison of two polytopes;
%    tolerance is not enforced due to call to mpt toolbox
%
% Syntax:
%    res = P1 ~= P2
%    res = ne(P1,P2)
%    res = ne(P1,P2,tol)
%
% Inputs:
%    P1 - mptPolytope object
%    P2 - mptPolytope object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      31-August-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = P1.P ~= P2.P;

%------------- END OF CODE --------------