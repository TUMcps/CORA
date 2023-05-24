function res = eq(P1,P2,varargin)
% eq - overloads the '==' operator for exact comparison of two polytopes;
%    tolerance is not enforced due to call to mpt toolbox
%
% Syntax:
%    res = P1 == P2
%    res = eq(P1,P2)
%    res = eq(P1,P2,tol)
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
% See also: none

% Author:       Matthias Althoff
% Written:      31-August-2015
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if isa(P2,'mptPolytope')
    res = P1.P == P2.P;
elseif isa(P2,'interval')
    P2 = mptPolytope(P2);
    res = P1.P == P2.P;
else
    throw(CORAerror('CORA:noops',P1,P2));
end

%------------- END OF CODE --------------