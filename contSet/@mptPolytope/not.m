function P = not(P)
% not - overloads '~' operator to compute the complement of a polytope
%    note: only supported for single halfspaces and inequalities; also, the
%    given implement does not support Ax < b (instead of <=), hence, we
%    dirtily shift the offset vector b by eps
%
% Syntax:
%    P = ~P;
%    P = not(P);
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    P - complement of mptPolytope object
%
% Example: 
%    P = polytope([1 0],1);
%    P_ = ~P;
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

% only supported for single halfspaces and no equality constraints
if ~isempty(P.P.Ae) || length(P.P.b) > 1
    throw(CORAerror('CORA:notSupported',['Complement operation is only '...
        'supported for a single inequality constraint.']));
end

% switch inequality sign: complement of Ax <= b is Ax > b <=> -Ax < -b
% to account for < instead of <= (which we have to use), shift by eps
P = mptPolytope(-P.P.A,-P.P.b-eps);

%------------- END OF CODE --------------