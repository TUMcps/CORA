function P_out = not(P)
% not - overloads '~' operator to compute the complement of a polytope
%    note: only supported for single halfspaces and inequalities; also, the
%    given implement does not support Ax < b (instead of <=), hence, we
%    dirtily shift the offset vector b by eps
%
% Syntax:
%    P_out = ~P;
%    P_out = not(P);
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    P_out - set complement of polytope object
%
% Example: 
%    P = polytope([1 0],1);
%    P_out = ~P;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       25-April-2023
% Last update:   16-December-2023 (MW, support empty sets)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% different cases
if isempty(P.Ae) && length(P.b) == 1
    % exactly one inequality constraint -> switch inequality sign:
    %    Ax <= b  ->  Ax > b <=> -Ax < -b
    % shift by eps to account for < instead of <= (which we have to use)
    P_out = polytope(-P.A,-P.b-eps);

elseif representsa_(P,'emptySet',eps)
    % empty set -> complement is R^n
    P_out = polytope.Inf(dim(P));

elseif representsa_(P,'fullspace',eps)
    % R^n -> complement is the empty set
    P_out = polytope.empty(dim(P));

else
    % other cases not supported because the complement cannot be expressed
    % as a polytope
    throw(CORAerror('CORA:notSupported',['Complement operation is only '...
        'supported for a single inequality constraint or empty sets.']));
end

% ------------------------------ END OF CODE ------------------------------
