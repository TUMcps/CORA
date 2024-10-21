function res = isequal(E,S,varargin)
% isequal - checks if an ellipsoid is equal to another set or point
%
% Syntax:
%    res = isequal(E,S)
%    res = isequal(E,S,tol)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet object, numeric
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    E1 = ellipsoid([1,0;0,1/2],[1;1]);
%    E2 = ellipsoid([1+1e-15,0;0,1/2],[1;1]);
%    res = isequal(E1,E2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   15-October-2019
%                19-March-2021 (use 'eq')
%                04-July-2022 (VG, class array case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,3);

% parse input arguments
tol = setDefaultValues({min(E.TOL,S.TOL)},varargin);

% check input arguments
inputArgsCheck({{E,'att',{'ellipsoid','numeric'},'scalar'};
                {S,'att',{'contSet','numeric'},'scalar'};
                {tol,'att','numeric',{'nonnan','nonnegative','scalar'}}});

% ensure that numeric is second input argument
[E,S] = reorderNumeric(E,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < E.precedence
    res = isequal(S,E,tol);
    return
end

% ambient dimensions must match
if ~equalDimCheck(E,S,true)
    res = false;
    return
end

% ellipsoid-ellipsoid case
if isa(S,'ellipsoid')
    res = aux_isequal_ellipsoid(E,S,tol);
    return
end

throw(CORAerror('CORA:noops',E,S));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isequal_ellipsoid(E,S,tol)

% assume false
res = false;

% check if dimensions are equal
if dim(E) ~= dim(S)
    return;
end

% check for emptyness
if (representsa_(E,'emptySet',eps) && ~representsa_(S,'emptySet',eps)) ...
        || (~representsa_(E,'emptySet',eps) && representsa_(S,'emptySet',eps))
    return;
elseif representsa_(E,'emptySet',eps) && representsa_(S,'emptySet',eps)
    res = true;
    return;
end

% compare shape matrix and center numerically
res = all(all(withinTol(E.Q,S.Q,tol))) && all(withinTol(E.q,S.q,tol));

end

% ------------------------------ END OF CODE ------------------------------
