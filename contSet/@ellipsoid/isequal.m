function res = isequal(E1,E2,varargin)
% isequal - checks if two ellipsoids are equal
%
% Syntax:
%    res = isequal(E1,E2)
%    res = isequal(E1,E2,tol)
%
% Inputs:
%    E1 - ellipsoid object
%    E2 - ellipsoid object
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
% See also: -

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   15-October-2019
%                19-March-2021 (use 'eq')
%                04-July-2022 (VG, class array case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% parse input arguments
tol = setDefaultValues({min(E1.TOL,E2.TOL)},varargin);

% check input arguments
inputArgsCheck({{E1,'att','ellipsoid','scalar'};
                {E2,'att','ellipsoid','scalar'};
                {tol,'att','numeric',{'nonnan','nonnegative','scalar'}}});

% assume false
res = false;

% check if dimensions are equal
if dim(E1) ~= dim(E2)
    return;
end

% check for emptyness
if (representsa_(E1,'emptySet',eps) && ~representsa_(E2,'emptySet',eps)) ...
        || (~representsa_(E1,'emptySet',eps) && representsa_(E2,'emptySet',eps))
    return;
elseif representsa_(E1,'emptySet',eps) && representsa_(E2,'emptySet',eps)
    res = true;
    return;
end

% compare shape matrix and center numerically
res = all(all(withinTol(E1.Q,E2.Q,tol))) && all(withinTol(E1.q,E2.q,tol));

% ------------------------------ END OF CODE ------------------------------
