function ls = levelSet(P)
% levelSet - converts polytope object to an equivalent levelSet object
%
% Syntax:
%    ls = levelSet(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    ls - levelSet object
%
% Example:
%    P = polytope([1 0],1);
%    ls = levelSet(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl
% Written:       23-March-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute halfspace representation if only vertices given
constraints(P);

if ~isempty(P.Ae_.val)
    throw(CORAerror('CORA:notSupported','Equality constraints not supported.'));
end

% read out variables
vars = sym('x',[dim(P), 1]);
A = P.A_.val; b = P.b_.val;

% init symbolic equations
eq = A*vars - b;
% stack <= comparison operator times the number of inequalities
compOps = repmat({'<='},length(P.b_.val),1);

% init resulting level set
ls = levelSet(eq,vars,compOps);

% ------------------------------ END OF CODE ------------------------------
