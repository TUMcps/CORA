function pgon = mtimes(pgon,M)
% mtimes - overloads the '*' operator
%
% Syntax:
%    pgon = M * pgon
%    pgon = pgon * M
%    pgon = mtimes(pgon, M)
%
% Inputs:
%    pgon - polygon
%    M - numeric
%
% Outputs:
%    han - plot handle
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% find class argument
[pgon, M] = findClassArg(pgon, M, 'polygon');

% check dimension of the matrix
if ~isnumeric(M) || (any(size(M) ~= [2, 2]) && ~isscalar(M))
    throw(CORAerror('CORA:notDefined', ...
        'Operation "polygon/mtimes" is only defined for square matrices of dimension 2 or scalars!'));
end

% multiplication with matrix
w = warning();
warning('off');

V = vertices_(pgon);
V = M * V;

pgon = polygon(V);

warning(w);

end

% ------------------------------ END OF CODE ------------------------------
