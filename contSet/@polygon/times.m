function pgon = times(pgon, M)
% times - overloads the '.*' operator
%
% Syntax:
%    pgon = M .* pgon
%    pgon = pgon .* M
%    pgon = times(pgon, M)
%
% Inputs:
%    pgon - polygon
%    M - numeric
%
% Outputs:
%    pgon - polygon
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
if ~isnumeric(M) && ~(all(size(M) == [2, 1]) || isscalar(M))
    throw(CORAerror('CORA:notSupported', ...
        'Operation "times" is only defined for vectors matrices of dimension 2 or scalars!'));
end

% use mtimes
pgon = diag(M) * pgon;

end

% ------------------------------ END OF CODE ------------------------------
