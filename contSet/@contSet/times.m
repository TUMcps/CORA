function S = times(factor1,factor2)
% times - Overloaded '.*' operator for the multiplication with a contSet
%
% Syntax:
%    E = times(A,E)
%
% Inputs:
%    A - numerical column vector
%    S - contSet object 
%
% Outputs:
%    S - contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/mtimes

% Authors:       Tobias Ladner
% Written:       06-April-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[S,A] = findClassArg(factor1, factor2, 'contSet');

% parse input
if ~isnumeric(A)
    throw(CORAerror('CORA:noops', S, A))
elseif ~isvector(A) && ~isempty(A)
    throw(CORAerror('CORA:notSupported', 'Multiplied vector has to be a column vector.'))
end

% transform column vector to diagonal matrix and use mtimes
A = diag(A);
S = A * S;

% ------------------------------ END OF CODE ------------------------------
