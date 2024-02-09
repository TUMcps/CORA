function simRes = times(factor1,factor2)
% times - Overloaded '.*' operator for the states of the simulated trajectories
%
% Syntax:
%    simRes = times(A,simRes)
%
% Inputs:
%    A - numeric vector
%    simRes - simResult object
%
% Outputs:
%    simRes - transformed simResult object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       02-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[simRes,A] = findClassArg(factor1, factor2, 'simResult');

% parse input
if ~isnumeric(A)
    throw(CORAerror('CORA:noops', simRes, A))
elseif ~isvector(A) && ~isempty(A)
    throw(CORAerror('CORA:notSupported', 'Multiplied vector has to be a column vector.'))
end

% transform column vector to diagonal matrix and use mtimes
A = diag(A);
simRes = A * simRes;

end

% ------------------------------ END OF CODE ------------------------------
