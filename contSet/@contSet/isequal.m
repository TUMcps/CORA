function res = isequal(varargin)
% isequal - checks if two sets are equal
%
% Syntax:
%    res = isequal(S1, S2)
%    res = isequal(S1, S2, tol)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",varargin{:}))

% ------------------------------ END OF CODE ------------------------------
